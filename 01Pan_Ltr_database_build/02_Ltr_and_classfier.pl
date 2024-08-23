#
# Module Dependence
#
use strict;
use FindBin;
use lib $FindBin::RealBin;
use Getopt::Long;
use POSIX qw(:sys_wait_h ceil floor);
use File::Copy;
use File::Spec;
use File::Path;
use File::Basename;
use Cwd qw(abs_path getcwd cwd);
use Data::Dumper;
use Pod::Text;
use Time::HiRes qw( gettimeofday tv_interval);

# RepeatModeler Libraries
use RepModelConfig;
use lib $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'};
use RepeatUtil;
use SeedAlignment;
use SeedAlignmentCollection;
use ThreadedTaskSimple;

# RepeatMasker Libraries
use SearchResult;
use SearchResultCollection;
use WUBlastSearchEngine;
use NCBIBlastSearchEngine;
use SeqDBI;
use SimpleBatcher;
use FastaDB;


# use Devel::Size qw(size total_size);

#
# Global Class Variables/Constants
#
my $CLASS = "RepeatModeler";
my $DEBUG = 0;
$DEBUG = 1 if ( $RepModelConfig::DEBUGALL == 1 );
$| = 1;    # Turn autoflush on
my %TimeBefore = ();
my @xx = glob 'RM*' ;
my $tmpDir = @xx[0];
my $LOG;
open $LOG, ">>", "$tmpDir/rmod.log";
#
# Version
#
my $version = $RepModelConfig::VERSION;

my $cmdLine = $0 . join( ' ', @ARGV );
my @getopt_args = (
                    '-version',                 # print out the version and exit
                    '-help',
                    '-dir=s',
                    '-debug',
                    '-quick',
                    '-database=s',
                    '-engine=s',
                    '-recoverDir=s',
                    '-pa=i',
                    '-threads=i',
                    '-srand=i',
                    '-rsSampleSize=s',
                    '-numAddlRounds=s',
                    '-genomeSampleSizeMax=s',
                    '-LTRStruct',
                    '-LTRMaxSeqLen=i',
                    '-onlyRS',
                    '-skipRS',
                    '-disableRMBlastStat',
                    '-onlyIns',                     # For internal use only
                    '-recon_batch_size=i',          # For debugging only 
                    '-recon_sample_size_start=i',   # For debugging only
);

# Add configuration parameters as additional command-line options
push @getopt_args, RepModelConfig::getCommandLineOptions();

#
# Get the supplied command line options, and set flags
#
my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

$DEBUG = 1 if ( exists $options{'debug'} );

sub usage {
  my $p = Pod::Text->new();
  $p->output_fh( *STDOUT );
  my $pod_str;
  open IN, "<$0"
      or die "Could not open self ($0) for generating documentation!";
  while ( <IN> ) {
    if ( /^=head1\s+CONFIGURATION OVERRIDES\s*$/ ) {
      my $c_pod = RepModelConfig::getPOD();
      if ( $c_pod ) {
        $pod_str .= $_ . $c_pod;
      }
    }
    else {
      $pod_str .= $_;
    }
  }
  close IN;
  print "$0 - $version\n";
  $p->parse_string_document( $pod_str );
  exit( 1 );
}

#
# Resolve configuration settings using the following precedence:
# command line first, then environment, followed by config
# file.
#
RepModelConfig::resolveConfiguration( \%options );
my $config           = $RepModelConfig::configuration;
my $NCBIDBCMD_PRGM   = $config->{'RMBLAST_DIR'}->{'value'} . "/blastdbcmd";
my $NCBIBLASTDB_PRGM = $config->{'RMBLAST_DIR'}->{'value'} . "/makeblastdb";
my $NCBIDBALIAS_PRGM =
    $config->{'RMBLAST_DIR'}->{'value'} . "/blastdb_aliastool";
my $DUSTMASKER_PRGM = $config->{'RMBLAST_DIR'}->{'value'} . "/dustmasker";
my $RMBLASTN_PRGM   = $config->{'RMBLAST_DIR'}->{'value'} . "/rmblastn";
my $XDFORMAT_PRGM   = $config->{'ABBLAST_DIR'}->{'value'} . "/xdformat";
my $XDGET_PRGM      = $config->{'ABBLAST_DIR'}->{'value'} . "/xdget";
my $WUBLASTN_PRGM   = $config->{'ABBLAST_DIR'}->{'value'} . "/blastn";
my $TRF_DIR         = $config->{'TRF_DIR'}->{'value'};
my $RSCOUT_DIR      = $config->{'RSCOUT_DIR'}->{'value'};
my $RECON_DIR       = $config->{'RECON_DIR'}->{'value'};
my $CDHIT_DIR       = $config->{'CDHIT_DIR'}->{'value'};

## validate params
# Minimal commands
foreach my $param ( "RMBLAST_DIR", "TRF_DIR", "RSCOUT_DIR", "RECON_DIR" ) {
  if ( !RepModelConfig::validateParam( $param ) ) {
    die
"\n  RepeatModeler dependency missing or incorrectly set for $param!\n  Rerun ./configure or check your command line to ensure that RepeatModeler\n  has access to and the correct version of this dependency.\n\n";
  }
}

# LTRPipeline dependencies
if ( $options{'LTRStruct'} ) {
  foreach my $param (
                      "GENOMETOOLS_DIR", "LTR_RETRIEVER_DIR",
                      "MAFFT_DIR",       "NINJA_DIR",
                      "CDHIT_DIR"
      )
  {
    if ( !RepModelConfig::validateParam( $param ) ) {
      die
"\n  LTRPipeline dependency missing or incorrectly set for $param!\n   Rerun ./configure or check your command line to ensure that RepeatModeler\n  has access to and the correct version of this dependency.\n\n";
    }
  }
  unless ( eval "require threads;" ) 
  {
    die "LTRRetriever used in the LTRPipeline depends on perl threads.  This version of\n".
        "perl does not appear to be compiled with the thread option.\n";
  }
}

if ( $options{'version'} ) {
  print "RepeatModeler version $version\n";
  exit;
}

# Print the internal POD documentation if something is missing
if ( !defined $options{'database'} || $options{'help'} ) {
  print "No database indicated\n\n";
  usage();
}

if ( $options{'pa'} ) {
  die "\nERROR: The -pa parameter has been deprecated.  Please use the newer\n" .
        "      \"-threads #\" parameter which precisely controls the maximum\n" .
        "       number of simultaneous threads this program will use. Note that\n" .
        "       the -pa parameter previously controlled the number of simultaneous\n" .
        "       batches each run with rmblast using 4 threads, therefore -pa 4 could\n" .
        "       potentially have used up to 16 threads.  Keep this in mind when\n" .
        "       comparing runtime differences between versions.\n\n";
}

if ( ! exists $options{'threads'} || $options{'threads'} < 2 ) {
  print "WARNING: RepeatModeler is a computationally intensive program.\n" .
        "         It is recommended that for anything other than debugging\n" .
        "         purposes the program should run with greater than eight\n" .
        "         threads (-threads #).\n\n";
}
# Make it easier to pass to functions
my $threads = $options{'threads'} ? $options{'threads'} : 1;

# Seed the random number generator if requested
my $seed = time;
$seed = $options{'srand'} if ( defined $options{'srand'} );
srand( $seed );

if ( exists $options{'disableRMBlastStat'} ) {
  $ENV{'BLAST_USAGE_REPORT'} = "false";
}

#
# Setup the search engines
#
my $searchEngineN;
my $engine = "rmblast";
my $engineVersion = "";
my $isRMBlastQueryThreadable = 0; # A flag that is set if RMBlastn supports query threading (2.13 and above)

# DEPRECATED: Using only rmblast now
if ( $options{'engine'} && $options{'engine'} =~ /wublast|abblast/i ) {
  die "ERROR: \"-engine $options{'engine'}\" is deprecated, this verison of RepeatModeler uses rmblast only.\n";
}

my $genomeDB = $options{'database'};
if ( -s $genomeDB ) {
  $genomeDB = File::Spec->rel2abs($genomeDB);
}elsif ( -s "$genomeDB.nsq" ){
  $genomeDB = File::Spec->rel2abs("$genomeDB.nsq");
}
$genomeDB =~ s/(.+)\.n[nihs][rndiq]$/$1/;
my $dbIndex;
$dbIndex = `$NCBIDBCMD_PRGM -db $genomeDB -entry all -outfmt "%i %l"`;

sub log_print {
  my $string = shift;
  print "$string";
  print $LOG "$string";
}
elapsedTime( "runtime" ); 
#chdir( $tmpDir );
if ( $options{'LTRStruct'} ) {
  elapsedTime( 1 );
  log_print "\n\nLTR Structural Analysis\n";
  log_print "=======================\n";

  # Save genome to sequence file
  open OUT, ">$tmpDir/tmpInputSeq"
      or die "Could not open $tmpDir/tmpInputSeq for writing!\n";
  open SEQ, "$NCBIDBCMD_PRGM -db $genomeDB -entry all -outfmt \"%f\"|"
      or die "Could not run: $NCBIDBCMD_PRGM\n";
  while ( <SEQ> ) {
    if ( /^>/ ) {

      # Replace >gi|232 with >gi-232 so that LTR_retriever is happy.  This is
      # a known issue with LTR_retriever-2.6
      s/\|/-/;
    }
    print OUT;
  }
  close SEQ;
  close OUT;
  if ( -s "$tmpDir/tmpInputSeq" ) {
    my $optionalArgs = " -giToID $genomeDB.translation";
    $optionalArgs .= " -LTRMaxSeqLen $options{'LTRMaxSeqLen'}"
        if ( exists $options{'LTRMaxSeqLen'} );
    $optionalArgs .= " -threads $options{'threads'}" if ( exists $options{'threads'} );
    $optionalArgs .= " -debug"             if ( $DEBUG );
    ## Run the LTR Pipeline
    print "Calling LTRPipeline as: $FindBin::RealBin/LTRPipeline $optionalArgs $tmpDir/tmpInputSeq\n" if ( $DEBUG );
    system( "$FindBin::RealBin/LTRPipeline $optionalArgs $tmpDir/tmpInputSeq" );
    unlink "$tmpDir/tmpInputSeq" unless ( $DEBUG );
    if ( -s "$tmpDir/tmpInputSeq-ltrs.fa" ) {

      # Combine results between both pipelines
      # run CDhit
      # /usr/local/cd-hit/cd-hit-est -aS 0.8 -c 0.8 -g 1 -G 0 -A 80 -M 10000
      #         -i all_consensi.fasta -o all_consensi.clusters1.fasta -T 20
      #   T is threads, 0 = unlimited
      #
      # >Cluster 0
      # 0	8283nt, >rnd-3_family-2#LINE... *
      # 1	3592nt, >rnd-4_family-1#LINE... at 1:3592:4587:8184/+/98.72%
      # 2	2990nt, >rnd-5_family-1#LINE... at 1:2990:4994:7985/+/97.46%
      # 3	1641nt, >rnd-5_family-22#LIN... at 2:1623:10:1632/+/87.13%
      #
      # Jullien clustered and removed RECON/RepeatScout hits that overlapped
      # Then she clustered a second time and collapsed clusters using Refiner.
      #
      # What we are doing:
      #   Cluster the sequences with RepeatModeler's RECON/RepeatScout pipeline.
      #   Identify non-singletons and collapse in favor of the longest sequence
      #   ( for RECON/RepeatScout clusters ) or in favor of LTRPipeline in mixed
      #   clusters.  Keep all LtrPipeline produced results.
      #
      log_print "  -- Clustering results with previous rounds...\n";

      # 1. Combine results from both pipelines into a file
      system(
"cat $tmpDir/tmpInputSeq-ltrs.fa $tmpDir/consensi.fa > $tmpDir/combined.fa" );
      system(
"cat $tmpDir/tmpInputSeq-ltrs.stk $tmpDir/families.stk > $tmpDir/combined.stk"
      );

      # 2. Cluster results
      my $cmd =
            "$CDHIT_DIR/cd-hit-est -aS 0.8 -c 0.8 -g 1 -G 0 -A 80 -M 10000 "
          . "-i $tmpDir/combined.fa -o $tmpDir/cd-hit-out -T $threads "
          . "> $tmpDir/cd-hit-stdout 2>&1";
      system( $cmd);
      unlink( "$tmpDir/cd-hit-stdout" ) if ( -e "$tmpDir/cd-hit-stdout" );
      unlink( "$tmpDir/cd-hit-out" )    if ( -e "$tmpDir/cd-hit-out" );

      # 3. Process clusters and remove redundancy
      my %redundant_families;
      my %putative_subfamilies;
      my $ltrFamCnt = 0;
      my $rrFamCnt  = 0;
      if ( -s "$tmpDir/cd-hit-out.clstr" ) {
        open IN, "<$tmpDir/cd-hit-out.clstr"
            or die
"RepeatModeler: Could not open $tmpDir/cd-hit-out.clstr for reading!\n";
        my @cluster = ();
        my $longest_ltr_id;
        my $longest_ltr_size;
        my $longest_rnd_id;
        my $longest_rnd_size;
        while ( <IN> ) {
          if ( /^>Cluster/ ) {
            if ( @cluster > 1 ) {

              # I started out by picking the longest LTRPipeline
              # derived family as the cluster rep.  Now we
              # consider all LTRPipeline candidates as dominant.
              # The only reason two more more LTRPipeline candidates
              # would show up in a cluster would be due to
              # overlap between LTR/INT sequences that didn't get
              # removed by LTR_retriever.
              #
              # Keeping the longest_ltr_id variable as an indicator
              # that there is at least one LTRPipeline candidate in the
              # cluster.
              if ( $longest_ltr_id ) {
                foreach my $id ( @cluster ) {

                  # Now...remove only overlapping RECON/RepeatScout
                  # candidates ( rnd-#_family-# ).
                  if ( $id !~ /^ltr-\d+_family-\d+/ ) {

            # Use the longest as the "reason" why we are removing it.
            #print "Removing $id because it's redundant with $longest_ltr_id\n";
                    $redundant_families{$id}++;
                  }
                }
              }
              elsif ( $longest_rnd_id ) {
                foreach my $id ( @cluster ) {
                  if ( $id ne $longest_rnd_id ) {

               #print "Labeling $id as putative subfamily of $longest_rnd_id\n";
                    $putative_subfamilies{$id} = $longest_rnd_id;
                  }
                }
              }
            }
            $longest_ltr_id   = undef;
            $longest_ltr_size = 0;
            $longest_rnd_id   = undef;
            $longest_rnd_size = 0;
            @cluster          = ();
            next;
          }
          if ( /^\d+\s+(\d+)nt,\s*>((rnd|ltr)-\d+_family-\d+)/ ) {
            my $size = $1;
            my $id   = $2;
            my $type = $3;
            if ( $type eq "ltr" ) {
              $ltrFamCnt++;
              if ( $size > $longest_ltr_size ) {
                $longest_ltr_id   = $id;
                $longest_ltr_size = $size;
              }
            }
            if ( $type eq "rnd" ) {
              $rrFamCnt++;
              if ( $size > $longest_rnd_size ) {
                $longest_rnd_id   = $id;
                $longest_rnd_size = $size;
              }
            }
            push @cluster, $id;
          }
        }
        close IN;

        # TODO: Do we really want to keep this?
        #unlink("$tmpDir/cd-hit-out.clstr" );
      }
      log_print "       - $rrFamCnt RepeatScout/RECON families\n";
      log_print "       - $ltrFamCnt LTRPipeline families\n";
      if ( keys( %redundant_families ) ) {
        system(
               "mv $tmpDir/consensi.fa $tmpDir/consensi.fa.recon_rscout_only" );
        system( "mv $tmpDir/combined.fa $tmpDir/consensi.fa.with_redundancy" );

        # Filter consensi.fa and families.stk
        open IN, "<$tmpDir/consensi.fa.with_redundancy"
            or die
"RepeatModeler: Could not open $tmpDir/consensi.fa.with_redundancy for reading";
        open OUT, ">$tmpDir/consensi.fa"
            or die
            "RepeatModeler: Could not open $tmpDir/consensi.fa for writing";
        my $id;
        my $data;
        while ( <IN> ) {
          if ( /^>(\S+)/ ) {
            my $tmpID = $1;
            if ( $data ) {
              if ( !exists $redundant_families{$id} ) {
                print OUT $data;
              }
            }
            if ( exists $putative_subfamilies{$tmpID} ) {
              my $tstr = $_;
              $tstr =~ s/[\n\r]//g;
              $data =
                  "$tstr [ putative subfamily of "
                  . $putative_subfamilies{$tmpID} . " ]\n";
            }
            else {
              $data = $_;
            }
            $id = $tmpID;
            next;
          }
          $data .= $_;
        }
        if ( $data ) {
          if ( !exists $redundant_families{$id} ) {
            print OUT $data;
          }
        }
        close IN;
        close OUT;
        system(
             "mv $tmpDir/families.stk $tmpDir/families.stk.recon_rscout_only" );
        system(
               "mv $tmpDir/combined.stk $tmpDir/families.stk.with_redundancy" );
        open IN, "<$tmpDir/families.stk.with_redundancy"
            or die
"RepeatModeler: Could not open $tmpDir/families.stk.with_redundancy for reading";
        open OUT, ">$tmpDir/families.stk"
            or die
            "RepeatModeler: Could not open $tmpDir/families.stk for writing";
        $id   = "";
        $data = "";

        while ( <IN> ) {
          if ( /^#=GF\s+ID\s+(\S+)/ ) {
            $id = $1;
          }

          if ( /^#=GF\s+DE\s+(\S.*)/ ) {
            if ( exists $putative_subfamilies{$id} ) {
              my $tstr = $_;
              $tstr =~ s/[\n\r]//g;
              $data .=
                  "$tstr [ putative subfamily of "
                  . $putative_subfamilies{$id} . " ]\n";
            }
            else {
              $data .= $_;
            }
          }
          else {
            $data .= $_;
          }

          if ( /^\/\// ) {
            if ( !exists $redundant_families{$id} ) {
              print OUT $data;
            }
            $data = "";
          }
        }
        close IN;
        close OUT;
        log_print "       - Removed "
            . scalar( keys( %redundant_families ) )
            . " redundant LTR families.\n";
        log_print "       - Final family count = "
            . (
            ( $rrFamCnt + $ltrFamCnt ) - scalar( keys( %redundant_families ) ) )
            . "\n";
      }
    }
  }
  else {
    log_print
"\nWARNING: Could not create input file for LTRPipeline from $genomeDB! Continuing using\nresults from RepeatScout/RECON pipeline only.\n";
  }
  log_print "LTRPipeline Time: " . elapsedTime( 1 ) . "\n";
}

my $numModels = 6;
elapsedTime( 1 );
if ( $numModels > 0 ) {
  log_print "\n\n";
  my $origDir = getcwd();
  chdir( $tmpDir );
  my $addlOpts = "";
  $addlOpts = "-threads $options{'threads'} " if ( $options{'threads'} );
  system(   "$FindBin::RealBin/RepeatClassifier "
          . "-consensi consensi.fa -stockholm families.stk" );
  chdir( $origDir );
  log_print "Classification Time: " . elapsedTime( 1 ) . "\n";
}

log_print "\n\nProgram Time: " . elapsedTime( "runtime" ) . "\n";
print "Working directory:  $tmpDir\n";
print "may be deleted unless there were problems with the run.\n";

if ( $numModels > 0 ) {
  system( "cp $tmpDir/consensi.fa.classified $genomeDB-families.fa" )
      if ( -s "$tmpDir/consensi.fa.classified" );
  system( "cp $tmpDir/families-classified.stk $genomeDB-families.stk" )
      if ( -s "$tmpDir/families-classified.stk" );
  system( "cp $tmpDir/rmod.log $genomeDB-rmod.log" )
      if ( -s "$tmpDir/rmod.log" );
  print "\nThe results have been saved to:\n";
  print
"  $genomeDB-families.fa  - Consensus sequences for each family identified.\n";
  print
"  $genomeDB-families.stk - Seed alignments for each family identified.\n";
  print
"  $genomeDB-rmod.log     - Execution log.  Useful for reproducing results.\n\n";
  print "The RepeatModeler stockholm file is formatted so that it can\n";
  print
"easily be submitted to the Dfam database.  Please consider contributing\n";
  print
      "curated families to this open database and be a part of this growing\n";
  print
      "community resource.  For more information contact help\@dfam.org.\n\n\n";
}
else {
  print "No families identified.  Perhaps the database is too small\n";
  print "or contains overly fragmented sequences.\n";
}
sub elapsedTime {
  my ( $TimeHistIdx ) = @_;
  if ( defined $TimeBefore{$TimeHistIdx} ) {
    my $DiffTime = time - $TimeBefore{$TimeHistIdx};
    $TimeBefore{$TimeHistIdx} = time;
    my $Min = int( $DiffTime / 60 );
    $DiffTime -= $Min * 60;
    my $Hours = int( $Min / 60 );
    $Min -= $Hours * 60;
    my $Sec = $DiffTime;
    my $timeStr = sprintf( "%02d:%02d:%02d", $Hours, $Min, $Sec );
    return "$timeStr (hh:mm:ss) Elapsed Time";
  }
  else {
    $TimeBefore{$TimeHistIdx} = time;
    return 0;
  }
}
exit;
