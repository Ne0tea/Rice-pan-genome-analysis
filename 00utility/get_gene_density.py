from Bio import SeqIO
import pandas as pd
class WaterMelon:
    def __init__(self,gtf,fasta,feature):
        self.gtf = gtf
        self.fasta = fasta
        self.feature = feature

    def chromoLen(self):
        chrLen = {}
        for rec in SeqIO.parse(self.fasta,'fasta'):
            chrLen[rec.id] = len(rec.seq)

        return chrLen

    def geneDensity(self,gene_window):
        self.gene_window = gene_window
        final_df = []
        df = pd.read_table(self.gtf,header=None,comment="#",sep="\t",
                           usecols=[0,2,3,4],
                           names="Chromosome Feature Start End".split())
        #df.columns = "Chromosome Source Feature Start End Score Strand Frame Attribute".split()
        if self.feature=='all':
            pass
        else:
            df = df[df.Feature==self.feature]
        chrLen = self.chromoLen()
        for chr_id in chrLen.keys():
            print(chr_id)
            df1 = df[df.Chromosome==chr_id]
            gene_start = [int(a) for a in df1.Start]
            gene_start.insert(0,0)
            gene_start.append(round(chrLen[chr_id]/self.gene_window)*self.gene_window+self.gene_window)
            #print(gene_start)
            bin_start = [int(a.left) for a in pd.cut(gene_start,bins=round(chrLen[chr_id]/self.gene_window)+1).value_counts().index]
            bin_start[0] = 0
            gene_count = list(pd.cut(gene_start,bins=round(chrLen[chr_id]/self.gene_window)+1).value_counts().values)
            #print(len(bin_start))
            #print(len(gene_count))
            #print("OK")
            final_df.append(pd.DataFrame({'chr_id':chr_id,'bin_start':bin_start,'gene_count':gene_count}))
            
        return pd.concat(final_df)