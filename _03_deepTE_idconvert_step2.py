import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

fr = open(input_file,'r')
fw = open(output_file,'w')

for line in fr:
    if len(line.strip().split("\t")[1].split("_")) >= 2 and line.strip().split("\t")[1].split("_")[1] == "DNA":
        prefix = line.strip().split("\t")[0].split("#")[0]
        suffix = "DNA/"+"-".join(line.strip().split("\t")[1].split("_")[2:])
        fw.write("%s\t%s\n"%(line.strip().split("\t")[0],prefix+"#"+suffix))
        
        
#print(i)
fw.close()
fr.close()
