import sys

input_fa = sys.argv[1]
input_dict = sys.argv[2]
output_fa = sys.argv[3]


fr01 = open(input_fa,'r')
fr02 = open(input_dict,'r')

id_dict = {}

for line in fr02:
    id_dict[line.strip().split("\t")[0]] = line.strip().split("\t")[1]
    
fw = open(output_fa,'w')
    
for line in fr01:
    if line.startswith(">") and line.strip().split(" ")[0].replace(">","") in id_dict.keys():
        first_part = id_dict[line.strip().split(" ")[0].replace(">","")]
        second_part = ' '.join(line.strip().split(" ")[1:])
        fw.write(">"+first_part+" "+second_part+"\n")
    else:
        fw.write(line)
        
fr01.close()
fr02.close()
fw.close()
