'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-12-21 22:53:02
LastEditors: Ne0tea
LastEditTime: 2024-03-03 13:03:39
'''
import os
import sys

def main(ful_file,c_file,porpotion):
    # 文件是否存在，以及是否为空
    porpotion=int(porpotion)
    if os.path.exists(ful_file):
        sz = os.path.getsize(ful_file)
        if not sz:
            ord=0
            chr_loc_count={}
            with open(c_file,'r') as cf:
                for i in cf:
                    line=i.strip().split()
                    count=int(line[3])/(int(line[2])-int(line[1]))*100
                    chr_loc_count[line[0].split('_')[0]+'_'+str(ord)]=count
                    ord+=1
                    if ord >= porpotion:
                        ord=0
            with open(ful_file,'w') as of:
                for i in chr_loc_count:
                    st=str(int(i.split('_')[1])*porpotion+1)
                    en=str(int(int(i.split('_')[1])+1)*porpotion)
                    w_line='\t'.join([i,st,en,str(chr_loc_count[i])])
                    of.write(w_line+'\n')
            return
        else:
            ord=0
            chr_loc_count={}
            with open(c_file,'r') as cf:
                for i in cf:
                    line=i.strip().split()
                    count=int(line[3])/(int(line[2])-int(line[1]))*100
                    chr_loc_count[line[0].split('_')[0]+'_'+str(ord)]=count
                    ord+=1
                    if ord >= porpotion:
                        ord=0
                    pass
            ful_loc_count={}
            with open(ful_file,'r') as ff:
                for i in ff:
                    line=i.strip().split()
                    count=float(line[3])/(int(line[2])-int(line[1]))*100
                    ful_loc_count[line[0]]=count
            with open('tmp_result.txt','w') as of:
                for i in ful_loc_count:
                    st=str(int(i.split('_')[1])*porpotion+1)
                    en=str(int(int(i.split('_')[1])+1)*porpotion)
                    w_line='\t'.join([i,st,en,str(ful_loc_count[i]+chr_loc_count[i])])
                    of.write(w_line+'\n')
            os.rename('tmp_result.txt', ful_file)
    else:
        ord=0
        chr_loc_count={}
        with open(c_file,'r') as cf:
            for i in cf:
                line=i.strip().split()
                count=int(line[3])/(int(line[2])-int(line[1]))*100
                chr_loc_count[line[0].split('_')[0]+'_'+str(ord)]=count
                ord+=1
                if ord >= porpotion:
                    ord=0
        with open(ful_file,'w') as of:
            for i in chr_loc_count:
                st=str(int(i.split('_')[1])*porpotion+1)
                en=str(int(int(i.split('_')[1])+1)*porpotion)
                w_line='\t'.join([i,st,en,str(chr_loc_count[i])])
                of.write(w_line+'\n')
        return

if __name__ == "__main__":
    ful_file=sys.argv[1]
    c_file=sys.argv[2]
    porpotion=sys.argv[3]
    main(ful_file,c_file,porpotion)
