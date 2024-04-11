library(ggplot2)
family_platte<-c("Gypsy"='#F2AE2C',"Copia"='#FFE38D',"L1"='#cbb440',"tRNA"='#fff3a2',"Caulimovirus"='#fdfe96',
                 "MULE-MuDR"='#bdebfb',"hAT"='#c7f7f9',"CMC-EnSpm"='#718dcc',"PIF-Harbinger"='#233597',
                 "TcMar-Stowaway"='#0d374f',"Helitron"='#3BA738','Others'='grey')
enrich_C_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/family_length_chr_centro.csv"
family_length_chr_C_df<-read.csv(enrich_C_file,header = 1,check.names=FALSE)
family_length_chr_C_df<-family_length_chr_C_df[-1,]
colnames(family_length_chr_C_df)[1]<-'X'
rownames(family_length_chr_C_df)<-family_length_chr_C_df$X
family_length_chr_C_long <- gather(family_length_chr_C_df, key = "Chr",value = 'length',-'X')
colnames(family_length_chr_C_long)<-c('Sfamily','Chr','Length')  

enrich_OC_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/family_length_chr_OUT_centro.csv"
family_length_chr_OC_df<-read.csv(enrich_OC_file,header = 1,check.names=FALSE)
family_length_chr_OC_df<-family_length_chr_OC_df[-1,]
colnames(family_length_chr_OC_df)[1]<-'X'
rownames(family_length_chr_OC_df)<-family_length_chr_OC_df$X
family_length_chr_OC_long <- gather(family_length_chr_OC_df, key = "Chr",value = 'length',-'X')
colnames(family_length_chr_OC_long)<-c('Sfamily','Chr','Length')  
if (0) {
  ref_fai_file='E:/Bio_analysis/Weedyrice/pan_weedyrice/all_asm_70.fa.fai'
  chr_length_df<-read.table(ref_fai_file,col.names = c("chr_sp","length","_","_1","_2"))
  chr_length_df <- separate(chr_length_df, chr_sp, into = c("chr", "sp"), sep = "_")
  chr_length_df <- chr_length_df %>%
    group_by(chr) %>%
    summarize(mean_value = sum(length, na.rm = TRUE))
  colnames(chr_length_df)<-c('Chr','length')
  chr_length_df$length<-chr_length_df$length-centro_df$length
  centro_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/centromere_region_70m_used.bed"
  centro_df<-read.table(centro_file,col.names = c("chr_sp","start","end"))
  centro_df <- separate(centro_df, chr_sp, into = c("chr", "sp"), sep = "_")
  centro_df$length<-centro_df$end-centro_df$start
  centro_df <- centro_df %>%
    group_by(chr) %>%
    summarize(mean_value = sum(length, na.rm = TRUE))
  colnames(centro_df)<-c('Chr','length') 
  gypsy_ratio<-merge(family_length_chr_C_long[family_length_chr_C_long$Sfamily=='Gypsy',],
                     family_length_chr_OC_long[family_length_chr_OC_long$Sfamily=='Gypsy',],by = c('Sfamily','Chr'))
  colnames(gypsy_ratio)[3:4]<-c('seq_C','seq_OC')
  gypsy_ratio<-merge(gypsy_ratio,centro_df,by='Chr')
  colnames(gypsy_ratio)[5]<-c('length_C')
  gypsy_ratio<-merge(gypsy_ratio,chr_length_df,by='Chr')
  colnames(gypsy_ratio)[6]<-c('length_OC')
  
  gypsy_ratio$seqC_ratio<-gypsy_ratio$seq_C/(gypsy_ratio$seq_C+gypsy_ratio$length_C)
  gypsy_ratio$lengthC_ratio<-gypsy_ratio$length_C/(gypsy_ratio$seq_C+gypsy_ratio$length_C)
  gypsy_ratio$seqOC_ratio<-gypsy_ratio$seq_OC/(gypsy_ratio$seq_OC+gypsy_ratio$length_OC)
  gypsy_ratio$lengthOC_ratio<-gypsy_ratio$length_OC/(gypsy_ratio$seq_OC+gypsy_ratio$length_OC)
  gypsy_ratio<-gypsy_ratio[,c('Chr','Sfamily','seqC_ratio','seqOC_ratio','lengthC_ratio','lengthOC_ratio')]
  }

family_length_chr_C_long<-family_length_chr_C_long[family_length_chr_C_long$Sfamily!='None',]
family_length_chr_OC_long<-family_length_chr_OC_long[family_length_chr_OC_long$Sfamily!='None',]
Sfamily_chr_C_ratio<-family_length_chr_C_long%>%group_by(Chr)%>%reframe(Sfamily=Sfamily,count = Length/sum(Length))
Sfamily_chr_OC_ratio<-family_length_chr_OC_long%>%group_by(Chr)%>%reframe(Sfamily=Sfamily,count = Length/sum(Length))
Sfamily_chr_C_ratio<-Sfamily_chr_C_ratio[Sfamily_chr_C_ratio$count>0,]
Sfamily_chr_OC_ratio<-Sfamily_chr_OC_ratio[Sfamily_chr_OC_ratio$count>0,]

for (i in Sfamily_chr_C_ratio$Chr) {
  c_line=Sfamily_chr_C_ratio[Sfamily_chr_C_ratio$Chr==i,]
  oc_line=Sfamily_chr_OC_ratio[Sfamily_chr_OC_ratio$Chr==i,]
  plot_df<-data.frame()
  last_count=0
  for (x in c_line[order(c_line$count, decreasing = TRUE),"Sfamily"]$Sfamily) {
    # print(x)
    count<-as.numeric(c_line[c_line$Sfamily==x,'count'])
    # print(last_count+count)
    plot_line<-c(Chr=i,name=x,type='C',ymax=last_count+count,ymin=last_count)
    
    last_count=last_count+count
    plot_df<-rbind(plot_df,plot_line)
    colnames(plot_df)<-c('Chr','name','type','ymax','ymin')
  }
  last_count=0
  for (x in oc_line[order(oc_line$count, decreasing = TRUE), "Sfamily"]$Sfamily) {
    count<-as.numeric(oc_line[oc_line$Sfamily==x,'count'])
    # print(last_count+count)
    plot_line<-c(Chr=i,name=x,type='OC',ymax=last_count+count,ymin=last_count)
    
    last_count=last_count+count
    plot_df<-rbind(plot_df,plot_line)
    colnames(plot_df)<-c('Chr','name','type','ymax','ymin')
  }
  plot_df[,c('ymax','ymin')]<-lapply(plot_df[,c('ymax','ymin')],as.numeric)
  xx<-ggplot() + 
    geom_rect(data = plot_df[plot_df$type=='C',],aes(fill=name, ymax=ymax, ymin=ymin, xmax=3, xmin=2),colour = "black") +
    geom_rect(data=plot_df[plot_df$type=='OC',],aes(fill=name, ymax=ymax, ymin=ymin, xmax=2, xmin=0),
              colour = "black",inherit.aes = FALSE) +
    scale_fill_manual(values = family_platte)+
    theme_bw()+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank())+
    coord_polar(theta="y")
  ggsave(paste(i,'_Gypsy_ratio.pdf',sep = ''),height = 3,width = 4,path="E:/Bio_analysis/Weedyrice/pan_weedyrice/1207_statresult") 
}

value<-c(16774021,27287185.5,18061795.5,9622945.5,1355974,
  6193262,7347303,85750766.5,12581505,9943,
  880038,136227.5,23280.5,27529.5,416039.5,
  44917.5,46210,20952,18685,8364,30661.5,185586119)
names<-c('Helitron', 'MULE-MuDR', 'CMC-EnSpm', 'PIF-Harbinger', 'TcMar-Stowaway', 
'L1', 'hAT', 'Gypsy', 'Copia', 'Maverick', 'tRNA', 'Unknown', 'ID', 'PiggyBac', 
'Caulimovirus', 'Pao', 'hAT-Tip100', 'Chlamys', 'RTE-BovB', 'RTE-RTE', 'hAT-Ac','Others')
pie(value, labels=names,
    radius = 1.0,clockwise=T,
    main = "Male individuals",col=family_platte[names])
names<-c('LTR', 'DNA', 'LINE', 'RC','Others')
value<-c(659228.5, 35589, 1371, 11227,3279563)
pie(value, labels=names, 
    radius = 1.0,clockwise=T,
    main = "Male individuals",col=c("LTR"='#F2AE2C',"LINE"='#FFE38D',"DNA"='#718dcc',"RC"='#3BA738','Others'='grey'))
names<-c('Gypsy', 'CMC-EnSpm', 'MULE-MuDR', 'Helitron', 'RTE-BovB', 
         'hAT', 'Copia', 'TcMar-Stowaway', 'L1', 'PIF-Harbinger', 
         'Pao', 'Chlamys','Others')
value<-c(655163.5, 11794, 15333.5, 11227, 26, 5421.5, 3996, 163, 1345, 2877, 69, 236,3279563)
pie(value, labels=names,
    radius = 1.0,clockwise=T,
    main = "Male individuals",col=family_platte[names])
names<-c('Gypsy', 'CMC-EnSpm', 'MULE-MuDR', 'Helitron', 'RTE-BovB', 
         'hAT', 'Copia', 'TcMar-Stowaway', 'L1', 'PIF-Harbinger', 
         'Pao', 'Chlamys','Others')
value<-c(1153945,9151000,6638391,7283194,20231571,1456846,851810)
names<-c('ATLANTYS', 'Gypsy', 'RETR', 'RIRE', 'SZ-22', 'other_SZ', 'TRUNCATOR')
pie(value, labels=names,
    radius = 1.0,clockwise=T,
    main = "Male individuals",col=c('ATLANTYS'='#718dcc', 'Gypsy'='#cbb440', 'RETR'='#FFE38D', 'RIRE'='#c7f7f9', 
                                    'SZ-22'='#F2AE2C', 'other_SZ'='#fff3a2', 'TRUNCATOR'='#233597'))