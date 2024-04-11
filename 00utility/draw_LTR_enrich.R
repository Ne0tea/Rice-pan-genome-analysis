#绘制百分比的堆叠图

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggsci)
library(ggsignif)

family_platte<-c("Gypsy"='#F2AE2C',"Copia"='#FFE38D',"L1"='#cbb440',"tRNA"='#fff3a2',"Caulimovirus"='#fdfe96',
                 "MULE-MuDR"='#bdebfb',"hAT"='#c7f7f9',"CMC-EnSpm"='#718dcc',"PIF-Harbinger"='#233597',
                 "TcMar-Stowaway"='#0d374f',"Helitron"='#3BA738')
phy_platte<-c('or4'='#76D273','or1'='#51C54E','or2'='#3BA738',"or3"='#215A20','trp'='#D0B34E','tmp'='#F2AE2C',
              'XI3'='#53C2C6','XI2'='#3AA8AC','XI1A'='#2E8285','XI1B'='#225C5F','aus'='#4E84C3','obar'='#C3DEAF','ogla'='#FFE38D')
####phylogentic relationship
phy_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/70_phylogenic.txt"
phy_df<-read.table(phy_file,col.names = c("sp","group"))
sp_group<-c(phy_df$group)
names(sp_group)<-c(phy_df$sp)
####centro length as background
centro_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/centromere_region_70m_used.bed"
centro_df<-read.table(centro_file,col.names = c("chr_sp","start","end"))
centro_df <- separate(centro_df, chr_sp, into = c("chr", "sp"), sep = "_")
centro_df$length<-centro_df$end-centro_df$start
centro_df <- centro_df %>%
  group_by(sp) %>%
  summarize(mean_value = sum(length, na.rm = TRUE))
colnames(centro_df)<-c('sp','length')
centro_df$sp<-factor(centro_df$sp, levels = rev(names(sp_group)))

enrich_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/family_length_centro.csv"
family_length_df<-read.csv(enrich_file,header = 1,check.names=FALSE)
colnames(family_length_df)[1]<-'X'
rownames(family_length_df)<-family_length_df$X
family_length_long <- gather(family_length_df, key = "sp",value = 'length',-'X')
colnames(family_length_long)<-c('Sfamily','Sp','Length')
family_length_long<-family_length_long[family_length_long$Sfamily!='None' & family_length_long$Sfamily!='Unknown',]
if (1) {
  family_length_long$Sp<-factor(family_length_long$Sp, levels = rev(names(sp_group)))
  # sp_plot<-ggplot(data=centro_df,aes(sp,length/48000))+
  sp_plot<-ggplot(data=centro_df,aes(sp,length/2))+
    geom_rect(xmin=0,xmax=9.5,ymin=-Inf,ymax=2,
              fill=phy_platte['XI1A'],alpha=0.3)+
    geom_rect(xmin=9.5,xmax=17.5,ymin=-Inf,ymax=2,
              fill=phy_platte['XI1B'],alpha=0.3)+
    geom_rect(xmin=17.5,xmax=24.5,ymin=-Inf,ymax=2,
              fill=phy_platte['XI3'],alpha=0.3)+
    geom_rect(xmin=24.5,xmax=27.5,ymin=-Inf,ymax=2,
              fill=phy_platte['XI2'],alpha=0.3)+
    geom_rect(xmin=27.5,xmax=33.5,ymin=-Inf,ymax=2,
              fill=phy_platte['aus'],alpha=0.3)+
    geom_rect(xmin=33.5,xmax=35.5,ymin=-Inf,ymax=2,
              fill=phy_platte['or2'],alpha=0.3)+
    geom_rect(xmin=35.5,xmax=37.5,ymin=-Inf,ymax=2,
              fill=phy_platte['or1'],alpha=0.3)+
    geom_rect(xmin=37.5,xmax=51.5,ymin=-Inf,ymax=2,
              fill=phy_platte['tmp'],alpha=0.1)+
    geom_rect(xmin=51.5,xmax=55.5,ymin=-Inf,ymax=2,
              fill=phy_platte['trp'],alpha=0.1)+
    geom_rect(xmin=55.5,xmax=61.5,ymin=-Inf,ymax=2,
              fill=phy_platte['or3'],alpha=0.1)+
    geom_rect(xmin=61.5,xmax=62.5,ymin=-Inf,ymax=2,
              fill=phy_platte['or4'],alpha=0.1)+
    geom_rect(xmin=62.5,xmax=66.5,ymin=-Inf,ymax=2,
              fill=phy_platte['ogla'],alpha=0.1)+
    geom_rect(xmin=66.5,xmax=70.5,ymin=-Inf,ymax=2,
              fill=alpha(phy_platte['obar'],0.1))+
    geom_bar(stat="identity", fill="lightgrey", width=0.7,linewidth=0.25)+
    # ggplot(data=freq_sp_df,aes(Sp,value,fill=name))+
    geom_bar(data=na.omit(family_length_long[family_length_long$Length!=0,]),
             aes(Sp,Length,fill=Sfamily),stat="identity",position="stack",width=0.7,linewidth=0.25,inherit.aes=FALSE,)+
    scale_fill_manual(values = family_platte)+
    theme_classic()+
    guides(fill = guide_legend(ncol=1))+
    coord_flip()
  ggsave('04_Sp_Superfamily_length.pdf',sp_plot,path = paste(dirname(enrich_file),'/1207_statresult',sep = ''),height = 12,width = 7) 
  XI<-names(sp_group[sp_group=='XI1A' | sp_group=='XI1B'])
  Geng<-names(sp_group[sp_group=='tmp' | sp_group=='trp'])
  XIvGe<-data.frame('var'=c(rep('XI',17),rep('Geng',18)),
                    'value'=c(family_length_long[family_length_long$Sp %in% XI & family_length_long$Sfamily=='Gypsy','Length'],
                                     family_length_long[family_length_long$Sp %in% Geng & family_length_long$Sfamily=='Gypsy','Length']))
  plot <- ggplot(XIvGe, aes(x = var, y = value,fill=var)) +
    geom_bar(stat = "summary", position = "dodge",color='black') +  # 柱状图
    geom_point(shape=21,position = position_jitter(width = 0.1), color = "black") +  # 散点图
    stat_summary(fun.data = 'mean_se', geom = "errorbar", colour = "black",
                 width = 0.2,position = position_dodge(1),linewidth = 0.75)+
    # labs(title = "柱状散点图", x = "Group", y = "Value") +  # 标题和轴标签
    geom_signif(comparisons = list(c('XI','Geng')),
                step_increase = 0.3,tip_length=0,
                map_signif_level = T,
                test = t.test)+
    theme_classic()+
    scale_fill_manual(values = c('Geng'='#F2AE2C','XI'='#4E84C3'))
  ggsave('04_XIvsGeng_significant.pdf',plot,path = paste(dirname(enrich_file),'/1207_statresult',sep = ''),height = 4,width = 4) 
}

if (1) {
  enrich_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/family_length_chr_centro.csv"
  family_length_chr_df<-read.csv(enrich_file,header = 1,check.names=FALSE)
  family_length_chr_df<-family_length_chr_df[-1,]
  colnames(family_length_chr_df)[1]<-'X'
  rownames(family_length_chr_df)<-family_length_chr_df$X
  family_length_chr_long <- gather(family_length_chr_df, key = "Chr",value = 'length',-'X')
  colnames(family_length_chr_long)<-c('Sfamily','Chr','Length')  
  
  centro_df<-read.table(centro_file,col.names = c("chr_sp","start","end"))
  centro_df <- separate(centro_df, chr_sp, into = c("chr", "sp"), sep = "_")
  centro_df$length<-centro_df$end-centro_df$start
  centro_df <- centro_df %>%
    group_by(chr) %>%
    summarize(mean_value = sum(length, na.rm = TRUE))
  colnames(centro_df)<-c('chr','length')
  
  sp_plot<-ggplot(data=centro_df,aes(chr,length))+  
  geom_bar(stat="identity", fill="lightgrey", width=0.7,linewidth=0.25)+
    # ggplot(data=freq_sp_df,aes(Sp,value,fill=name))+
    geom_bar(data=family_length_chr_long[family_length_chr_long$Length!=0,],aes(Chr,Length,fill=Sfamily),
             stat="identity",position="stack",width=0.7,linewidth=0.25,inherit.aes=FALSE,)+
    scale_fill_manual(values = family_platte)+
    theme_classic()+
    guides(fill = guide_legend(ncol=1))+
    coord_flip()
  ggsave('04_Chr_Superfamily_length_C.pdf',sp_plot,path = paste(dirname(enrich_file),'/1207_statresult',sep = ''),height = 5,width = 7) 
}

ref_fai_file='E:/Bio_analysis/Weedyrice/pan_weedyrice/all_asm_70.fa.fai'
enrich_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/family_length_chr_OUT_centro.csv"
family_length_chr_df<-read.csv(enrich_file,header = 1,check.names=FALSE)
family_length_chr_df<-family_length_chr_df[-1,]
colnames(family_length_chr_df)[1]<-'X'
rownames(family_length_chr_df)<-family_length_chr_df$X
family_length_chr_long <- gather(family_length_chr_df, key = "Chr",value = 'length',-'X')
colnames(family_length_chr_long)<-c('Sfamily','Chr','Length')  

chr_length_df<-read.table(ref_fai_file,col.names = c("chr_sp","length","_","_1","_2"))
chr_length_df <- separate(chr_length_df, chr_sp, into = c("chr", "sp"), sep = "_")
chr_length_df <- chr_length_df %>%
  group_by(chr) %>%
  summarize(mean_value = sum(length, na.rm = TRUE))
colnames(chr_length_df)<-c('chr','length')
chr_length_df$length<-chr_length_df$length-centro_df$length

sp_plot<-ggplot(data=chr_length_df,aes(chr,length))+  
  geom_bar(stat="identity", fill="lightgrey", width=0.7,linewidth=0.25)+
  # ggplot(data=freq_sp_df,aes(Sp,value,fill=name))+
  geom_bar(data=family_length_chr_long[family_length_chr_long$Length!=0 ,],aes(Chr,Length,fill=Sfamily),
           stat="identity",position="stack",width=0.7,linewidth=0.25,inherit.aes=FALSE,)+
  scale_fill_manual(values = family_platte)+
  theme_classic()+
  guides(fill = guide_legend(ncol=1))+
  coord_flip()
ggsave('04_Chr_Superfamily_length_OC.pdf',sp_plot,path = paste(dirname(enrich_file),'/1207_statresult',sep = ''),height = 5,width = 7)
