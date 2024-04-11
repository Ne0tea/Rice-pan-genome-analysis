library(ggplot2)
library(tidyr)
library(dplyr)
library(ggsci)

family_platte<-c("Gypsy"='#F2AE2C',"Copia"='#FFE38D',"L1"='#cbb440',"tRNA"='#fff3a2',"Caulimovirus"='#fdfe96',
                 "MULE-MuDR"='#bdebfb',"hAT"='#c7f7f9',"CMC-EnSpm"='#718dcc',"PIF-Harbinger"='#233597',
                 "TcMar-Stowaway"='#0d374f',"Helitron"='#3BA738')

phy_platte<-c('or4'='#76D273','or1'='#51C54E','or2'='#3BA738',"or3"='#215A20','trp'='#D0B34E',
              'tmp'='#F2AE2C','XI3'='#53C2C6','XI2'='#3AA8AC','XI1A'='#2E8285','XI1B'='#225C5F',
              'aus'='#4E84C3','obar'='#C3DEAF','ogla'='#FFE38D')

phy_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/70_phylogenic.txt"
phy_df<-read.table(phy_file,col.names = c("sp","group"))
sp_group<-c(phy_df$group)
names(sp_group)<-c(phy_df$sp)

Sp_length_file<-"E:/Bio_analysis/Weedyrice/Sfamily_spchr.txt"
Sp_length_df<-read.table(Sp_length_file,check.names=FALSE)
colnames(Sp_length_df)<-c('Sp','Chr','superfamily','length','clength','loc')
Sp_length_df$phy<-sp_group[Sp_length_df$Sp]
Sp_length_df$Sp<-factor(Sp_length_df$Sp, levels = rev(names(sp_group)))
Sp_length_df<-Sp_length_df[Sp_length_df$loc=='C',]
for (i in unique(Sp_length_df$Chr)) {
  sf_plot<-ggplot(data=Sp_length_df[Sp_length_df$Chr==i,])+
    geom_point(aes(x=clength,y=length,color=phy))+
    # xlim(c(0,100000))+
    # ylim(0, 200000)+
    scale_color_manual(values = phy_platte)+
    theme_classic()+
    guides(fill = guide_legend(ncol=1))
  ggsave(paste(i,'_phy_ratio.pdf',sep = ''),height = 7,width = 8,path="E:/Bio_analysis/Weedyrice/pan_weedyrice/1207_statresult/Centro_region")   
  
  
  ratio_df <- Sp_length_df[Sp_length_df$Chr==i,] %>%
    group_by(Sp) %>%
    summarize(total_sum = sum(length))
  ratio_df<-merge(unique(Sp_length_df[Sp_length_df$Chr==i,c('Sp','clength')]),ratio_df)
  ratio_df$Sp<-factor(ratio_df$Sp, levels = rev(names(sp_group)))
  sp_plot<-ggplot()+  
    geom_bar(data=ratio_df,
             aes(x=Sp,y=clength),stat="identity", fill="lightgrey", width=0.7,linewidth=0.25,inherit.aes=FALSE)+
    geom_line(data=ratio_df,aes(x=Sp,y=total_sum/clength*max(clength),group=1),color='grey',linewidth=0.75,inherit.aes=FALSE)+
    geom_point(data=ratio_df,aes(x=Sp,y=total_sum/clength*max(clength)),shape=19,size=1.25,inherit.aes=FALSE)+
    # ggplot(data=freq_sp_df,aes(Sp,value,fill=name))+
    geom_bar(data=Sp_length_df[Sp_length_df$Chr==i,],aes(Sp,length,fill=superfamily),
             stat="identity",position="stack",width=0.7,linewidth=0.25)+
    scale_fill_manual(values = family_platte)+
    theme_classic()+
    guides(fill = guide_legend(ncol=1))+
    coord_flip()
  ggsave(paste(i,'_Sfamily_sp.pdf',sep = ''),height = 7,width = 8,path="E:/Bio_analysis/Weedyrice/pan_weedyrice/1207_statresult/Centro_region") 
  }
sf_plot<-ggplot(data=Sp_length_df)+
  geom_point(aes(x=clength,y=length,color=phy))+
  # xlim(c(0,100000))+
  ylim(0, 200000)+
  scale_color_manual(values = phy_platte)+
  theme_classic()+
  guides(fill = guide_legend(ncol=1))
sf_plot



