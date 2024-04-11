library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
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
Sp_length_C_df<-Sp_length_df[Sp_length_df$loc=='C',]
Sp_length_OC_df<-Sp_length_df[Sp_length_df$loc=='OC',]

Sfamily_ratio_OC <- Sp_length_OC_df %>%
  group_by(Sp, Chr) %>%
  reframe(TE_length = sum(length)) %>%
  group_by(Sp) %>%
  reframe(TE_ratio = TE_length/sum(TE_length),Chr=Chr) %>%
  pivot_wider(names_from = Chr, values_from = TE_ratio)
Sfamily_ratio_OC[is.na(Sfamily_ratio_OC)]<-0

Sfamily_ratio_C <- Sp_length_C_df %>%
  group_by(Sp, Chr) %>%
  reframe(TE_length = sum(length)) %>%
  group_by(Sp) %>%
  reframe(TE_ratio = TE_length/sum(TE_length),Chr=Chr) %>%
  pivot_wider(names_from = Chr, values_from = TE_ratio)
Sfamily_ratio_C[is.na(Sfamily_ratio_C)]<-0

freq_col <- colorRamp2(
  breaks = c(0, 0.1, 0.2), 
  colors = c("#fefae0", "#dda15e", "#bc6c25")
)
Sfamily_ratio_C<-gather(Sfamily_ratio_C,key = Chr,value = value,-Sp)
Sfamily_ratio_OC<-gather(Sfamily_ratio_OC,key = Chr,value = value,-Sp)
# data4<-rbind(Sfamily_ratio_C,Sfamily_ratio_OC)

for (i in Sfamily_ratio_C$Chr) {
  chr_plot<-ggplot()+
    geom_line(data = Sfamily_ratio_C[Sfamily_ratio_C$Chr==i,],aes(x=Sp,y=value),group=1)+
    geom_point(data = Sfamily_ratio_C[Sfamily_ratio_C$Chr==i,],aes(x=Sp,y=value),shape=21,fill='#FFE38D')+
    geom_line(data = Sfamily_ratio_OC[Sfamily_ratio_OC$Chr==i,],aes(x=Sp,y=value),group=1)+
    geom_point(data = Sfamily_ratio_OC[Sfamily_ratio_OC$Chr==i,],aes(x=Sp,y=value),shape=21,fill='#718dcc')+
    labs (y=i,x='')+
    theme_bw() +
    # facet_grid(Chr~.)+
    coord_flip()
  ggsave(paste(i,'.pdf',sep=''),chr_plot,height = 7.5,width =2.5,path="E:/Bio_analysis/Weedyrice/pan_weedyrice/1207_statresult/Centro_region")
  
}
