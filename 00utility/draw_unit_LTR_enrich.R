# library(ComplexHeatmap)
# library(circlize)
# freq_sp_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/1207_statresult/01_centro_soloLTR_freq.csv"
# freq_sp_df<-read.csv(freq_sp_file)
# rownames(freq_sp_df)<-freq_sp_df$X
# freq_sp_df<-freq_sp_df[,-1]

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggsci)

##以物种分类，富集unit分布
enrich_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/1207_statresult/02_centro_fullLTR_freq_ful_sp_long.txt"
# enrich_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/1207_statresult/02_centro_fullLTR_length_ful_sp_long.txt"
freq_sp_df<-read.table(enrich_file,header = 1)

phy_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/70_phylogenic.txt"
phy_df<-read.table(phy_file,col.names = c("sp","group"))
sp_group<-c(phy_df$group)
names(sp_group)<-c(phy_df$sp)

centro_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/70_pan_centro.bed"
centro_df<-read.table(centro_file,col.names = c("chr_sp","start","end",'no'))
centro_df <- separate(centro_df, chr_sp, into = c("chr", "sp"), sep = "_")
centro_df$length<-centro_df$end-centro_df$start
centro_df <- centro_df %>%
  group_by(sp) %>%
  summarize(mean_value = sum(length, na.rm = TRUE))
colnames(centro_df)<-c('sp','length')
centro_df$sp<-factor(centro_df$sp, levels = rev(names(sp_group)))
freq_sp_df$Sp<-factor(freq_sp_df$Sp, levels = rev(names(sp_group)))
# sp_plot<-ggplot(data=centro_df,aes(sp,length/48000))+
sp_plot<-ggplot(data=centro_df,aes(sp,length/3))+
          geom_bar(stat="identity", fill="lightgrey", width=0.7,linewidth=0.25)+
          # ggplot(data=freq_sp_df,aes(Sp,value,fill=name))+
          geom_bar(data=freq_sp_df[freq_sp_df$value!=0,],aes(Sp,value,fill=name),stat="identity",position="stack", 
                   width=0.7,linewidth=0.25,inherit.aes=FALSE,)+
          theme_classic()+
          guides(fill = guide_legend(ncol=1))+
          coord_flip()
ggsave('03_Sp_fullLTR_freq.pdf',sp_plot,path = dirname(enrich_file),height = 12,width = 7)


##以Chr分类，富集unit分布
# enrich_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/1207_statresult/02_centro_fullLTR_freq_ful_chr_long.txt"
enrich_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/1207_statresult/02_centro_fullLTR_length_ful_chr_long.txt"
freq_chr_df<-read.table(enrich_file,header = 1)
centro_df<-read.table(centro_file,col.names = c("chr_sp","start","end",'no'))
centro_df <- separate(centro_df, chr_sp, into = c("chr", "sp"), sep = "_")
centro_df$length<-centro_df$end-centro_df$start
centro_df <- centro_df %>%
  group_by(chr) %>%
  summarize(mean_value = sum(length, na.rm = TRUE))
colnames(centro_df)<-c('chr','length')
# chr_plot<-ggplot(data=centro_df,aes(chr,length/8400))+
chr_plot<-ggplot(data=centro_df,aes(chr,length/3))+
  geom_bar(stat="identity", fill="lightgrey", width=0.7,linewidth=0.25)+
  # ggplot(data=freq_sp_df,aes(Sp,value,fill=name))+
  geom_bar(data=freq_chr_df,aes(Sp,value,fill=name),stat="identity",position="stack",
           width=0.7,linewidth=0.25,inherit.aes=FALSE)+
  theme_classic()+
  theme(
    legend.position = "bottom",  # 将图例放置在底部
    legend.box = "horizontal" ,  # 将图例横向排列
    legend.text.align = 1,  # 图例文本水平对齐
    legend.title.align = 1,  # 图例标题旋转并放置在图标下方
    legend.text = element_text(angle = 90, vjust = -1)  # 图例标题旋转90度并向下对齐
  )+
  guides(fill = guide_legend(nrow=1))+
  coord_flip()
ggsave('03_Chr_fullLTR_length.pdf',chr_plot,path = dirname(enrich_file),height = 8,width = 7)


#补充绘图脚本
if (0) {
  ##根据unit插入类型（sata或者sd）分布数量
  freq_sp_type<-freq_sp_df[,c('Sp','name','statue','value')]
  freq_sp_type_df <- spread(freq_sp_type, key = "statue",
                            value = "value")
  freq_sp_type_df[is.na(freq_sp_type_df)]=0
  freq_sp_type_df<-freq_sp_type_df[,c('Sp','name','intact','sata')]
  pal<-pal_d3("category20")(length(unique(freq_sp_type_df$name)))
  ggplot(freq_sp_type_df,aes(x = intact,y=sata,color=name))+geom_point()+scale_color_manual(values =rev(pal))
  
}
