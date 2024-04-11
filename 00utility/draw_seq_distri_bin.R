#绘制百分比的堆叠图

library(ggplot2)
library(tidyverse)
library(cowplot)
library(dplyr)

depth_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/ful_sp_SZ-22.depth"
depth_df<-read.table(depth_file,check.names=FALSE)
colnames(depth_df)<-c('Chr','Start','End','value')
my_umap_Plot <- function(i){
  depth_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/ful_sp_SZ-22.depth"
  depth_R="E:/Bio_analysis/Weedyrice/pan_weedyrice/ful_sp_RIRE7.depth"
  depth_CRM="E:/Bio_analysis/Weedyrice/pan_weedyrice/ful_sp_CRM.depth"
  depth_df<-read.table(depth_file,check.names=FALSE)
  colnames(depth_df)<-c('Chr','Start','End','value')
  depth_df$type<-'SZ-22'
  
  depth_Rf<-read.table(depth_R,check.names=FALSE)
  colnames(depth_Rf)<-c('Chr','Start','End','value')
  depth_Rf$type<-'RIRE7'
  
  depth_CRMf<-read.table(depth_CRM,check.names=FALSE)
  colnames(depth_CRMf)<-c('Chr','Start','End','value')
  depth_CRMf$type<-'CRM'
  
  depth_df<-rbind(depth_df,depth_Rf,depth_CRMf)
  draw_data<-depth_df[grep(i, depth_df$Chr),]
  draw_data<-separate(data = draw_data, col = Chr, into = c("Chr", "loc"), sep = "_")
  draw_data$loc<-factor(draw_data$loc,levels = unique(draw_data$loc))
  selected_groups <- draw_data[c(1, which(round(seq(1,nrow(draw_data))/200, 5) %in% c(0.20000,0.40000, 0.60000,0.80000)), nrow(draw_data)),'loc']
  ggplot(data=draw_data,aes(loc,value,fill=type))+
    geom_bar(stat="identity", width=1,linewidth=0.5,position = 'stack')+
    scale_x_discrete(name = i,labels = c('0','75','80','160','165','240'),breaks=selected_groups) +  # 将离散型 X 轴设置为空标签
    # scale_x_continuous(breaks = 1:length(draw_data$Chr), labels = draw_data$Chr) +  # 设置连续型 X 轴
    scale_fill_manual(values = c('SZ-22'='#bc8a57','RIRE7'='#b7c7a0','CRM'='#88ada6'))+
    # facet_zoom(xlim = c(paste(i,'_40',sep = ''), paste(i,'_80',sep = '')))+
    ylab('Count/100bp')+
    theme_classic()+
    theme(legend.position = 'none')
}
for (i in unique(separate(data = depth_df, col = Chr, into = c("Chr", "xx"), sep = "_")[1])){
  umap_list <- lapply(i,my_umap_Plot)
}
draw_data<-depth_df[grep(i, depth_df$Chr),]
xx<-ggplot(data=draw_data,aes(Chr,value,fill=type))+
  geom_bar(stat="identity", width=1,linewidth=0.5,position = 'stack')+
  scale_x_discrete(name = i,labels = c('0','75','80','160','165','240'),breaks=selected_groups) +  # 将离散型 X 轴设置为空标签
  # scale_x_continuous(breaks = 1:length(draw_data$Chr), labels = draw_data$Chr) +  # 设置连续型 X 轴
  scale_fill_manual(values = c('SZ-22'='#bc8a57','RIRE7'='#b7c7a0','CRM'='#88ada6'))+
  # facet_zoom(xlim = c(paste(i,'_40',sep = ''), paste(i,'_80',sep = '')))+
  theme_classic()+
  theme(legend.position = 'none')
plot_grid(plotlist = umap_list, align = "h",nrow = 3)


depth_df<-separate(data = depth_df, col = Chr, into = c("Chr", "loc"), sep = "_")
depth_df$loc<-factor(depth_df$loc,levels = unique(depth_df$loc))
selected_groups <- depth_df[c(1, which(round(seq(1,nrow(depth_df))/200, 5) %in% c(0.20000,0.40000, 0.60000,0.80000)), nrow(depth_df)),'loc']
ggplot(data=depth_df,aes(loc,value,fill=type))+
  geom_bar(stat="identity", width=1,linewidth=0.5,position = 'stack')+
  scale_x_discrete(labels = c('0','20','40','60','80','100'),breaks=selected_groups) +  # 将离散型 X 轴设置为空标签
  scale_fill_manual(values = c('SZ-22'='#bc8a57','RIRE7'='#b7c7a0','CRM'='#88ada6'))+
  xlab(NULL)+ylab('Count/100bp')+
  guides(fill = guide_legend(title = ''))+
  theme_classic()+
  theme(legend.position = 'top')

#计算12chr平均值
new_depth_df <- depth_df %>%
  group_by(loc,type) %>%
  summarize(mvalue = mean(value))

#计算type占bin的比值
new_depth_df <- new_depth_df %>%
  group_by(loc) %>%
  mutate(E = sum(mvalue)) %>%
  ungroup() %>%
  mutate(proportion = ifelse(is.nan(mvalue / E), 0, mvalue / E))

new_depth_df$loc <- as.numeric(new_depth_df$loc)
selected_groups <- c(1, which(round(unique(new_depth_df['loc'])/200, 5)$loc %in% c(0.20000,0.40000, 0.60000,0.80000)),240)
#绘制每个bin中type所占比例（mean（12chr））
ggplot(new_depth_df, aes(x = loc, y = proportion, color = type), group = 1) +
  geom_line() +
  geom_point() +
  # scale_x_discrete(labels = c('0','20','40','60','80','100'),breaks=selected_groups) +
  scale_color_manual(values = c('SZ-22'='#bc8a57','RIRE7'='#b7c7a0','CRM'='#88ada6'))+
  labs(x = "", y = "Count proportion", color = "") +
  theme_classic()+theme(legend.position = 'top')
