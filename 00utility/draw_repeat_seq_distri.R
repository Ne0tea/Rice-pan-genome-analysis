library(ComplexHeatmap)
library(circlize)
phy_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/70_phylogenic.txt"
phy_df<-read.table(phy_file,col.names = c("sp","group"))
sp_group<-c(phy_df$group)
names(sp_group)<-c(phy_df$sp)
print(sp_group)

length_sp_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/intact_stat/centro_region/parse_seq_step_tmp/length_sp_df.csv"
length_sp_df<-read.csv(length_sp_file)
rownames(length_sp_df)<-length_sp_df$sp
length_sp_df<-length_sp_df[,-1]

freq_sp_file="E:/Bio_analysis/Weedyrice/pan_weedyrice/intact_stat/centro_region/parse_seq_step_tmp/freq_sp_df.csv"
freq_sp_df<-read.csv(freq_sp_file)
rownames(freq_sp_df)<-freq_sp_df$sp
freq_sp_df<-freq_sp_df[,-1]

freq_col <- colorRamp2(
  breaks = c(0, 5, 10), 
  colors = c("#f1faee", "#a8dadc", "#457b9d")
)
length_col <- colorRamp2(
  breaks = c(0, 10000, 20000), 
  colors = c("#f1faee", "#a8dadc", "#457b9d")
)
for (i in seq(1:12)) {
  if (i<10) {
    chr<-paste('Chr0',i,sep = '')
  }else {
    chr<-paste('Chr',i,sep = '')
  }
  freq_sp_file=paste("E:/Bio_analysis/Weedyrice/pan_weedyrice/intact_stat/centro_region/repeat_anno_group_by_chr/length_",chr,"_sp_df.csv",sep = '')
  freq_sp_df<-read.csv(freq_sp_file)
  rownames(freq_sp_df)<-freq_sp_df$sp
  freq_sp_df<-freq_sp_df[,-1]
  draw_matrix<-as.matrix(freq_sp_df)
  inner_sp<-intersect(names(sp_group),rownames(draw_matrix))
  draw_matrix<-draw_matrix[inner_sp,]
  pdf(file=paste("E:/Bio_analysis/Weedyrice/pan_weedyrice/intact_stat/centro_region/repeat_anno_group_by_chr/length_",chr,"_sp_df.pdf",sep = ''),width = 20,height = 10)
  ht_list<-Heatmap(draw_matrix,col=length_col,
                   # width = unit(1, "bigpts")*ncol(draw_t), 
                   # height = unit(10, "bigpts")*nrow(draw_t),
                   cluster_columns=FALSE,
                   cluster_rows=FALSE,
                   row_names_gp=gpar(fontsize = 10),
                   column_names_gp=gpar(fontsize = 8),
                   # right_annotation = right_an,
                   # top_annotation = column_ha, 
                   # show_heatmap_legend=FALSE,
                   # show_column_names = FALSE,
                   rect_gp = gpar(col = "grey",lwd = 0.5),
                   row_names_side='right',
                   row_title_rot=0,
                   row_title_gp=gpar(fontsize=25),
                   column_title_rot=45,
                   column_title_gp = gpar(fontsize=10),
                   row_split = factor(unname(sp_group[inner_sp]),levels = unique(unname(sp_group[inner_sp]))),
                   cluster_column_slices = FALSE,
                   column_gap = unit(1, "mm")
  )
  draw(ht_list)
  print('Chr is done!')
  dev.off()
}

pdf(file=paste("E:/Bio_analysis/Weedyrice/pan_weedyrice/length_ful_sp_df.pdf",sep = ''),width = 20,height = 10)
###length 和 freq所使用的legend范围不一致
length_col <- colorRamp2(
  breaks = c(0, 10000, 20000), 
  colors = c("#f1faee", "#a8dadc", "#457b9d")
)
draw_matrix<-as.matrix(length_sp_df)
draw_matrix<-draw_matrix[names(sp_group),]
ht_list<-Heatmap(draw_matrix,col=length_col,
                 cluster_columns=FALSE,
                 cluster_rows=FALSE,
                 row_names_gp=gpar(fontsize = 10),
                 column_names_gp=gpar(fontsize = 8),
                 # right_annotation = right_an,
                 # top_annotation = column_ha, 
                 # show_heatmap_legend=FALSE,
                 # show_column_names = FALSE,
                 rect_gp = gpar(col = "grey",lwd = 0.5),
                 row_names_side='right',
                 # row_order=c(names(sp_group)),
                 row_title_rot=0,
                 row_title_gp=gpar(fontsize=25),
                 column_title_rot=45,
                 column_title_gp = gpar(fontsize=10),
                 row_split = factor(sp_group,levels = unique(sp_group)),
                 cluster_column_slices = FALSE,
                 column_gap = unit(1, "mm")
                 # use_raster = TRUE,
                 # raster_device = "png"
)
draw(ht_list)
dev.off()

pdf(file=paste("E:/Bio_analysis/Weedyrice/pan_weedyrice/freq_ful_sp_df.pdf",sep = ''),width = 20,height = 10)
###length 和 freq所使用的legend范围不一致
freq_col <- colorRamp2(
  breaks = c(0, 20, 40), 
  colors = c("#f1faee", "#a8dadc", "#457b9d")
)
draw_matrix<-as.matrix(freq_sp_df)
draw_matrix<-draw_matrix[names(sp_group),]
ht_list<-Heatmap(draw_matrix,col=freq_col,
                 # width = unit(1, "bigpts")*ncol(draw_t), 
                 # height = unit(10, "bigpts")*nrow(draw_t),
                 cluster_columns=FALSE,
                 cluster_rows=FALSE,
                 row_names_gp=gpar(fontsize = 10),
                 column_names_gp=gpar(fontsize = 8),
                 # right_annotation = right_an,
                 # top_annotation = column_ha, 
                 # show_heatmap_legend=FALSE,
                 # show_column_names = FALSE,
                 rect_gp = gpar(col = "grey",lwd = 0.5),
                 row_names_side='right',
                 row_order=c(names(sp_group)),
                 row_title_rot=0,
                 row_title_gp=gpar(fontsize=25),
                 column_title_rot=45,
                 column_title_gp = gpar(fontsize=10),
                 # column_order = sg_ord,
                 # bottom_annotation = length_an,
                 # column_title = " ",
                 # border = TRUE,
                 # heatmap_legend_param = list(
                 #   # at = c(-2, 0, 2),
                 #   # labels = c("low", "zero", "high"),
                 #   # title = paste("Hapstat",i),
                 #   # legend_height = unit(5, "cm"),
                 #   title_position = "leftcenter",nrow=1,
                 #   grid_height = unit(2, "cm"), 
                 #   grid_width = unit(2, "cm"),
                 #   # legend_height=unit(4, "cm"),
                 #   # legend_width=unit(2, "cm"),
                 #   title_gp = gpar(fontsize = 20,cex=3),
                 #   labels_gp = gpar(fontsize = 30),
                 #   gap=unit(2,"cm")
                 # ),
                 row_split = factor(sp_group,levels = unique(sp_group)),
                 cluster_column_slices = FALSE,
                 column_gap = unit(1, "mm")
                 # use_raster = TRUE,
                 # raster_device = "png"
)
draw(ht_list)
dev.off()