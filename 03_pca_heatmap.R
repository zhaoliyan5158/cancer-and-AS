#生信技能树 31-GEO
rm(list = ls())  
load(file = "step1output.Rdata")
load(file = "step2output.Rdata")

# 1.PCA 图----
#按照PCA示例数据来调整自己的数据，调整为行→样本，列→基因
  #t（）转置函数，不论之前数据是什么类型，转置后都变为矩阵
dat=as.data.frame(t(exp))
library(FactoMineR)
library(factoextra) 
dat.pca <- PCA(dat, graph = FALSE)
#这个函数不能用于有NA或者无穷的数值
sum(is.na(dat))
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = Group, # color by groups
                         palette = c("#00AFBB", "#E7B800","#00688b","#00BFFF"),
                         addEllipses = TRUE, # Concentration ellipses
                         legend.title = "Groups"
)+theme(
  legend.text = element_text(size = 12),    # 图例文字大小
  legend.title = element_text(size = 14),  # 图例标题大小
  axis.text = element_text(size = 10),     # 坐标轴刻度文字大小
  axis.title = element_text(size = 12),    # 坐标轴标题大小
  plot.title = element_text(size = 14)     # 图标题大小
)
pca_plot
ggsave(plot = pca_plot,filename = paste0(gse_number,"_PCA.png"))
save(pca_plot,file = "pca_plot.Rdata")

# 2.top 1000 sd 热图---- 
#调出表达矩阵中方差较大的1000个基因名挑出来
cg=names(tail(sort(apply(exp,1,sd)),1000))
n=exp[cg,]

# 直接画热图，对比不鲜明
# 热图的调整 31-GEO-45分
library(pheatmap)
annotation_col=data.frame(group=Group)
rownames(annotation_col)=colnames(n) 
pheatmap(n,
         show_colnames =F,
         show_rownames = F,
         annotation_col=annotation_col
)

# 按行标准化scale：只关注一个基因在不同分组之间的表达变化，不关心一个基因与另一个基因变化的趋势
#breaks = seq(-3,3,length.out = 100) breaks 设置颜色范围，排除极端值的影响
pheatmap(n,
         show_colnames =F,
         show_rownames = F,
         annotation_col=annotation_col,
         scale = "row",
         breaks = seq(-3,3,length.out = 100)
         ) 
dev.off()

# 关于scale的进一步学习：zz.scale.R