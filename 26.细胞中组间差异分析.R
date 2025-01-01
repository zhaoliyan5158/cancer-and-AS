#组间单细胞差异分析
rm(list = ls(all.names = TRUE))
library(dplyr)
library(Seurat)
library(patchwork)
load(file = "step3-1output.Rdata")
scRNA=mergelayers
names(scRNA@meta.data)
unique(scRNA$group)
#创建新的细胞标识符：这段代码首先合并celltype和group两个列来创建一个新的列celltype.group，用于区分不同的细胞类型和分组。
scRNA$celltype.group <- paste(scRNA$celltype, scRNA$group, sep = "_")
Idents(scRNA) <- "celltype.group"
#寻找差异表达的标记基因：使用FindMarkers函数寻找在MAC_0 week和MAC_8 weeks两个分组间差异表达的基因。
##ident.1是病态样本要比较组(Disease)，ident.2是被比较组，健康样本(Control)
CELLDEG <- FindMarkers(scRNA, ident.1 = 'MAC_8 weeks', ident.2 = 'MAC_0 week', verbose = FALSE, test.use = 'wilcox', min.pct = 0.1)
#获取细胞类型：
cellForDeg <- levels(scRNA$cellType)
#for(i in 1:length(cellForDeg)){
  CELLDEG <- FindMarkers(scRNA, ident.1 = paste0(cellForDeg[i], '_GC'), ident.2 = paste0(cellForDeg[i], '_Normal'), verbose = FALSE)
  write.csv(CELLDEG, paste0(cellForDeg[i], ".csv"))}
FeaturePlot(pbmc, features = "Daxx")
DotPlot(pbmc, features = "Daxx") + RotatedAxis()

library(Seurat)
library(ggplot2)

# 假设你的Seurat对象名为pbmc
# 提取Daxx基因表达数据和细胞群体信息
# 假设`celltype`和`group`是两列分别代表细胞类型和分组信息
data_for_plot <- FetchData(scRNA, vars = c("Irf7", "celltype", "group"))



ggplot(data_for_plot, aes(x = celltype, y = Irf7, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("0 week" = "blue", "8 weeks" = "orange", "16 weeks" = "yellow", "26 weeks" = "red")) +
  theme_bw() +
  labs(title = "Expression of IRF7 across Cell Types", x = "Cell Type", y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0.1, max(data_for_plot$Irf7, na.rm = TRUE))  # 设置Y轴的范围




ggplot(data_for_plot, aes(x = group, y = Irf7, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("0 week" = "blue", "8 weeks" = "orange", "16 weeks" = "yellow", "26 weeks" = "red")) +
  theme_bw() +
  labs(title = "Expression of IRF7 across Cell Types", x = "Cell Type", y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0.1, max(data_for_plot$Irf7, na.rm = TRUE))+ # 设置Y轴的范围
  facet_wrap(~ celltype, scales = "free_y")