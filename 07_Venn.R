#deg1=AAA与control；deg2=AOD与control
rm(list = ls()) 
load(file = "gsea_outputAOD.Rdata")
load(file = "gsea_outputAAA.Rdata")
library(ggplot2)
if (!require(ggvenn)) install.packages("ggvenn")
library(ggvenn)
# 提取上调和下调基因：AAA中上调45，下调186共231；AOD：上调90，下调1039，共1129
up_genes1 <- deg1$symbol[deg1$change == "up"]
down_genes1 <- deg1$symbol[deg1$change == "down"]
up_genes2 <- deg2$symbol[deg2$change == "up"]
down_genes2 <- deg2$symbol[deg2$change == "down"]
# 组合基因集
gene_sets <- list(
  "Up in Set 1" = up_genes1,
  "Down in Set 1" = down_genes1,
  "Up in Set 2" = up_genes2,
  "Down in Set 2" = down_genes2
)
# 绘制韦恩图
ggplot() +
  ggvenn(gene_sets, show_percentage = TRUE) +
  theme_minimal()  # 使用简洁的主题
# 自定义颜色和其他选项
ggplot() +
  ggvenn(gene_sets, show_percentage = F) +
  scale_fill_manual(values = c("#9CD7C8", "#9ACBCD", "#F592B3", "#F9918C")) +
  theme(
    axis.text = element_blank(),    # 移除坐标轴文字
    axis.ticks = element_blank(),   # 移除坐标轴刻度
    panel.grid = element_blank(),   # 移除网格线
    plot.background = element_blank() # 移除背景
  )+annotate("text", x = 2, y = 1, label = "Down in AAA", size = 5, color = "black")
