#差异基因
rm(list = ls()) 

library(ggplot2)
if (!require(ggvenn)) install.packages("ggvenn")
library(ggvenn)
library(openxlsx)

data1 <- read.csv("cancer.csv")
data2 <- read.csv("tuquoise-brown.csv")
data3 <- read.xlsx("Up_Regulated_Genes_GSE100927.xlsx")
# 使用unique函数去除重复的基因名
data1 <- unique(data1$symbol)
data1 <-as.data.frame(data1)
colnames(data1)[1] <- "symbol"

#韦恩图
library(VennDiagram)
# 转换为列表
list_of_sets <- list(
  symbols_genes_tumor <- data1$symbol,
  symbols_genes_up_GSE100927 <- data3$symbol,
  symbols_genes_WGCNA <- data2$symbol
)

library(grid)
# 绘制韦恩图
venn_plot <- venn.diagram(
  x = list_of_sets,
  filename = NULL,
  category.names = c("Cancer_drivers", "AS_DEGs", "AS_WGCNA"),
  fill = c("#C86193", "#00BA38", "#619CFF"), # 填充颜色
  col = "transparent", # 将所有边框颜色设置为透明
  cat.cex = 2, # 调整类别名称的字体大小，1为默认大小，小于1则字体变小，大于1则字体变大
  cex = 2 # 同时调整数字的字体大小
)
grid.draw(venn_plot)














# 使用Reduce()和intersect()找出所有数据集的共同基因
common_genes <- Reduce(intersect, list(data1$symbol, data3$symbol, data2$symbol))
# 将向量转换为数据框
common_genes_df <- data.frame(Symbol = common_genes)
# 保存数据框为CSV文件
write.csv(common_genes_df, "common_genes.csv", row.names = FALSE)














