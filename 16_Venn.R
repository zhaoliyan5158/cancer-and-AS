#差异基因
rm(list = ls()) 

library(ggplot2)
if (!require(ggvenn)) install.packages("ggvenn")
library(ggvenn)
library(openxlsx)

data1 <- read.csv("lasso.csv")
data2 <- read.csv("random forest.csv")
data3 <- read.csv("svm.csv")


#韦恩图
library(VennDiagram)
# 转换为列表
list_of_sets <- list(
  symbols_1 <- data1$symbol,
  symbols_2<- data2$symbol,
  symbols_3 <- data3$symbol
)

library(grid)
par(cex=1.5)
venn_plot <- venn.diagram(
  x = list_of_sets,
  filename = NULL,  # 设置为NULL以便在R环境中直接显示
  category.names = c("LASSO", "Random Forest", "SVM-REF"),
  fill = c("#299D8F", "#E9C46A", "#D87659"),  # 填充颜色
  col = "transparent",  # 边框颜色设置为透明
  cex = 1.5,  # 调整标签的字体大小
  fontface = "bold",  # 字体风格
  cat.col = c("#299D8F", "#E9C46A", "#D87659"),  # 类别名称的颜色与填充颜色相同
  cat.cex = 1.2,  # 类别名称的字体大小
  cat.pos = c(-20, 20, 135),  # 调整类别名称的位置（角度）
  cat.dist = 0.05,  # 类别名称与Venn图距离
  label.col = "black",  # 数字标签的颜色
  margin = 0.1  # 图形边缘留白
  )
grid.draw(venn_plot)




# 使用Reduce()和intersect()找出所有数据集的共同基因
common_genes <- Reduce(intersect, list(data1$symbol, data3$symbol, data2$symbol))
# 将向量转换为数据框
common_genes_df <- data.frame(Symbol = common_genes)
# 保存数据框为CSV文件
write.csv(common_genes_df, "common_genes.csv", row.names = FALSE)














