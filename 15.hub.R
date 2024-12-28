#数据整理

library(corrplot)
library(RColorBrewer)
k <- read.csv("common_genes_exp.csv", row.names = 1)  # 读取CSV文件，并设置第一列为行名
common_genes <- read.csv("common_genes.csv", header = TRUE)
common_genes_vector <- common_genes$Symbol  # 假设基因名列名为 GeneName
#使用 common_genes_vector 从 k 中筛选出对应的基因：
filtered_data <- k[rownames(k) %in% common_genes_vector,]
#使用 t() 函数转置数据：列名为基因名
transposed_data <- t(filtered_data)
#计算筛选后数据的相关矩阵：

cor_matrix <- cor(transposed_data)

#使用 corrplot 绘制热图，确保热图颜色和您提供的图像类似：
corrplot(cor_matrix, method = "color", type = "upper", order = "hclust",
         tl.col = "black", tl.srt = 45, # 标签颜色和角度
         col = colorRampPalette(c("#6C96CC", "white", "#C92321"))(200), # 颜色范围
         addCoef.col = "black") # 显示相关系数
#另一种画法
corrplot(corr = cor_matrix, 
         type = "upper",
         method = 'pie',
         col = colorRampPalette(colors = brewer.pal(11, "Spectral"))(200),
         outline = TRUE,
         tl.col = "black")
#另一种画法

corrplot(corr = cor_matrix, 
         type = "upper",
         method = 'pie',
         col = colorRampPalette(colors = brewer.pal(11, "Spectral"))(200),
         outline = TRUE,
         tl.col = "black",
         addCoef.col = "black",  # 显示相关系数的颜色
         number.cex = 0.7,       # 调整相关系数字体大小
         addCoefasPercent = FALSE, # 显示为数值而非百分比
         diag = T)           # 是否在对角线上显示

