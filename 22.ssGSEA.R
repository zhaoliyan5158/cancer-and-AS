#ssGSEA
rm(list = ls()) 
library(openxlsx)
library(GSVA)
library(msigdbr)
library(corrplot)
library(dplyr)
data <- read.csv("exp.csv", row.names=1)
# 从数据框中删除第一列
data <- data [,-1]
fenK <- read.xlsx("samples.xlsx")   
#转换组别信息group
# 转换 group 列中的 0 和 1 为 "SCM" 和 "Healthy"
group <- ifelse(fenK$group == 1, "control", "AS")

#获取 Hallmark 基因集
# 加载 msigdbr 包中的 Hallmark 基因集
msigdb_data <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_gene_sets <- split(msigdb_data$gene_symbol, msigdb_data$gs_name)
# 继续使用传统的调用方式
ssgsea_scores <- gsva(as.matrix(data), hallmark_gene_sets, method = "ssgsea", kcdf = "Gaussian")

# 转置 ssGSEA 得分
ssgsea_scores <- t(ssgsea_scores)


#计算标志性基因表达与 ssGSEA 得分的相关性
#基因集
signature_genes <- c("TNF", "FHOD1", "IRF7","ZSWIM3")

# 提取标志性基因表达数据
signature_expression <- data[signature_genes, ]

# 计算标志性基因表达与 ssGSEA 得分的相关性矩阵
cor_matrix <- cor(t(signature_expression), ssgsea_scores)
##去掉HALLMARK前缀
# 获取当前列名和行名
original_colnames <- colnames(cor_matrix)
original_rownames <- rownames(cor_matrix)

# 去掉 "HALLMARK_" 前缀
new_colnames <- gsub("HALLMARK_", "", original_colnames)
new_rownames <- gsub("HALLMARK_", "", original_rownames)

# 赋值回矩阵
colnames(cor_matrix) <- new_colnames
rownames(cor_matrix) <- new_rownames


#可视化相关性矩阵
# 使用 corrplot 包可视化相关性矩阵
corrplot(cor_matrix, method = "color", type = "full", tl.col = "black", tl.cex = 0.8, 
         addCoef.col = "black", number.cex = 0.7, col = colorRampPalette(c("blue", "white", "red"))(200))







