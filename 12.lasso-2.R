library(tidyverse) 
library(broom)      # 用于整理模型输出的包
library(glmnet)     # 用于拟合广义线性模型和弹性网正则化的包
set.seed(123)
k <- read.csv("common_genes_exp.csv", row.names = 1)  # 读取CSV文件，并设置第一列为行名
#geneK <- read.csv("Treg 48.csv")                             # 读取另一个CSV文件
#k <- k[row.names(k) %in% geneK$gene, ]                       # 从k中筛选出geneK中存在的行
k <- as.data.frame(t(k))                                     # 转置k并转换为数据框
library(openxlsx)
fenK <- read.xlsx("samples.xlsx")                   # 读取另一个CSV文件，设置第一列为行名
# 将第一列设置为行名
rownames(fenK) <- fenK[[1]]
fenK <- fenK[-1]
k <- k[row.names(fenK), ]                                    # 从k中筛选出fenK中存在的行
k <- cbind(fenK$lasso, k)                                    # 将fenK中的lasso列绑定到k的前面
states <- as.matrix(k)                                       # 将k转换为矩阵
x <- states[, -1]                                            # x是去掉第一列后的所有列(自变量)
y <- states[, 1]     


fit=glmnet(x, y, family = "binomial", alpha=1)

cvfit=cv.glmnet(x, y, family="binomial", alpha=1, type.measure='deviance', nfolds = 10)
plot(cvfit, cex.axis = 1.5, cex.lab = 1.5) # 调整横坐标和纵坐标字体大小，以及Y轴标签字体大小

dev.off()
pdf(file="cvfit-1.pdf",width=6,height=5.5)


#????ɸѡ??????????
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)


cvfit$lambda.min#查看数值，后续要用
c(cvfit$lambda.min, cvfit$lambda.1se)
lasso <- glmnet(x, y, family="binomial", alpha=1, nlambda = 100)
coef(lasso, s=c(0.004095378,0.027578895))#两个数字记得修改为上面出现的数字
print(class(lasso))
plot(lasso, label=TRUE, cex.axis = 1.5, cex.lab = 1.5, cex=2)



plot(lasso,xvar="lambda", label=T,cex.axis = 1.5, cex.lab = 1.5, cex=2)

lasso.coef <- predict(lasso, s=0.4, type="coefficients")#回归系数
plot(lasso, xvar="dev", label=T)#解释偏差和回归系数关系
lasso.y <- predict(lasso, newx=x, type="response", s=0.4)#拟合
plot(lasso.y, y, xlab="Predicted", ylab="Actual", main="Lasso Regression")

