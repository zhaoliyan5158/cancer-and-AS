#数据整理
library(tidyverse) 
library(broom)      # 用于整理模型输出的包
library(glmnet)     # 用于拟合广义线性模型和弹性网正则化的包
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
k <- cbind(fenK$lasso, k) 
states <- as.matrix(k)    # 将k转换为矩阵
#修改列名
names(k)[1] <- "fen.lasso"  # 如果您知道列的位置

library(caret)

library(shape)
display.progress = function (index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
}

##k<-read.csv("k.csv",row.names = 1)#不需要
set.seed(19991018) # 设置种子
control <- rfeControl(functions = rfFuncs, # 选择随机森林；详细可参考http://topepo.github.io/caret/recursive-feature-elimination.html#rfe
                      method = "LGOCV", # 选择交叉验证法
                      number = 10) # 10折交叉验证
tmp <- k
candidate.gene<-colnames(k)

results <- rfe(x = tmp[,-1], 
               y =as.factor(k$fen.lasso), 
               metric = "Accuracy",
               sizes = 1:(length(candidate.gene)-2), # 步长为1，速度较慢请耐心（约2小时）#-2要注意，是基因数-1，有一列为分组信息，所以-2
               rfeControl = control)

final.gene <- predictors(results) # 取出最终基因
write.table(final.gene,"output_selected features.txt",sep = "\t",row.names = F,col.names = F,quote = F)

accres <- results$results # 取出迭代结果
write.table(accres, "output_accuracy result.txt", sep = "\t", row.names = F,col.names = T,quote = F)


# 设置颜色
jco <- c("#2874C5","#EABF00")

# 图1：随机森林准确性图
pdf(file = "accuracy.pdf", width = 6, height = 4.5)
par(bty="o", mgp = c(2,0.5,0), mar = c(3.1,4.1,2.1,2.1),tcl=-.25,las = 1)
index <- which.max(accres$Accuracy) # 取出准确率最大时的索引（基因个数）
## 画圈圈
par(cex=1.5) # 设置全局文本大小为默认的1.5倍
plot(accres$Variables,
     accres$Accuracy,
     ylab = "",
     xlab = "Number of genes",
     col = "steelblue")
## 添加连线
lines(accres$Variables,accres$Accuracy,col = "steelblue")
## 定位最大值
points(index, accres[index,"Accuracy"],
       col = "steelblue",
       pch = 19,
       cex = 1.2)
## 补Y轴坐标（在plot时候写会和axis文字重叠）
#mtext("Accuracy (Repeated Cross-Validation)",side = 2,line = 2.5, las = 3)
mtext("Accuracy (Repeated Cross-Validation)", side = 2, line = 2.5, las = 3, cex = 1.5)
## 添加好看滴箭头
Arrows(x0 = index-7, x1 = index-2,
       y0 = accres[index,"Accuracy"], y1 = accres[index,"Accuracy"],
       arr.length = 0.2,
       lwd = 2,
       col = "black",
       arr.type = "triangle")
## 添加基因数目信息
text(x = index - 7,
     y = accres[index,"Accuracy"],
     labels = paste0("N=",index),
     pos = 2)




invisible(dev.off())

 