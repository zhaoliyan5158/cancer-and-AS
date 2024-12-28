#数据整理
library(tidyverse) 
library(broom)      # 用于整理模型输出的包
library(glmnet)     # 用于拟合广义线性模型和弹性网正则化的包
k <- read.csv("common_genes_exp.csv", row.names = 1)  # 读取CSV文件，并设置第一列为行名
k <- as.data.frame(t(k))                                     # 转置k并转换为数据框
library(openxlsx)
fenK <- read.xlsx("samples.xlsx")                   # 读取另一个CSV文件，设置第一列为行名
# 将第一列设置为行名
rownames(fenK) <- fenK[[1]]
fenK <- fenK[-1]
k <- k[row.names(fenK), ]                                    # 从k中筛选出fenK中存在的行
k <- cbind(fenK$lasso, k) 
a<-k
library(tidyverse)
library(glmnet)
source('msvmRFE.R')   #文件夹内自带
library(VennDiagram)
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)
train<-a
input <- train

#采用五折交叉验证 (k-fold crossValidation）
svmRFE(input, k = 5, halve.above = 100) #分割数据，分配随机数
nfold = 5
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100) #特征选择
top.features = WriteFeatures(results, input, save=F) #查看主要变量
head(top.features)

#把SVM-REF找到的特征保存到文件
write.csv(top.features,"feature_svm.csv")

# 运行时间主要取决于选择变量的个数，一般的电脑还是不要选择太多变量
# 选前48个变量进行SVM模型构建，体验一下

featsweep = lapply(1:214, FeatSweep.wrap, results, input) #214个变量就修改为214,耗时13小时
featsweep

load("featsweep.RData")
# 选前300个变量进行SVM模型构建，然后导入已经运行好的结果
#featsweep = lapply(1:300, FeatSweep.wrap, results, input) #300个变量
save(featsweep,file = "featsweep.RData")


#画图
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

#dev.new(width=4, height=4, bg='white')
#pdf("B_svm-error.pdf",width = 5,height = 5)
PlotErrors(errors, no.info=no.info) #查看错误率
par(cex=1.5) # 设置全局文本大小为默认的1.5倍
PlotErrors(errors, no.info=TRUE)
#dev.off()

#dev.new(width=4, height=4, bg='white')
#pdf("B_svm-accuracy.pdf",width = 5,height = 5)
Plotaccuracy(1-errors,no.info=no.info) #查看准确率

#dev.off()

# 图中红色圆圈所在的位置，即错误率最低点
which.min(errors) 
top<-top.features[1:which.min(errors), "FeatureName"]
write.csv(top,"top.csv")
