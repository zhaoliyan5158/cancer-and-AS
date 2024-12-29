#数据整理
# 加载rms包
library(rms)
k <- read.csv("common_genes_exp.csv", row.names = 1)  # 读取CSV文件，并设置第一列为行名
common_genes <- read.csv("common_genes.csv", header = TRUE)
common_genes_vector <- common_genes$Symbol  # 假设基因名列名为 GeneName
#使用 common_genes_vector 从 k 中筛选出对应的基因：
filtered_data <- k[rownames(k) %in% common_genes_vector,]
m <- as.data.frame(t(filtered_data))                                     # 转置k并转换为数据框
library(openxlsx)
fenK <- read.xlsx("samples.xlsx")                   # 读取另一个CSV文件，设置第一列为行名
# 将第一列设置为行名
rownames(fenK) <- fenK[[1]]
fenK <- fenK[-1]
m <- m[row.names(fenK), ]                                    # 从k中筛选出fenK中存在的行
m <- cbind(fenK$lasso, m) 
states <- as.matrix(m)    # 将k转换为矩阵
#修改列名
names(m)[1] <- "status"  # 如果您知道列的位置
rt<-m
###
# 加载所需的库
library(rms)
library(survival)

# 设置列名，这里可能是对特定列进行操作，列名根据实际数据确定
colnames(rt)

# 设置数据集的距离计算方式，可能是为了某些模型或函数设置
dist <- datadist(rt)
options(datadist = "dist")

# 逻辑回归模型拟合，使用了 rms 包的 lrm 函数
#fit <- lrm(status~FHOD1+IRF7+TNF+ZSWIM3, data = rt, x = T, y = T)#四个hub基因,拟合不出

fit <- lrm(status~FHOD1+IRF7+TNF, data = rt, x = T, y = T)#3个hub基因
fit

# 设置全局图形参数，增加字体大小和加粗
par(cex = 1.5, font = 2)  # cex控制字体大小，font=2表示加粗

# 创建nomogram
nom <- nomogram(fit, fun = plogis, 
                fun.at = c(0.001, 0.1, 0.4, 0.7, 0.95, 0.999, 1),
                lp = F, funlabel = "Linear Predictor")

# 绘制nomogram
plot(nom)

# 使用 rms 包的 c-index 计算模型预测的准确度
cindex <- rcorrcens(status~predict(fit), data = rt)

# 加载 PROC 包，通常用于 ROC 曲线的绘制
library(pROC)

# 初始化图形布局
par(mfrow=c(ceiling(length(genes)/4), 2))  # 根据基因数量调整

# 使用ROC函数计算ROC曲线
gfit <- roc(status~predict(fit), data = rt,ci = TRUE)#ci是显示95%置信区间


# 绘制ROC曲线，不自动添加坐标轴标签
library(pROC)

# 设置图形布局为1行2列
par(mfrow=c(1, 2))
plot(gfit,
     print.auc=TRUE,  # 在图中显示AUC值
     main="Nomogram in training set",  # 图的主标题
     col="#FF0000",  # ROC曲线颜色，确保只有一个颜色设置
     print.thres=FALSE,  # 不打印阈值
     identity.col="black",  # 对角线（完全随机线）颜色
     identity.lty=1,  # 对角线的线型
     identity.lwd=1,  # 对角线的线宽
     legacy.axes=TRUE,  # 使用传统坐标轴样式
     xlab="", ylab="",
     cex.axis=1.2,  # 调整坐标轴刻度标签的字体大小
     font.axis=2)  # 清空自动生成的坐标轴标签

# 手动添加X轴和Y轴标签，字体大小加大并加粗
mtext("1-Specificities", side=1, line=3.5, cex=1.5, font=2)
mtext("Sensitivity", side=2, line=2.5, cex=1.5, font=2)


# 校准曲线的计算
cal <- calibrate(fit, method="boot", B=1000)

# 绘制校准曲线
plot(cal,
     xlab="Nomogram-predicted probability of nonadherence",
     ylab="Actual diagnosed nonadherence (proportion)",
     sub=F)
#每个基因的ROC曲线怎么绘制
# 获取基因列表
genes <- colnames(rt)[-1]  # 假设第一列是状态列


library(pROC)

# 设置图形布局为1行2列
par(mfrow=c(1, 2))

for (gene in genes) {
  # 计算ROC曲线，并开启置信区间的计算
  roc_obj <- roc(rt$status, rt[[gene]], ci=TRUE, boot.n=1000)  # 启用置信区间计算
  
  # 绘制ROC曲线
  plot(roc_obj, print.auc=TRUE,
       main=paste("ROC for", gene), 
       col="#1c61b6", xlab="", ylab="",
       print.thres=FALSE, print.thres.col="black",  # 阈值标记的颜色
       identity.col="#1c61b6",  # 对角线（完全随机线）颜色
       identity.lty=1,  # 对角线的线型
       identity.lwd=1,  # 对角线的线宽
       legacy.axes=TRUE,  # 使用传统坐标轴样式
       cex.axis=1.2,  # 调整坐标轴刻度标签的字体大小
       font.axis=2)  # 调整坐标轴刻度标签的字体样式为加粗
  
  # 调整X轴和Y轴标签的字体大小并加粗
  mtext("1-Specificities", side=1, line=3.5, cex=1.5, font=2)
  mtext("Sensitivity", side=2, line=2.5, cex=1.5, font=2)
}


