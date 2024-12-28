
#数据准备：行名为基因名，列名为样本名
rm(list = ls())  
load(file = "step1output.Rdata")
load(file = "step2output.Rdata")
library(dplyr)
library(WGCNA)
###数据准备####

#1.加probe_id列，把行名变成一列
class(exp)
exp_df <- as.data.frame(exp)
exp_df <- mutate(exp_df,probe_id=rownames(exp_df))
head(exp_df)
#2.加上探针注释
#有多个探针对应同一个基因，怎么办？1.取最大值；2.平均值；3.随机选取，这里用随机去重
ids2= ids2[!duplicated(ids2$GENE_SYMBOL),]
#更改列名
colnames(ids2) <- c("probe_id", "symbol")
ids= ids2
#其他去重方式在zz.去重.R
exp_symbol <- inner_join(exp_df,ids,by="probe_id")
head(exp_symbol)
nrow(exp_symbol)
# 首先将 'symbol' 列的值设置为行名
row.names(exp_symbol) <- exp_symbol$symbol
# 然后从数据框中移除 'symbol' 列
exp_symbol$symbol <- NULL  
# 然后从数据框中移除 'probe_id' 列
exp_symbol$probe_id <- NULL  
#数据框转换为矩阵
exp_symbol<- as.matrix(exp_symbol)
mydata<-exp_symbol#备份数据
datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
#筛选筛选方差前25%的基因：剔除那些随性状改变变化不大的基因
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[4])]
datExpr1<-data.matrix(expro.upper)
#判断是否有不好的sample或gene
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
#如果这里返回的结果是TRUE，说明所有基因都通过了检查。

#如果你用全部基因作为输入，很有可能返回FALSE，说明存在不好的基因或sample。

#下面的代码就会去除那些不好的基因或sample。
#去除不好的sample或gene
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
#判断是否有离群样本
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar=c(5, 4, 4, 2) + 0.5) # 调整边距
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2,
     cex.axis = 2, cex.main = 2)
#如果你的数据有离群样本需要去除，就运行下面这段代码。
#去除离群样本
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2, cex.axis = 2, cex.main = 2) +
  #想用哪里切，就把“h = 110”和“cutHeight = 110”中换成你的cutoff
  abline(h = 80, col = "red") 


clust = cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
keepSamples = (clust==1)#保留非离群样本，2是保留离群样本
datExpr = datExpr1[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
dim(datExpr)

datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#选择构建网络的合适阈值
#通过这步计算，找出scale free topology modle fit接近0.9的最小power（soft threshold），用于下一步构建网络。
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold1.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 1.5
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
#构建网络，找出gene module
net = blockwiseModules(datExpr, power = 14,#选择变得平缓的数值
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)
table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.pdf",width = 10, height = 5)
par(cex.axis = 1.2, cex.lab = 1.2) # 设置坐标轴标签和图例的字体大小为默认大小的1.2倍
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
#把gene module输出到文件
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"csv",sep = "."),quote = F)
}
#表型与模块的相关性
samples=read.csv('samples.csv',row.names = 1)
#samples=read.xlsx("samples.xlsx", sheet = 1)
moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.pdf",width = 6, height = 6)
par(mar=c(5, 4, 4, 2) + 2) # 调整边距
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 1,  yColorWidth=0.02, 
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 1, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()


