#31-GEO-56分
rm(list = ls()) 
load(file = "step2output.Rdata")
#差异分析，用limma包来做
#需要表达矩阵和Group，不需要改
library(limma)
#model.matrix根据Group构建模型矩阵，以下是两组的代码
design=model.matrix(~Group)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)

#deg：AS vs 对照组


#1.加probe_id列，把行名变成一列
library(dplyr)
deg <- mutate(deg,probe_id=rownames(deg))
head(deg)
#2.加上探针注释
#有多个探针对应同一个基因，怎么办？1.取最大值；2.平均值；3.随机选取，这里用随机去重
ids2= ids2[!duplicated(ids2$GENE_SYMBOL),]
#更改列名
colnames(ids2) <- c("probe_id", "symbol")
#其他去重方式在zz.去重.R
deg <- inner_join(deg,ids2,by="probe_id")
head(deg)
nrow(deg)

#3.加change列,标记上下调基因
logFC_t=0.5
P.Value_t = 0.05
k1 = (deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value < P.Value_t)&(deg$logFC > logFC_t)
deg <- mutate(deg,change = ifelse(k1,"down",ifelse(k2,"up","stable")))
table(deg$change)
library(openxlsx)
write.csv(deg, "Differentially expressed genes.csv")

#4.加ENTREZID列，用于富集分析（symbol转entrezid，然后inner_join）
library(clusterProfiler)
library(org.Hs.eg.db)
s2e1 <- bitr(deg$symbol, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)#人类
#其他物种http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
dim(deg)
deg<- inner_join(deg,s2e1,by=c("symbol"="SYMBOL"))
dim(deg)
length(unique(deg$symbol))
save(Group,deg,logFC_t,P.Value_t,gse_number,file = "step4output1.Rdata")
# 提取change列为'up'的基因名
genes_up <- deg %>% filter(change == "up") %>% dplyr::select(symbol)

# 提取change列为'down'的基因名
genes_down <- deg %>% filter(change == "down") %>% dplyr::select(symbol)
library(openxlsx)
write.xlsx(genes_up, "Up_Regulated_Genes_GSE38713.xlsx")
write.xlsx(genes_down, "Down_Regulated_Genes_GSE38713.xlsx")
