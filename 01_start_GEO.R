#数据下载#生信技能树 30-GEO
rm(list = ls())
library(GEOquery)
gse_number = "GSE100927"
#下载并读取getGEO函数
eSet <- getGEO(gse_number,    
               destdir = '.', 
               getGPL = F)
class(eSet)
length(eSet)
#把列表变成对象（ExpressionSet），取子集
eSet = eSet[[1]] 
#(1)提取表达矩阵exp:作者研发的一个函数，如下
#也等于eSet@assayData[["exprs"]]
exp <- exprs(eSet)
exp[1:4,1:4]
#一定判断数据有没有取过log，没有→用下面函数取log
#判断数据有没有取log：一般取过后数据在0-15之间。+1是为了避免0取log后值无穷小
#exp = log2(exp+1)
#基础包画箱线图，看看数据的质量，看纵坐标，有没有异常样本
boxplot(exp)
#判断需不需要标准化
#exp = limma::normalizeBetweenArrays(exp)
#(2)提取临床信息
pd <- pData(eSet)
#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#(4)提取芯片平台编号
gpl_number <- eSet@annotation
save(gse_number,pd,exp,gpl_number,file = "step1output.Rdata")

