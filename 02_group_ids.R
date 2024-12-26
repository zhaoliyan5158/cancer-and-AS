# Group(实验分组)和ids(探针注释)
rm(list = ls())  
load(file = "step1output.Rdata")
library(stringr)
# 1.Group----
#如果里面含有
Group=ifelse(str_detect(pd$source_name_ch1,"control"),"control", "UC" )            
              
#设置参考水平，指定levels，对照组在前，处理组在后
Group = factor(Group,
               levels = c("control","UC"))
Group

#2.探针注释的获取-----------------
#捷径   选一种方法
library(tinyarray)
find_anno(gpl_number)
ids <- AnnoProbe::idmap('GPL570')
#方法1 BioconductorR包(最常用)
gpl_number 
#http://www.bio-info-trainee.com/1399.html
if(!require(hgu133plus2.db))BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
ls("package:hgu133plus2.db")
ids <- toTable(hgu133plus2SYMBOL)
head(ids)
# 方法2 读取GPL平台的soft文件，按列取子集
##https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570
if(T){
  #注：soft文件列名不统一，活学活用，有的表格里没有symbol列，也有的GPL平台没有提供注释
  a = getGEO(gpl_number,destdir = ".")
  b = a@dataTable@table
  colnames(b)
  ids2 = b[,c("ID","Gene Symbol")]
  colnames(ids2) = c("probe_id","symbol")
  ids2 = ids2[ids2$symbol!="" & !str_detect(ids2$symbol,"///"),]
}


# 方法3 官网下载注释文件并读取
##http://www.affymetrix.com/support/technical/byproduct.affx?product=hg-u133-plus
# 方法4 自主注释 
#https://mp.weixin.qq.com/s/mrtjpN8yDKUdCSvSUuUwcA
save(exp,Group,ids,gse_number,file = "step2output.Rdata")
