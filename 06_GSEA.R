#AS与control
rm(list = ls()) 
#load(file = "step1output.Rdata")
load(file = "step4output1.Rdata")
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggridges)
#按照logFC降序排列
data_all_sort <- deg %>% 
  arrange(desc(logFC))
#
geneList = data_all_sort$logFC 
names(geneList) <- data_all_sort$ENTREZID 
KEGG_database="hsa"
#gsea结果富集：gsea→result中
gsea_UC<- gseKEGG(geneList, organism = KEGG_database, pvalueCutoff = 0.05)
#将ENTREZID 转换为基因名字
gsea_UC<- setReadable(gsea_UC, OrgDb=org.Hs.eg.db,keyType = 'ENTREZID')
#气泡图

dotplot(gsea_UC)
#山体图label_format = 50改变字的重叠
ridgeplot(gsea_UC,label_format = 50)
#GSEA标准图:gsea,1→代表画几条
gseaplot2(gsea_UC, 1:5, pvalue_table =F, base_size = 14)
#后5条通路
# 筛选下调通路 (NES < 0)
downregulated <- gsea_UC@result[gsea_UC@result$NES < 0, ]

# 按 NES 值降序排列（负值绝对值最大排在前）
downregulated <- downregulated[order(downregulated$NES), ]

# 提取 NES 负值最小（绝对值最大）的后 5 条通路
last5 <- rownames(head(downregulated, 5))

# 查看筛选结果
print(last5)

# 绘制下调的后 5 条通路
gseaplot2(gsea_UC, last5, pvalue_table = FALSE, base_size = 14)

dev.off()

save(gsea_UC,deg,file = "gsea_outputUC.Rdata")
load(file = "gsea_outputUC.Rdata")
library(openxlsx)
write.xlsx(gsea_UC@result, file = "GSEA_Results_AS.xlsx")
table(Group)
