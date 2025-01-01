#差异分析
rm(list = ls(all.names = TRUE))
library(dplyr)
library(Seurat)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)

load(file = "step3-2output.Rdata")
pbmc=mergelayers
#差异分析：单细胞数据怎么分析基因Daxx在两组之间的差异表达情况，在全部细胞簇中都进行分析
Idents(pbmc)#查看默认分组
DimPlot(pbmc)
names(pbmc@meta.data)
unique(pbmc$group)
DimPlot(pbmc,split.by = 'group',label = T)
#差异分析（细胞间）
Idents(pbmc) <- "celltype"
cell_deg <- FindAllMarkers(object = pbmc, 
                               only.pos = FALSE,
                               test.use = "wilcox",
                               slot = "data", 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)

library(clusterProfiler)  # 用于功能富集分析
library(enrichplot)       # 用于富集结果的可视化
library(org.Mm.eg.db)     # Homo sapiens 的注释包，包含基因与ENTREZ ID的映射

#colnames(cell_deg)[ncol(cell_deg)] <- "SYMBOL"

#映射基因名到ENTREZ ID：：
s2e1 <- bitr(cell_deg$gene, 
             fromType = "SYMBOL",
             toType = "ENTREZID",
             OrgDb = org.Mm.eg.db)#小鼠
#其他物种http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
dim(cell_deg)
deg<- inner_join(cell_deg,s2e1,by=c("gene"="SYMBOL"))
dim(deg)
length(unique(deg$gene))
save(deg,file = "step4output1.Rdata")
load(file = "step4output1.Rdata")
#筛选出需要的内容
# 使用 dplyr 的 filter 函数进行筛选
filtered_genes_MAC <- dplyr::filter(deg, cluster == "MAC", p_val_adj < 0.05, avg_log2FC > 1)

filtered_genes_SEM <- dplyr::filter(deg, cluster == "SEM", p_val_adj < 0.05, avg_log2FC > 1)
filtered_genes_SMC <- dplyr::filter(deg, cluster == "SMC", p_val_adj < 0.05, avg_log2FC > 1)
filtered_genes_F <- dplyr::filter(deg, cluster == "Fibrochondrocyte", p_val_adj < 0.05, avg_log2FC > 1)
#KEGG富集分析-MAC
kek_MAC<- enrichKEGG(gene= filtered_genes_MAC$ENTREZID,
                  organism   = 'mmu',
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2)

# KEGG富集结果可视化
dotplot(kek_MAC,showCategory = 15) + ggtitle("KEGG Enrichment Analysis")
#将KEGG富集分析结果转换为数据框
kek_df <- as.data.frame(kek_MAC)

kek_df_Mac<- setReadable(kek_MAC, OrgDb=org.Mm.eg.db,keyType = 'ENTREZID')

#将KEGG富集分析结果转换为数据框
kek_df_Mac <- as.data.frame(kek_df_Mac)
# 将数据框写入CSV文件
write.csv(kek_df_Mac, file = "kegg_Mac.csv", row.names = FALSE)

#KEGG富集分析-SEM
kek_SEM<- enrichKEGG(gene= filtered_genes_SEM$ENTREZID,
                     organism   = 'mmu',
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.2)

# KEGG富集结果可视化
dotplot(kek_SEM,showCategory = 15) + ggtitle("KEGG Enrichment Analysis")

kek_df_SEM<- setReadable(kek_SEM, OrgDb=org.Mm.eg.db,keyType = 'ENTREZID')

#将KEGG富集分析结果转换为数据框
kek_df_SEM<- as.data.frame(kek_df_SEM)
# 将数据框写入CSV文件
write.csv(kek_df_SEM, file = "kegg_SEM.csv", row.names = FALSE)


#KEGG富集分析-smc
kek_SMC<- enrichKEGG(gene= filtered_genes_SMC$ENTREZID,
                     organism   = 'mmu',
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.2)

# KEGG富集结果可视化
dotplot(kek_SMC,showCategory = 15) + ggtitle("KEGG Enrichment Analysis")

kek_df_SMC<- setReadable(kek_SMC, OrgDb=org.Mm.eg.db,keyType = 'ENTREZID')

#将KEGG富集分析结果转换为数据框
kek_df_SMC <- as.data.frame(kek_df_SMC)
# 将数据框写入CSV文件
write.csv(kek_df_SMC, file = "kegg_SMC.csv", row.names = FALSE)
##KEGG富集分析-F
kek_F<- enrichKEGG(gene= filtered_genes_F$ENTREZID,
                     organism   = 'mmu',
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.2)

# KEGG富集结果可视化
dotplot(kek_F,showCategory = 15) + ggtitle("KEGG Enrichment Analysis")

kek_df_F<- setReadable(kek_F, OrgDb=org.Mm.eg.db,keyType = 'ENTREZID')

#将KEGG富集分析结果转换为数据框
kek_df_F<- as.data.frame(kek_df_F)
# 将数据框写入CSV文件
write.csv(kek_df_F, file = "kegg_F.csv", row.names = FALSE)

#KEGG 柱形图
barplot(kek_SMC,
        x = "Count",
        color = "p.adjust",
        showCategory = 10)#展示top10通路
####3保存数据
save(filtered_genes_MAC,filtered_genes_SMC,filtered_genes_SEM,filtered_genes_F,file = "step4-1output.Rdata")
#三种细胞富集结果合并

kek_df_Mac$CellType <- "MAC"
kek_df_SEM$CellType <- "SEM"
kek_df_SMC$CellType <- "SMC"
kek_df_F$CellType <- "Fibrochondrocyte"
#合并
combined_df <- bind_rows(kek_df_Mac, kek_df_SEM, kek_df_SMC,kek_df_F)
# 清理 Description 列，删除 '- Mus musculus (house mouse)'
combined_df<- combined_df %>%
  mutate(Description = gsub(" - Mus musculus \\(house mouse\\)$", "", Description))


#画图
#筛选前五个
top5_df <- combined_df %>%
  group_by(CellType) %>%
  top_n(n = 5, wt = Count) %>%
  ungroup()  # 取消分组，以避免后续操作受到分组的影响
library(ggplot2)


# 绘制图表
ggplot(top5_df, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  facet_grid(CellType ~ ., space = 'free_y', scales = 'free_y') +
  scale_fill_gradient(low = "#C94741", high = "#F7B799", name = "p.adjust", 
                      guide = guide_colorbar(title = "p.adjust", title.position = "top", title.hjust = 0.5)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    axis.text = element_text(size = 14, face = "bold", color = "black"),
    axis.title.x = element_text(size = 16, face = "bold", color = "black"),
    axis.title.y = element_blank(),  # 移除纵轴标题
    legend.text = element_text(size = 14, face = "bold", color = "black"),
    legend.title = element_text(size = 16, face = "bold", color = "black"),
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    strip.text = element_text(size = 14, face = "bold", color = "black")
  )

save(top5_df,kek_df_Mac,kek_df_SEM,kek_df_SMC,file = "step5-1output.Rdata")








pbmc$celltype.group <- paste(pbmc$celltype, pbmc$group, sep = "_")
pbmc$celltype <- Idents(pbmc)
Idents(pbmc) <- "celltype.group"

mydeg_Macrophages <- FindMarkers(pbmc,ident.1 = 'Macrophages_AAA',ident.2 = 'Macrophages_Sham', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
#解放生产力 通过循环自动计算差异基因(组间)
cellfordeg<-levels(pbmc$celltype)
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(pbmc, ident.1 = paste0(cellfordeg[i],"_Sham"), ident.2 = paste0(cellfordeg[i],"_AAA"), verbose = FALSE)
  write.csv(CELLDEG,paste0(cellfordeg[i],".CSV"))
}
list.files()
# 用VlnPlot可视化Daxx基因在不同细胞类型和组别之间的表达
VlnPlot(pbmc, features = "Daxx", group.by = "celltype.group")

FeaturePlot(pbmc, features = "Daxx")
DotPlot(pbmc, features = "Daxx") + RotatedAxis()

library(Seurat)
library(ggplot2)

# 假设你的Seurat对象名为pbmc
# 提取Daxx基因表达数据和细胞群体信息
data_for_plot <- FetchData(pbmc, vars = c("Daxx", "celltype.group"))
ggplot(data_for_plot, aes(x = celltype.group, y = Daxx)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Expression of Daxx across Cell Groups", x = "Cell Group", y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

