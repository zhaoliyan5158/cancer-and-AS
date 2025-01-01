library(dplyr)
library(Seurat)
library(patchwork)
library(Seurat)
library(hdf5r)
library(ggplot2)
library(Matrix)
library(Seurat)
#单细胞txt格式文件读取
dir<- "D:/wendy/GSE155514-scRNA-AS/GSE155514-vsmc/"

samples=list.files( dir ,pattern = 'gz')
samples 
library(data.table)
ctList = lapply(samples,function(pro){ 
  print(pro)
  ct=fread(file.path( dir ,pro),data.table = F)
  ct[1:4,1:4]
  rownames(ct)=ct[,1]
  colnames(ct) = paste(gsub('.txt.gz','',pro),
                       colnames(ct) ,sep = '_')
  ct=ct[,-1] 
  return(ct)
})
lapply(ctList, dim)
tmp =table(unlist(lapply(ctList, rownames)))
cg = names(tmp)[tmp==length(samples)]
bigct = do.call(cbind,
                lapply(ctList,function(ct){ 
                  ct = ct[cg,] 
                  return(ct)
                }))
combined_seurat_object=CreateSeuratObject(counts =  bigct, 
                           min.cells = 5,
                           min.features = 300)

# 输出一些基本信息以确认数据已被正确加载和合并
print(combined_seurat_object)
#查看表达矩阵
exp = combined_seurat_object@assays[["RNA"]]@layers[["counts"]]
dim(exp)

#[1] 12307 15916
#质控：计算线粒体RNA，高线粒体基因表达可能表明细胞受损。#这里是小鼠
combined_seurat_object[["percent.mt"]] <- PercentageFeatureSet(combined_seurat_object, pattern = "^mt-")
#人的MT,鼠源的用mt,如果线粒体基因为0，注意人和鼠的区别
head(combined_seurat_object@meta.data)

VlnPlot(combined_seurat_object, 
        features = c("nFeature_RNA",
                     "nCount_RNA", 
                     "percent.mt"), 
        ncol = 3)
#ncol = 3没行展示3个图
#相关性
plot1 <- FeatureScatter(combined_seurat_object, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(combined_seurat_object, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")

plot1 + plot2
#过滤：基因数量过少是死细胞，过多可能是双细胞，需要过滤
combined_seurat_object <- subset(combined_seurat_object, 
               subset = nFeature_RNA > 200 & 
                 nFeature_RNA < 4000 & 
                 percent.mt < 5)

dim(combined_seurat_object)
## 4.找高变基因(HVG)
#NormalizeData每一个细胞标准化：标准化数据存在data中
combined_seurat_object<- NormalizeData(combined_seurat_object,normalization.method = "LogNormalize", scale.factor = 10000)
#使用方差稳定转换（VST）方法选取2000个变异性最高的特征（基因）。这些特征（基因）通常用于后续的主成分分析（PCA）和数据聚类。
combined_seurat_object <- FindVariableFeatures(combined_seurat_object,selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(combined_seurat_object), 10)
#这里选了2000个，把前十个在图上标记出来。
plot1 <- VariableFeaturePlot(combined_seurat_object)
plot2 <- LabelPoints(plot = plot1, 
                     points = top10, 
                     repel = TRUE)
plot1 + plot2
dev.off()
#保存数据
save(combined_seurat_object,file = "step1output.Rdata")

load(file = "step1output.Rdata")
### 5. 标准化和降维

all.genes <- rownames(combined_seurat_object)
#ScaleData每一个基因标准化，防止后续细胞分群中，某基因表达离群，造成分群差异

combined_seurat_object <- ScaleData(combined_seurat_object, features = all.genes)#运行后显示内存不足
#只使高变基因标准化
#combined_seurat_object <- ScaleData(combined_seurat_object)


# 检查处理结果
print(combined_seurat_object)

#有问题 combined_seurat_object[["RNA"]]@scale.data[30:34,1:3]
combined_seurat_object@assays[["RNA"]]@layers[["scale.data"]][30:34,1:3]
# 5.1 线性降维PCA 
combined_seurat_object<- RunPCA(combined_seurat_object, features = VariableFeatures(combined_seurat_object))


# 查看前5个主成分由哪些feature组成
print(combined_seurat_object[["pca"]], dims = 1:5, nfeatures = 5)
#前两个点图形式
VizDimLoadings(combined_seurat_object, dims = 1:2, reduction = "pca")
#每个主成分对应基因的热图  选取前15个
DimHeatmap(combined_seurat_object, dims = 1:15, cells = 500)
combined_seurat_object1<-combined_seurat_object#备份

# 应该选多少个主成分进行后续分析
ElbowPlot(combined_seurat_object)

# 限速步骤
combined_seurat_object <- JackStraw(combined_seurat_object, num.replicate = 100)
combined_seurat_object<- ScoreJackStraw(combined_seurat_object, dims = 1:15)
JackStrawPlot(combined_seurat_object, dims = 1:15)


#PC1和2
DimPlot(combined_seurat_object, reduction = "pca")+ NoLegend()
# 结合JackStrawPlot和ElbowPlot，挑选5个PC，所以这里dims定义为1:5
combined_seurat_object <- FindNeighbors(combined_seurat_object, dims = 1:5)
combined_seurat_object <- FindClusters(combined_seurat_object, resolution = 0.2) #分辨率resolutio越大获得细胞类群越多
# 结果聚成几类，用Idents查看：获得5类
length(levels(Idents(combined_seurat_object)))
#### 5.2  UMAP 和t-sne两个二选1
#PCA是线性降维，这两个是非线性降维。
#UMAP
combined_seurat_object <- RunUMAP(combined_seurat_object, dims = 1:5)
DimPlot(combined_seurat_object, reduction = "umap")



### 6.找marker基因

#多个样本数才用到：合并分层
head(combined_seurat_object@assays)
mergelayers <- JoinLayers(combined_seurat_object)
#加入分组信息
#查看默认分组
table(mergelayers@meta.data$orig.ident)
#读取临床信息
patients_metadata<-data.table::fread("clinic.csv",header = T)
head(patients_metadata)
#从seurt对象提取orig.ident
metadata <- FetchData(mergelayers,"orig.ident")
metadata$cell_id<-rownames(metadata)
#dplyr中left_join合并
library(dplyr)
metadata=left_join(x=metadata,y=patients_metadata,by="orig.ident")
#重新加行名
rownames(metadata)<-metadata$cell_id
#通过AddMetadata直接加入4列
mergelayers<-AddMetaData(mergelayers,metadata =metadata )
table(mergelayers@meta.data$group)
#保存
saveRDS(mergelayers,'mergelayers.rds')
#读取
mergelayers<-readRDS('mergelayers.rds')

####找marker基因（按照细胞簇）
mergelayers_markers <- FindAllMarkers(mergelayers, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#min.pct = 0.25：只计算至少在(两簇细胞总数的)25%的细胞中有表达的基因
####3保存数据
save(mergelayers,mergelayers_markers ,file = "step2output.Rdata")

load(file = "step2output.Rdata")

top10<- mergelayers_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)#前2个Maker
#### 6.2 marker基因的热图：查看maker基因是否是细胞簇的标记基因
DoHeatmap(mergelayers, features = top10$gene) + NoLegend()

VlnPlot(mergelayers, features = top10$gene[1:3],pt.size=0)

DimPlot(mergelayers,label = T)#cluster

#查看有哪些基因
gene_names <- Features(mergelayers)
gene_names_df <- data.frame(GeneNames = gene_names)

library(openxlsx)
write.xlsx(gene_names_df, "Gene_Names.xlsx")
#singR注释
#BiocManager::install("SingleR")

library(SingleR)
library(celldex)
##载入人类参考数据集
#ref_ms <- MouseRNAseqData()
#save(ref_hs,file = "ref.ms.Rdata")
#load("ref.hs.Rdata")
##载入小鼠参考数据集
ref_mm <- MouseRNAseqData()
save(ref_mm,file = "ref_mm.Rdata")
load("ref_mm.Rdata")
###把rna的转录表达数据提取
testdata <- GetAssayData(mergelayers, layer = "data")

clusters <- mergelayers$seurat_clusters
cellpred <- SingleR(test = testdata, ref = ref_mm, 
                    labels = ref_mm$label.main, 
                    clusters = clusters,assay.type.test = "logcounts",
                    assay.type.ref = "logcounts")

##添加到metadata当中
celltype = data.frame(ClusterID=rownames(cellpred), 
                      celltype=cellpred$labels, stringsAsFactors = FALSE)
mergelayers@meta.data$celltype = "NA"
#AS全部细胞分类

#toms_markers5<- c('Acta2', 'Tagln', 'Myh11', 
                   'Myl9', 'Cald1', 'Pln', 
                   'Rgs5', 'Dcn', 'Col3a1',
                   'Apod', 'Cfd', 'Fbln1', 
                   'Lum', 'Bgn', 'Fn1', 'Eln', 
                   'Timp1', 'Vcan', 'Mgp', 'Mmp1',
                   'Mmp2', 'Mmp3', 'Mmp9', 'Mmp19', 
                   'Ccl2', 'Cxcl1', 'Cxcl2', 'Cxcl3', 
                   'Gdf10', 'Sox9', 'Lgals3', 'Klf4', 
                   'Vcam1', 'Runx2', 'Spp1', 'Trpv4', 'S100b')
#toms_markers<- c('Csrp2', 'Pln', 'Ramp1', 
                  'C12orf75', 'Mfap4', 'Hp1bp3', 
                  'Actn1', 'Srsf5', 'Hnrnpa2b1',
                  'Map1b', 'Dlx6as1', 'Sost', 
                  'Dlx5', 'Lmo2', 'Atp1b1', 'Tnfrsf11b', 
                  'Omd', 'Prss23', 'Aspn', 'Col8a1',
                  'Aif1', 'Fcer1g', 'Tyrobp', 'Laptm5', 
                  'Ctss', 'Steap4', 'Ccdc102b', 'Apoe', 
                  'Ndufa4l2', 'Igfbp4')
toms_markers<- c('Acta2', 'Tagln', 'Myh11',
                 'Myl9', 'Cald1', 
                 'Ly6a', 'Vcam1', 'Ly6c1',
               'Fn1', 'Col1a1', 'Col1a2', 'Acan','Sox9',
              'Lgals3','Cd68')
# 创建点图
library(ggplot2)
#install.packages("remotes")
#remotes::install_github("satijalab/seurat", ref = "develop")
#Idents(mergelayers) <- "seurat_clusters"
dotplot <- DotPlot(mergelayers, features = toms_markers) + 
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# 绘制点图
print(dotplot)






library(dplyr)

#修改注释结果
celltype <- celltype %>%
  mutate(celltype = case_when(
    ClusterID == "0" ~ "SMC",
    ClusterID == "1" ~ "Fibrochondrocyte",
    ClusterID == "2" ~ "MAC",
    ClusterID == "3" ~ "SEM",
    ClusterID == "4" ~ "SMC",
   TRUE ~ celltype  # 对于不匹配上述条件的细胞保持原标签不变
  ))

for(i in 1:nrow(celltype)){
  mergelayers@meta.data[which(mergelayers@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
library(patchwork)
library(ggplot2)
library(Seurat) # 确保加载了Seurat库

####3保存数据
save(mergelayers,mergelayers_markers ,file = "step3-1output.Rdata")

load(file = "step3-1output.Rdata")

# 生成两个DimPlot图形对象，并设置字体大小和加粗
p1 <- DimPlot(mergelayers, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  theme(text = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 13, face = "bold")) # 设置字体大小为13，加粗

p2 <- DimPlot(mergelayers, reduction = "umap", group.by = "celltype", label = TRUE) +
  theme(text = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 13, face = "bold")) # 设置字体大小为13，加粗

# 使用patchwork操作符来布局图形
p1 + p2


# 确保 'group' 是因子并设置水平顺序

mergelayers$group <- factor(mergelayers$group, levels = c("0 week", "8 weeks", "16 weeks", "26 weeks"))
DimPlot(mergelayers, reduction = "umap", group.by = "celltype", split.by = "group",label = T)+
  theme(text = element_text(size = 13, face = "bold"),
axis.text = element_text(size = 13, face = "bold"))
#DimPlot(mergelayers, reduction = "umap", split.by = "group",label = T)#不显示celltype标签
#确定活动条目
Idents(mergelayers) <- mergelayers$celltype

dotplot <- DotPlot(mergelayers, features = toms_markers) + 
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),  # 可以调整Y轴文本大小
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill=NA, size=1)  # 添加框线
  )

# 绘制点图
print(dotplot)
#不同分组UMAP

# 绘制特定基因的表达图
library(Seurat)
library(cowplot)

# 定义一个函数来创建主题
my_theme <- function() {
  theme(
    text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 16, face = "bold")
  )
}

# 指定基因列表
genes <- c("Fhod1", "Irf7")

# 生成基因表达图列表
gene_plots <- lapply(genes, function(gene) {
  FeaturePlot(mergelayers, features = gene, reduction = "umap", split.by = "group", cols = c("lightgrey", "red")) 
})

# 组合所有基因表达图
combined_gene_plots <- plot_grid(plotlist = gene_plots, ncol = 1)  # 单列显示

# 打印组合图
print(combined_gene_plots)
dev.off()

# 假设您的Seurat对象中的细胞类型存储在 `celltype` 字段，条件存储在 `condition` 字段
mergelayers$ident <- paste(mergelayers$celltype, mergelayers$group, sep = "_")
# 绘制点图
library(Seurat)

dot_plot <- DotPlot(mergelayers, features = genes, group.by = "celltype") +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),  # 如果不需要X轴标题，保持这个设置
    axis.title.y = element_blank(),  # 如果不需要Y轴标题，保持这个设置
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, face = "bold", color = "black"),  # X轴标签倾斜并加粗，颜色设置为黑色
    axis.text.y = element_text(size = 14, face = "bold", color = "black"),  # Y轴标签加粗，颜色设置为黑色
    axis.line = element_line(color = "black"),  # 添加轴线并设置颜色为黑色
    axis.ticks = element_line(color = "black")  # 确保刻度线显示并设置颜色为黑色
  )

print(dot_plot)

ggsave("dot_plot.pdf", plot = dot_plot, width =4, height = 4)

####3保存数据
save(mergelayers,mergelayers_markers ,file = "step3-2output.Rdata")

load(file = "step3-2output.Rdata")


# 按细胞类型查找标记基因
mergelayers_markers_celltype <- FindAllMarkers(mergelayers, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, group.by = "celltype")

# 提取每个细胞类型的Top 5基因
top5_markers <- mergelayers_markers_celltype %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  distinct(gene, .keep_all = TRUE)  # 确保基因名是唯一的
#确定活动条目
dot_plot <- DotPlot(mergelayers, features = top5_markers$gene, group.by = "celltype") +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),  # 如果不需要X轴标题，保持这个设置
    axis.title.y = element_blank(),  # 如果不需要Y轴标题，保持这个设置
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, face = "bold", color = "black"),  # X轴标签倾斜并加粗，颜色设置为黑色
    axis.text.y = element_text(size = 14, face = "bold", color = "black"),  # Y轴标签加粗，颜色设置为黑色
    axis.line = element_line(color = "black"),  # 添加轴线并设置颜色为黑色
    axis.ticks = element_line(color = "black")  # 确保刻度线显示并设置颜色为黑色
  )
print(dot_plot)


ggsave("top5.pdf", plot = dot_plot, width = 6, height = 4)





















