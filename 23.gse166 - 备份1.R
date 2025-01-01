library(dplyr)
library(Seurat)
library(patchwork)
#setwd("D:\\wendy\\GSE166676\\GSE166676")

# 定义数据的根目录
root_path <- "D:/wendy/GSE166676/GSE166676/"

# 定义每个样品数据的文件夹
sample_folders <- c("GSM5077727", "GSM5077728", "GSM5077729", "GSM5077730", "GSM5077731", "GSM5077732")

# 初始化一个列表来存储每个样品的Seurat对象
seurat_objects <- list()

# 循环读取每个样品的数据
for (i in 1:length(sample_folders)) {
  # 构造完整的数据路径
  data_path <- paste0(root_path, sample_folders[i], "/")
  
  # 使用Read10X函数读取数据
  data <- Read10X(data_path)
  
  # 创建Seurat对象
  seurat_objects[[i]] <- CreateSeuratObject(counts = data, project = sample_folders[i])
}

# 合并所有样品的数据为一个Seurat对象
combined_seurat_object <- Reduce(function(x, y) merge(x, y), seurat_objects)

# 输出一些基本信息以确认数据已被正确加载和合并
print(combined_seurat_object)
#查看表达矩阵
exp = combined_seurat_object@assays[["RNA"]]@layers[["counts.GSM5077727.SeuratProject.SeuratProject.SeuratProject.SeuratProject"]]
dim(exp)

#质控
combined_seurat_object[["percent.mt"]] <- PercentageFeatureSet(combined_seurat_object, pattern = "^MT-")
#鼠源的用mt
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
#过滤
combined_seurat_object <- subset(combined_seurat_object, 
               subset = nFeature_RNA > 200 & 
                 nFeature_RNA < 5000 & 
                 percent.mt < 5)
dim(combined_seurat_object)
## 4.找高变基因(HVG)
#NormalizeData每一个细胞标准化
combined_seurat_object<- NormalizeData(combined_seurat_object)
combined_seurat_object <- FindVariableFeatures(combined_seurat_object)
top10 <- head(VariableFeatures(combined_seurat_object), 10);top10
#这里选了2000个，把前十个在图上标记出来。
plot1 <- VariableFeaturePlot(combined_seurat_object)
plot2 <- LabelPoints(plot = plot1, 
                     points = top10, 
                     repel = TRUE)
plot1 + plot2
dev.off()
### 5. 标准化和降维

all.genes <- rownames(combined_seurat_object)
#ScaleData每一个基因标准化，防止后续细胞分群中，某基因表达离群，造成分群差异
combined_seurat_object <- ScaleData(combined_seurat_object, features = all.genes)
#有问题 combined_seurat_object[["RNA"]]@scale.data[30:34,1:3]
combined_seurat_object@assays[["RNA"]]@layers[["scale.data"]][30:34,1:3]
# 5.1 线性降维PCA 
combined_seurat_object<- RunPCA(combined_seurat_object, features = VariableFeatures(combined_seurat_object))

# 查看前5个主成分由哪些feature组成
print(combined_seurat_object[["pca"]], dims = 1:5, nfeatures = 5)
#前两个 点图形式
VizDimLoadings(combined_seurat_object, dims = 1:2, reduction = "pca")
#每个主成分对应基因的热图  选取前15个
DimHeatmap(combined_seurat_object, dims = 1:15, cells = 500)
# 应该选多少个主成分进行后续分析
ElbowPlot(combined_seurat_object)
# 限速步骤
f = "jc.Rdata"
if(!file.exists(f)){
  combined_seurat_object <- JackStraw(combined_seurat_object, num.replicate = 100)
  combined_seurat_object<- ScoreJackStraw(combined_seurat_object, dims = 1:15)
  save(combined_seurat_object,file = f)
}
load(f)
JackStrawPlot(combined_seurat_object, dims = 1:15)

#PC1和2
DimPlot(combined_seurat_object, reduction = "pca")+ NoLegend()
# 结合JackStrawPlot和ElbowPlot，挑选20个PC，所以这里dims定义为1:20
combined_seurat_object <- FindNeighbors(combined_seurat_object, dims = 1:15)
combined_seurat_object <- FindClusters(combined_seurat_object, resolution = 0.5) #分辨率resolutio越大获得细胞类群越多
# 结果聚成几类，用Idents查看：获得19类
length(levels(Idents(combined_seurat_object)))
#### 5.2  UMAP 和t-sne两个二选1
#PCA是线性降维，这两个是非线性降维。
#UMAP
combined_seurat_object <- RunUMAP(combined_seurat_object, dims = 1:15)
DimPlot(combined_seurat_object, reduction = "umap")
#tsen
library(Rtsne)
library(Seurat)
combined_seurat_object <- RunTSNE(combined_seurat_object, dims = 1:15)
DimPlot(combined_seurat_object, reduction = "tsne")



### 6.找marker基因

#多个样本数才用到：合并分层
head(combined_seurat_object@assays)
mergelayers <- JoinLayers(combined_seurat_object)
#鉴定cluster3和cluster 0，1之间的差异基因。如果不指定 ident.2 则鉴定cluster3 与其余clusters的差异基因。

markers_cluster5 <- FindMarkers(mergelayers, ident.1 = 3, ident.2 = c(0, 1), min.pct = 0.25)
head(markers_cluster5, n = 5)

#与所有剩余细胞相比，找到每个簇的标记
mergelayers_markers <- FindAllMarkers(mergelayers, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#min.pct = 0.25：只计算至少在(两簇细胞总数的)25%的细胞中有表达的基因
top10<- mergelayers_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)#前2个Maker
#### 6.2 marker基因的热图
DoHeatmap(mergelayers, features = top10$gene) + NoLegend()
#### 6.1 比较某个基因在几个cluster之间的表达量
VlnPlot(mergelayers, features = c("ACTA2", "MYH11","PDGFRb"))

#可以拿count数据画
#VlnPlot(combined_seurat_object, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
#在umap图上标记
#FeaturePlot(combined_seurat_object, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

## 7. 根据marker基因确定细胞：细胞注释
#new.cluster.ids <- c("Naive CD4 T", 
                     "CD14+ Mono", 
                     "Memory CD4 T",
                     "B", 
                     "CD8 T", 
                     "Vascular smooth muscle cell", #6
                     "NK", 
                     "DC", 
                     "Platelet")

#names(new.cluster.ids) <- levels(combined_seurat_object)
#combined_seurat_object<- RenameIdents(combined_seurat_object, new.cluster.ids)
#DimPlot(combined_seurat_object, 
        reduction = "umap", 
        label = TRUE, 
        pt.size = 0.5) + NoLegend()

#保存ids
#saveRDS(combined_seurat_object,'combined_seurat_object.rds')
#读取
#combined_seurat_object<-readRDS('combined_seurat_object.rds')
#singRz注释
BiocManager::install("SingleR")

library(SingleR)
library(celldex)
##载入人类参考数据集
ref_hs <- HumanPrimaryCellAtlasData()
save(ref_hs,file = "ref.hs.Rdata")
load("ref.hs.Rdata")
##载入小鼠参考数据集
#ref_mm <- MouseRNAseqData()
###把rna的转录表达数据提取
testdata <- GetAssayData(mergelayers , slot="data")
clusters <- mergelayers $seurat_clusters
cellpred <- SingleR(test = testdata, ref = ref_hs, 
                    labels = ref_hs$label.main, 
                    clusters = clusters,assay.type.test = "logcounts",
                    assay.type.ref = "logcounts")
##添加到metadata当中
celltype = data.frame(ClusterID=rownames(cellpred), 
                      celltype=cellpred$labels, stringsAsFactors = FALSE)
mergelayers@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  mergelayers@meta.data[which(mergelayers@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

DimPlot(mergelayers, reduction = "umap",label = T)
DimPlot(mergelayers, reduction = "umap", group.by = "celltype",label = F)
#更改Active Idents
#重新划分细胞亚群
new.cluster.ids <- c("Epithelial_cells",
                     "T_cells",
                     "Monocyte", 
                     "Macrophage", 
                     "Monocyte", 
                     "B_cells",
                     "Vascular smooth muscle_cells",
                     "T_cells",
                     "Neurons",
                     "Epithelial_cells",
                     "Monocyte",
                     "Keratinocytes",
                     "Endothelial_cells",
                     "B_cell") 

names(new.cluster.ids) <- levels(mergelayers)
mergelayers <- RenameIdents(mergelayers, new.cluster.ids)
DimPlot(mergelayers, reduction = "umap",label = T)
save(mergelayers,file="mergelayers.relabel.Rdata")
load("mergelayers.relabel.Rdata")

#保存ids
saveRDS(mergelayers,'mergelayers.rds')
#读取
mergelayers<-readRDS('mergelayers.rds')



