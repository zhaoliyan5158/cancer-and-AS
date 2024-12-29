#ssGSEA分析免疫细胞浸润分析（TCGA样本） 散装生信哥
rm(list = ls())
library(tidyverse)   # 数据清洗与处理
library(pheatmap)    # 绘制热图
library(ggplot2)     # 数据可视化
library(ggpubr)      # 增强的ggplot图形
library(paletteer)   # 颜色管理
library(corrplot)     # 相关性矩阵图
library(rstatix)     # 辅助统计分析
library(GSVA)        # 基因集变异分析
library(stringr) 

load(file = "step2output.Rdata")
class(exp)
exp=as.data.frame(exp)
ids2=as.data.frame(ids2)

# 将探针名设置为行名
rownames(ids2) <- ids2$ID

# 根据探针名合并两个数据框
exp <- merge(ids2, exp, by = "row.names", all.y = TRUE)
# 删除不再需要的列
exp1 <- exp[, -c(1,2)]  # 假设探针名和基因名是前两列
#删除没有基因名的行
exp1 <- exp1 %>% 
  filter(!is.na(GENE_SYMBOL) & GENE_SYMBOL != "")
# 对重复的基因名取平均值
exp1 <- exp1  %>%
  group_by(GENE_SYMBOL) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
# 将GENE_SYMBOL列设置为行名
rownames(exp1) <- exp1$GENE_SYMBOL
# 删除不再需要的列（注意删除后行名是否在）
exp2 <- exp1[ , -1]  # 删除第一列
# 保存行名
row_names <- exp1$GENE_SYMBOL

# 删除列
exp2 <- exp1[ , -1]

# 重新设置行名
rownames(exp2) <- row_names
save(exp2,Group,file = "exp_symbol.Rdata")#全部基因表达矩阵带基因名
load(file = "exp_symbol.Rdata")

# 从网络资源加载额外的基因集数据
cellMarker1 <- read.delim(file("CellReports.txt", encoding = "UTF-16LE"), header = FALSE, sep = "\t")
cellMarker1 <- cellMarker1 %>% column_to_rownames("V1") %>% t()

# 假定 cellMarker1 是已经读取的数据框

a <- cellMarker1  # 将数据框赋值给变量 a
a <- a[1:nrow(a), ]  # 这行代码实际上没有改变a的内容，只是重复赋值

set <- colnames(a)  # 获取数据框a的所有列名，存储在变量set中

geneSet <- list()  # 初始化一个空列表用来存储基因集

for (i in set) {  # 遍历每一个列名
  x <- as.character(a[, i])  # 提取当前列的数据，并转换成字符向量
  x <- x[nchar(x) != 0]  # 移除那些为空的元素（即字符长度为0的元素）
  geneSet[[i]] <- x  # 将处理后的非空基因名列表存储在geneSet列表中，以列名为键
}
# 使用GSVA包的gsva函数进行基因集变异分析
res <- gsva(as.matrix(exp2),  # 将exp数据框转换为矩阵形式
            geneSet,         # geneSet列表，包含基因集
            method = "ssgsea",  # 分析方法选择ssGSEA
            kcdf = "Gaussian",  # 核密度估计函数选择高斯分布
            mx.diff = FALSE    # 最大差异化参数设置为FALSE
)
#分组信息
load(file = "step1output.Rdata")
#如果里面含有
group_list<-ifelse(str_detect(pd$source_name_ch1,"Control"),"control", "AS" )      
table(group_list)
annotation <- data.frame(group_list)
rownames(annotation) <- colnames(res)
head(annotation)
# 创建一个新的数据框 resm，初始与 res 相同
resm <- res

# 循环遍历 res 的每一列
for (i in colnames(res)) {
  # 对每一列应用最小-最大归一化
  resm[, i] <- (res[, i] - min(res[, i])) / (max(res[, i]) - min(res[, i]))
}
#热图
pheatmap(res,
         show_colnames = F,
         annotation_col = annotation,
         fontsize = 10)
#箱线图
dt <- resm %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  gather(key = cell_type, value = value, -sample)
head(dt)
dtt <- dt %>%
  group_by(sample) %>%
  mutate(proportion = round(value / sum(value), 3))
head(dtt)
dtt$cell_type <- factor(dtt$cell_type, levels = unique(rownames(res)))

mytheme <- theme(
  axis.title = element_text(size = 12),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
  legend.position = "bottom",
  legend.text = element_text(size = 10)
)

# 安装 paletter 包
# install.packages("paletter")
# 加载 paletter 包
library(paletteer)

# 查看 paletter 包中所有可用的调色板名称
d_palettes <- palettes_d_names

# 从 paletter 包中选择一个调色板 (khroma::smoothrainbow)
col <- paletteer_d("khroma::smoothrainbow", n = 28)


# 使用 ggplot2 绘制箱线图
p <- ggplot(dtt, aes(x = cell_type, y = proportion, fill = cell_type)) +
  geom_boxplot(color = "black", alpha = 0.6, outlier.shape = 21, outlier.size = 1.2) +
  scale_fill_manual(values = col) +
  labs(x = "cell type", y = "proportion") +
  theme_bw() + mytheme

# 打印图形
print(p)

##### 分组箱线图 + 显著性标记 #####
# 不同分组（肿瘤样本和正常样本）不同类型免疫细胞的表达情况


# 新增分组列
# 将annotation中的行名变成一列
annotation <- annotation %>% rownames_to_column(var = "sample")

# 合并dtt和annotation数据框
dtt <- left_join(dtt, annotation, by = "sample")

# 查看合并后的数据框
head(dtt)
# 重命名列名
dtt <- dtt %>% rename(group = group_list)
# 确保 group 列按所需顺序排列：设置因子水平
dtt$group <- factor(dtt$group, levels = c("control", "AS"))
# 分组箱线图展示
p2 <- ggplot(dtt,
             aes(x = cell_type, y = proportion, fill = group)) +
  geom_boxplot(color = "black", alpha = 0.6, outlier.shape = 21, outlier.size = 1.2) +
  scale_fill_manual(values = c("#4979b6", "#d9332a")) +
  labs(x = "cell type", y = "proportion") +
  theme_bw() + mytheme + theme(axis.text.x = element_text(angle = 45))

# 打印图形
p2
## 使用 t test 或 wilcox test 进行两两比较（t 检验为例）：
t <- t_test(group_by(dtt, cell_type), proportion ~ group)
tj <- adjust_pvalue(t, method = 'fdr')  # p 值矫正；
tj

# 根据 p.adj 添加显著性标记符号；
tj <- add_significance(tj, 'p.adj')
tj

# 在图表中添加 p 值或显著性标记；
lab <- add_xy_position(tj, x = 'cell_type', dodge = 0.65)
# ggpubr 绘图：
library(ggpubr)

p3 <- ggboxplot(dtt, x = "cell_type", y = "proportion", 
                fill = "group", color = "black",
                width = 0.7, alpha = 0.6,
                outlier.shape = 21, outlier.size = 1.2) +
  scale_fill_manual(values = c("#4979b6", "#d9332a")) +
  labs(x = "", y = "Expression") +
  theme_bw() + 
  mytheme + 
  theme(axis.text.x = element_text(angle = 50, face = "bold", size = 10, color = "black",hjust =1,vjust = 1),  # 加粗并调整X轴字体大小，设置颜色为黑色
        axis.text.y = element_text(face = "bold", size = 12, color = "black"),  # 加粗并调整Y轴字体大小，设置颜色为黑色
        axis.title.x = element_text(face = "bold", size = 12, color = "black"),  # 加粗并调整X轴标题字体大小，设置颜色为黑色
        axis.title.y = element_text(face = "bold", size = 12, color = "black"),  # 加粗并调整Y轴标题字体大小，设置颜色为黑色
        legend.position = "right",  # 图例放在图的右边
        legend.box = "vertical",  # 图例分组竖着放
        legend.text = element_text(face = "bold", color = "black"),  # 设置图例文本颜色为黑色并加粗
        legend.title = element_text(face = "bold", color = "black")) + # 设置图例标题颜色为黑色并加粗
  stat_pvalue_manual(lab, label = "p.adj.signif", label.size = 3, bracket.size = 0.5, tip.length = 0.01)

p3

### 免疫细胞相关性热图 ###
# 查看 resm 矩阵格式，再看一下矩阵格式：每行为一类细胞，每列一个样本；
View(resm)

# 计算相关性系数；
resmcor <- cor(t(resm))
View(resmcor)

# 绘制相关性热图
corrplot(resmcor,
         method = "square",
         order = "hclust",
         tl.cex = 0.6,
         tl.col = "black")

# 添加显著性标记；
resmorp <- cor.mtest(resmcor, conf.level = .95)  # 使用 cor.mtest 做显著性检验；

# 提取 p 值矩阵；
p.mat <- resmorp$p
View(p.mat)
#相关性热图中展示显著性标记
library(corrplot)

# 生成相关性热图
corrplot(resmcor, 
         method = "color", 
         type = "full", 
         order = "hclust", 
         tl.cex = 1.0,  # 调整标签字体大小并加粗，这里设置为1.0，您可以根据需要调整
         tl.col = "black",  # 设置标签颜色为黑色
         p.mat = resmorp$p, 
         sig.level = c(0.001, 0.01, 0.05), 
         outline = "white", 
         insig = "label_sig", 
         pch.cex = 1.0,  # 调整显著性标记的大小
         pch.col = "black",
         diag = TRUE)  # 显示对角线

# 在相同的热图上添加相关系数的数字
corrplot(resmcor, 
         method = "number",  # 显示相关系数的数字
         type = "lower",      # 仅处理左下角的矩阵部分
         add = TRUE,          # 将此图层添加到现有的corrplot图上
         number.cex = 0.6,    # 调整数字大小，这里设置为1.0，您可以根据需要调整
         order = "hclust",    # 使用层次聚类排序，使相关性模式更清晰
         tl.pos = "n",        # 不显示任何标签
         diag = FALSE)        # 不显示对角线

dev.off()
# 免疫细胞与目标基因相关性热图
identical(colnames(resm), colnames(exp2))

# 指定基因列表
genes <- c("FHOD1", "IRF7", "TNF","ZSWIM3")

# 从exp2数据框中提取指定基因的数据
# 转换矩阵为数据框
exp_genes <- exp2[genes,]
rownames(exp_genes) <- genes
# 将resm和exp_genes合并为一个新的数据框
rb <- rbind(resm, exp_genes)
rownames(rb)

# 计算合并后数据框的转置的相关性矩阵
rbcor <- cor(t(rb))

# 对相关性矩阵进行显著性检验
rbcorp <- cor.mtest(rbcor, conf.level = .95)

# 提取p值矩阵
p.mat2 <- rbcorp$p
# 提取相关性矩阵或p值矩阵的特定行和列
split <- rbcor[1:nrow(resm), (nrow(rb) - length(genes) + 1):nrow(rb)]
View(split)  # 查看切割后的部分

# 同样的操作，这次针对p值矩阵
splitp <- p.mat2[1:nrow(resm),#行取免疫细胞所在行
                 (nrow(rb) - length(genes) + 1):nrow(rb)]#列取目的基因所在行
View(splitp)

mark <- ifelse(
  splitp < 0.001, "***", 
  ifelse(splitp < 0.01, "**", 
         ifelse(splitp < 0.05, "*", ""))
)
mark <- matrix(mark, nrow = nrow(splitp))


# 创建标记矩阵，显示显著性水平(备选)
###mark <- matrix(case_when(
  #splitp < 0.001 ~ "***",
 # splitp < 0.01 ~ "**",
  ##splitp < 0.05 ~ "*",
  #T ~ ""), 
 # nrow = nrow(splitp))

# 创建颜色映射，用于热图或相关性图的颜色表示
col2 <- colorRampPalette(c("greenyellow", "white", "red"))(95)
pheatmap(t(split),
         display_numbers = t(mark),
         number_color = "black",
         fontsize_number = 13,  # 调整数字的字体大小
         color = col2,
         border_color = "white",
         treeheight_col = 0,  # 不显示列的树状图
         treeheight_row = 0,  # 不显示行的树状图
         fontsize_col = 12,  # 增大列标签的字体大小
         fontface_col = "bold",  # 将列标签字体设置为加粗
         fontsize_row = 12,  # 增大行标签的字体大小
         fontface_row = "bold",  # 将行标签字体设置为加粗
         angle_col = 90)  # 列标签旋转90度











