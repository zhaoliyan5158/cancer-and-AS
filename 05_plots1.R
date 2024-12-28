#AS与control
rm(list = ls()) 
load(file = "step1output.Rdata")
load(file = "step4output1.Rdata")
#1.火山图----
library(dplyr)
library(ggplot2)
dat  = deg[!duplicated(deg$symbol),]
#####标记上下调基因
# 筛选上调和下调基因
up_genes <- dat %>%
  filter(change == "up") %>%
  arrange(desc(P.Value)) %>%
  top_n(-10, wt = P.Value)

down_genes <- dat %>%
  filter(change == "down") %>%
  arrange(desc(P.Value)) %>%
  top_n(-10, wt = P.Value)

# 合并上调和下调基因
labelled_genes <- bind_rows(up_genes, down_genes)
# 绘制火山图
p <- ggplot(data = dat, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.4, size = 3.5, aes(color = change)) +
  geom_text(data = labelled_genes, aes(label = symbol), 
            vjust = 1.5, hjust = 0.5, size = 4, check_overlap = TRUE) +  # 调整标注基因的文字大小
  ylab("-log10(Pvalue)") +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-logFC_t, logFC_t), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(P.Value_t), lty = 4, col = "black", lwd = 0.8) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 14),    # 图例文字大小
        legend.title = element_text(size = 14),  # 图例标题大小
        axis.text = element_text(size = 12),     # 坐标轴刻度文字大小
        axis.title = element_text(size = 14),    # 坐标轴标题大小
        plot.title = element_text(size = 16)    # 图标题大小)   

)


# 显示图形
print(p)


#给指定的将基因名标记在图上
for_label <- dat%>% 
  filter(symbol %in% c("HADHA","LRRFIP1"))

volcano_plot <- p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )
volcano_plot
ggsave(plot = p,filename = paste0(gse_number,"_volcano_UC.png"))

#2.差异基因热图----

load(file = 'step2output.Rdata')
# 表达矩阵行名替换
exp = exp[dat$probe_id,]
rownames(exp) = dat$symbol
#画全部差异基因的热图 32-GEO-12分钟,如果不画全部基因T改为F
if(F){
  #全部差异基因
  cg = dat$symbol[dat$change !="stable"]
  length(cg)
}else{
  #取前30上调和前30下调
  x=dat$logFC[dat$change !="stable"] 
  names(x)=dat$symbol[dat$change !="stable"] 
  cg=names(c(head(sort(x),20),tail(sort(x),20)))
  length(cg)
}
n=exp[cg,]
dim(n)

#差异基因热图
library(pheatmap)

# 创建列注释数据框
annotation_col <- data.frame(Type = Group)  # Group 是一个向量，包含样本的分组信息
rownames(annotation_col) <- colnames(n)  # 确保注释的行名与数据矩阵的列名匹配
# 定义热图颜色梯度
colors <- colorRampPalette(c("blue", "white", "red"))(100)  # 创建从蓝色到红色的100种颜色

# 为分组标签定义颜色
annotation_colors <- list(Type = c("control" = "blue", "AS" ="red" ))

# 生成热图
library(pheatmap)
library(grid)

# 首先生成热图
heat <- pheatmap(n,
                 show_colnames = FALSE,
                 show_rownames = TRUE,
                 annotation_col = annotation_col,
                 annotation_colors = annotation_colors,
                 color = colors,
                 scale = "row",
                 clustering_distance_rows = "correlation",
                 clustering_distance_cols = "correlation",
                 clustering_method = "complete",
                 fontsize_row = 10,
                 fontsize_col = 12,
                 fontsize = 14)

# 然后获取热图的 gtable 对象
gtable <- heat$gtable

# 调整图例的位置
# 注意：这里的 "legend" 和 "annotation_legend" 需要根据实际的 gtable 结构进行调整
gtable$grobs[[which(gtable$layout$name == "legend")]]$vp <- viewport(x = 0.4, y = 0.48)
gtable$grobs[[which(gtable$layout$name == "annotation_legend")]]$vp <- viewport(x = -0.2, y = 0.1)

# 使用 grid.draw 绘制调整后的热图
grid.draw(gtable)

dev.off()
# 3.感兴趣基因的箱线图----
g = c(head(cg,3),tail(cg,3))
library(tidyr)
library(tibble)
library(dplyr)
dat = t(exp[g,]) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(group = Group)

pdat = dat%>% 
  pivot_longer(cols = 2:(ncol(dat)-1),
               names_to = "gene",
               values_to = "count")

pdat$gene = factor(pdat$gene,levels = cg,ordered = T)
pdat$change = ifelse(pdat$gene %in% head(cg,10),"down","up")
library(ggplot2)
library(paletteer)
box_plot = ggplot(pdat,aes(gene,count))+
  geom_boxplot(aes(fill = group))+
  #scale_fill_manual(values = c("blue","red"))+
  scale_fill_paletteer_d("basetheme::minimal")+
  geom_jitter()+
  theme_bw()+
  facet_wrap(~change,scales = "free")
box_plot
ggsave(box_plot,filename = paste0(gse_number,"_boxplot.png"))

# 4.感兴趣基因的相关性----
library(corrplot)
M = cor(t(exp[g,]))
pheatmap(M)

my_color = rev(paletteer_d("RColorBrewer::RdYlBu"))
my_color = colorRampPalette(my_color)(10)
corrplot(M, type="upper",
         method="pie",
         order="hclust", 
         col=my_color,
         tl.col="black", 
         tl.srt=45)
library(cowplot)
cor_plot <- recordPlot() 

# 拼图
load("pca_plot.Rdata")
library(patchwork)
library(ggplotify)
(pca_plot + volcano_plot +as.ggplot(heatmap_plot))/box_plot

plot_grid(cor_plot,heatmap_plot$gtable)

