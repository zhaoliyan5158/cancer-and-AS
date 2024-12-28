#Hub基因在两组之间的箱线图
rm(list = ls()) 
library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)
library(paletteer)
library(tidyverse)
library(ggsignif)

load(file = "step6output.Rdata")

g = c("TNF","ZSWIM3","FHOD1","IRF7")
cg = c("TNF","ZSWIM3","FHOD1","IRF7")
dat = t(exp[g,]) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(group = Group)

pdat = dat%>% 
  pivot_longer(cols = 2:(ncol(dat)-1),
               names_to = "gene",
               values_to = "count")

pdat$gene = factor(pdat$gene,levels = cg,ordered = T)


# 只筛选出 TNF 基因的数据
tnf_data <- pdat %>% filter(gene == "TNF")

# 计算 TNF 的 P 值
tnf_plot <- ggplot(tnf_data, aes(x = group, y = count)) +
  geom_boxplot(aes(fill = group)) +
  geom_signif(comparisons = list(c("control", "AS")),
              annotations = sprintf("P = %.2e", tnf_p_value),
              y_position = max(tnf_data$count) + 0.2,
              tip_length = 0.02) +
  scale_fill_paletteer_d("basetheme::minimal") +
  theme_bw() +
  xlab("") + 
  ylab("TNF") +  
  theme(
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold"),  # 加粗X轴标签字体
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "none"  # 移除图例
  )

print(tnf_plot)
ggsave("tnf.pdf", plot = tnf_plot, width = 4, height = 4, device = "pdf")



# 只筛选出 FHOD1 基因的数据
FHOD1_data <- pdat %>% filter(gene == "FHOD1")

# 计算 FHOD1 的 P 值
FHOD1_p_value <- FHOD1_data %>%
  summarise(
    p_value = {
      control <- count[group == "control"]
      as <- count[group == "AS"]
      t_test <- t.test(control, as)
      t_test$p.value
    }
  ) %>%
  pull(p_value)

FHOD1_plot <- ggplot(FHOD1_data, aes(x = group, y = count)) +
  geom_boxplot(aes(fill = group)) +
  geom_signif(comparisons = list(c("control", "AS")),
              annotations = sprintf("P = %.2e", FHOD1_p_value),
              y_position = max(FHOD1_data$count) + 0.2, # 使用正确的数据集来调整y轴位置
              tip_length = 0.02) +
  scale_fill_paletteer_d("basetheme::minimal") +
  theme_bw() +
  xlab("") + # 移除 X 轴标签
  ylab("FHOD1") +  # 设置 Y 轴标签
  theme(
    panel.border = element_blank(),   # 移除四周的框线
    axis.line = element_line(color = "black"),  # 只保留下方和左侧的轴线
    panel.grid.major = element_blank(),  # 去除主要网格线
    panel.grid.minor = element_blank(),   # 去除次要网格线
    legend.position = "none",  # 移除图例
    axis.text.x = element_text(size = 12, face = "bold"),  # 加粗并设置X轴标签字体大小
    axis.text.y = element_text(size = 12),  # 设置Y轴标签字体大小
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

# 打印绘制的图形
print(FHOD1_plot)
# 保存图形为PDF文件
ggsave("FHOD1.pdf", plot = FHOD1_plot, width = 4, height = 4, device = "pdf")
# 只筛选出 IRF7 基因的数据
irf7_data <- pdat %>% filter(gene == "IRF7")

# 计算 IRF7 的 P 值
irf7_p_value <- irf7_data %>%
  summarise(
    p_value = {
      control <- count[group == "control"]
      as <- count[group == "AS"]
      t_test <- t.test(control, as)
      t_test$p.value
    }
  ) %>%
  pull(p_value)

irf7_plot <- ggplot(irf7_data, aes(x = group, y = count)) +
  geom_boxplot(aes(fill = group)) +
  geom_signif(comparisons = list(c("control", "AS")),
              annotations = sprintf("P = %.2e", irf7_p_value),
              y_position = max(irf7_data$count) + 0.5, # 增加y轴位置，使P值显示在图顶上方
              tip_length = 0.02) +
  scale_fill_paletteer_d("basetheme::minimal") +
  theme_bw() +
  xlab("") + # 移除 X 轴标签
  ylab("IRF7") +  # 设置 Y 轴标签
  theme(
    panel.border = element_blank(),   # 移除四周的框线
    axis.line = element_line(color = "black"),  # 只保留下方和左侧的轴线
    panel.grid.major = element_blank(),  # 去除主要网格线
    panel.grid.minor = element_blank(),   # 去除次要网格线
    legend.position = "none",  # 移除图例
    axis.text.x = element_text(size = 12, face = "bold"),  # 加粗并设置X轴标签字体大小
    axis.text.y = element_text(size = 12),  # 设置Y轴标签字体大小
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

# 打印绘制的图形
print(irf7_plot)

# 保存图形为PDF文件
ggsave("irf7.pdf", plot = irf7_plot, width = 4, height = 4, device = "pdf")

# 只筛选出 ZSWIM3 基因的数据
ZSWIM3_data <- pdat %>% filter(gene == "ZSWIM3")

# 计算 ZSWIM3 的 P 值
ZSWIM3_p_value <- ZSWIM3_data %>%
  summarise(
    p_value = {
      control <- count[group == "control"]
      as <- count[group == "AS"]
      t_test <- t.test(control, as)
      t_test$p.value
    }
  ) %>%
  pull(p_value)

ZSWIM3_plot <- ggplot(ZSWIM3_data, aes(x = group, y = count)) +
  geom_boxplot(aes(fill = group)) +
  geom_signif(comparisons = list(c("control", "AS")),
              annotations = sprintf("P = %.2e", ZSWIM3_p_value), # 使用正确的P值变量
              y_position = max(ZSWIM3_data$count) + 0.2, # 增加y轴位置，使P值显示在图顶上
              tip_length = 0.02) +
  scale_fill_paletteer_d("basetheme::minimal") +
  theme_bw() +
  xlab("") + # 移除 X 轴标签
  ylab("ZSWIM3") +  # 设置 Y 轴标签
  theme(
    panel.border = element_blank(),   # 移除四周的框线
    axis.line = element_line(color = "black"),  # 只保留下方和左侧的轴线
    panel.grid.major = element_blank(),  # 去除主要网格线
    panel.grid.minor = element_blank(),   # 去除次要网格线
    legend.position = "none",  # 移除图例
    axis.text.x = element_text(size = 12, face = "bold"),  # 加粗并设置X轴标签字体大小
    axis.text.y = element_text(size = 12),  # 设置Y轴标签字体大小
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

# 打印绘制的图形
print(ZSWIM3_plot)

# 保存图形为PDF文件
ggsave("ZSWIM3.pdf", plot = ZSWIM3_plot, width = 4, height = 4, device = "pdf")



