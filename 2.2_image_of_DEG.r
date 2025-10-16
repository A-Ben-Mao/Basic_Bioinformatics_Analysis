# 设置工作目录
original_dir <- "/Users/paperz/Desktop/ABenMao生信分析/analysis"
setwd(original_dir)

# 加载R包
library(tidyverse)
library(pheatmap)
library(ggpubr)
library(ggthemes)

# 读取fpkm文件
exp <- read.table("/Users/paperz/Desktop/ABenMao生信分析/analysis/formatted_data/TCGA_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

# 读取差异基因结果文件
load("/Users/paperz/Desktop/ABenMao生信分析/analysis/differential_gene_analysis/TCGA_DEG.rda")

# 差异基因筛选
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 0, padj < 0.05)

# 基因上下调分类
logFC_cutoff <- 1
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)

#### TCGA差异分析热图 ####
# 热图数据准备
cg = rownames(DEG)[DEG$change !="NOT"]
exp_diff <- exp[cg,]

# 样本分组注释
group_list=factor(ifelse(substr(colnames(exp),14,16) == "01A","T","N"),levels = c("N","T"))
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(exp_diff)

# 创建并切换到新工作目录
output_dir <- "differential_gene_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
setwd(output_dir)  # 切换工作目录到目标文件夹

# 绘制热图
pheatmap(exp_diff,
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =T, #是否聚类
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)
dev.off()

#### TCGA差异分析火山图 ####
# padj值对数转换
DEG$logP <- -log10(DEG$padj)

# 绘制普通火山图
ggscatter(DEG, x = "log2FoldChange", y = "logP", xlab = "log2FoldChange",
          ylab = "-log10(Adjust P-value)",
          color = "change",
          palette = c("blue", "gray", "red"),
          size = 1) +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")

dev.off()

# 添加Top基因信息
DEG$Label = ""   #新加一列label
DEG <- DEG[order(DEG$padj), ]   #对差异基因的p值进行从小到大的排序
DEG$Gene <- rownames(DEG)

# 选取上下调各5个最显著基因
up.genes <- head(DEG$Gene[DEG$change == "UP"], 5)
down.genes <- head(DEG$Gene[DEG$change == "DOWN"], 5)
DEG.top5.genes <- c(up.genes, down.genes)
DEG$Label[match(DEG.top5.genes, DEG$Gene)] <- DEG.top5.genes

# 绘制标注Top基因的火山图
ggscatter(DEG, x = "log2FoldChange", y = "logP",
          color = "change",
          palette = c("blue", "gray", "red"),
          size = 1,
          label = DEG$Label,
          font.label = 8,
          repel = T,
          xlab = "log2FoldChange",
          ylab = "-log10(Adjust P-value)") +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")

dev.off()

# 清空环境
rm(list = ls()) # 移除所有对象
gc() # 清理并释放内存
