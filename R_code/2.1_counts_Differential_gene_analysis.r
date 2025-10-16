# 有原始 Counts 数据 -> 可选用DESeq2。
# 只有 FPKM/TPM 数据 -> 使用limma包（需要先进行 log2 转换）。

# 设置工作目录
original_dir <- "/Users/paperz/Desktop/ABenMao生信分析/analysis"
setwd(original_dir)

# 差异表达分析（DESeq2）
library(tidyverse)
library(DESeq2)

# 读取处理好的counts数据
counts <- read.table("/Users/paperz/Desktop/ABenMao生信分析/analysis/formatted_data/TCGA_counts_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

# 差异基因分析（对处理后的counts文件进行分析）
# 过滤低表达基因
counts = counts[apply(counts, 1, function(x) sum(x > 1) > 32), ]

# 创建分组信息，这里可根据情况修改
conditions <- data.frame(
  sample = colnames(counts),
  group = factor(ifelse(substr(colnames(counts), 14, 16) == "01A", "T", "N"),
                 levels = c("N", "T"))
) %>% 
  column_to_rownames("sample")

# 构建DESeq2对象
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = conditions,
  design = ~ group
)

# 创建并切换到新工作目录
output_dir <- "differential_gene_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
setwd(output_dir)  # 切换工作目录到目标文件夹

# 运行差异分析
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
save(res,file = "TCGA_DEG.rda")

# 筛选显著差异基因
# res_deseq2 <- as.data.frame(res)%>% 
#   arrange(padj) %>% 
#   dplyr::filter(abs(log2FoldChange) > 3, padj < 0.05)#根据自己需要

# 清空环境
rm(list = ls()) # 移除所有对象
gc() # 清理并释放内存

