# 设置工作目录
original_dir <- "/Users/paperz/Desktop/ABenMao生信分析/analysis"
setwd(original_dir)

# 加载R包
library(survival)
library(tidyverse)

#使用时间依赖ROC同款整理文件
fpkm <- read.table("/Users/paperz/Desktop/ABenMao生信分析/analysis/formatted_data/TCGA_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
fpkm <- fpkm %>% t() %>% as.data.frame() # 转置

# 计算中位数并分组
med <- median(fpkm$DKK1, na.rm = TRUE)  # 处理缺失值
fpkm$group <- ifelse(fpkm$DKK1 > med, "High", "Low")
fpkm$group <- factor(fpkm$group, levels = c("Low", "High"))

# 创建新数据框
group_df <- data.frame(
  sample = rownames(fpkm),  # 提取行名作为样本ID
  group = fpkm$group,
  row.names = NULL  # 避免将行名作为单独列
)
rownames(group_df) <- group_df[,1]
group_df <- group_df[, -1, drop = FALSE]

# 检查分组情况
table(group_df$Group)

# 保存为RData文件
save(group_df, file = "group_data.RData")

