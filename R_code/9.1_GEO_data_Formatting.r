# GEO数据下载
# 官网：https://www.ncbi.nlm.nih.gov/geo/

# 本示例使用GSE85841
# 多阅读文献，参考科研人员使用过的GEO数据集

# 设置工作目录
original_dir <- "/Users/paperz/Desktop/ABenMao生信分析/analysis"
setwd(original_dir)  # 切换工作目录到目标文件夹

# 加载R包
library(tidyverse)
library(GEOquery)
library(limma) 

# 创建并切换到新工作目录
output_dir <- "GEO_data"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
setwd(output_dir)  # 切换工作目录到目标文件夹

#### 文件获取及格式化 ####
# 下载数据，如果文件夹中有会直接读入
# chooseBioCmirror()
gset = getGEO('GSE85841', destdir=".", AnnotGPL = F, getGPL = F) # 建议直接官网下载
gset[[1]] # 基因集信息

# 样本分组处理
pdata <- pData(gset[[1]])
table(pdata$source_name_ch1) # 查看患者分组
group_list <- ifelse(str_detect(pdata$source_name_ch1, "lung adenocarcinoma"), "tumor",
                     "normal")
group_list = factor(group_list,
                    levels = c("normal","tumor"))

# 表达矩阵处理
# 获得原始矩阵
exp <- exprs(gset[[1]])
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
dev.off()

# 大多数基因的表达量在不同样本中是不变的（即非差异表达基因占多数）。  
# 因此，如果所有样本的整体分布（如中位数）不一致，
# 更可能是技术误差导致的，而非真实的生物学差异。矫正目的是让技术误差最小。

# 数据校正
exp=normalizeBetweenArrays(exp)
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
range(exp)
exp <- log2(exp+1) # 进行log转换，这一步根据原始数据来判断是否需要对数转换
range(exp)
dev.off()

#### 基因ID手动注释 ####
# 读取GPL注释文件
GPL <- read.delim("~/Desktop/ABenMao生信分析/reference_data/GPL20115-26806的副本.txt", row.names=1)

# 数据交集处理
comname <- intersect(rownames(exp),rownames(GPL))
exp <- exp[comname,]
GPL <- GPL[comname,]

# 合并表达数据与注释信息
exp1 <- as.data.frame(exp)
exp1 <- cbind(exp,GPL)
exp1 <- exp1[!duplicated(exp1$GeneSymbol),] # 注意修改列名
rownames(exp1) <- exp1$GeneSymbol
exp1 <- exp1[,(1:16)] # 筛选患者列，注意修改数据
write.table(exp1, file = "GEO_GSE85841.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

