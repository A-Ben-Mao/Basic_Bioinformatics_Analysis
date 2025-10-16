# 一般使用fpkm或者tpm值，此二者可以互相缓缓

# 设置工作目录
original_dir <- "/Users/paperz/Desktop/ABenMao生信分析/analysis"
setwd(original_dir)

# 加载R包
library(utils)
library(estimate)
library(tidyverse)

# 创建并切换到新工作目录
output_dir <- "Imm_estimate"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
setwd(output_dir)  # 切换工作目录到目标文件夹

# 读取肿瘤患者01A表达谱
expr <- read.table("/Users/paperz/Desktop/ABenMao生信分析/analysis/formatted_data/TCGA_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
input_file <- "TCGA_fpkm_mRNA_01A.txt"
write.table(expr, file = input_file, sep = "\t", quote = F, col.names = NA)

#### 计算免疫评分 ####
# 转换格式
# 筛选ESTIMATE算法所需的1041个特征基因
# 转换数据为GCT格式（基因表达矩阵标准格式）
filterCommonGenes(input.f = "TCGA_fpkm_mRNA_01A.txt",   #输入文件名
                  output.f = "TCGA_fpkm_mRNA_01A.gct",   #输出文件名
                  id = "GeneSymbol")   #行名为gene symbol

# 计算ESTIMATE评分
# StromalScore：间质评分
# ImmuneScore：免疫评分（核心目标）
# ESTIMATEScore：综合评分
# TumorPurity：肿瘤纯度
estimateScore("TCGA_fpkm_mRNA_01A.gct",   #刚才的输出文件名
              "TCGA_fpkm_mRNA_01A_estimate_score.txt",   #新的输出文件名（即估计的结果文件）
              platform="affymetrix")   #默认平台

# 读取计算结果
result <- read.table("TCGA_fpkm_mRNA_01A_estimate_score.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
result <- result[,-1]  #去除首列重复的内容
colnames(result) <- result[1,]  #将首行转变为行title
result <- as.data.frame(t(result[-1,])) #转置

# ImmuneScore越高表示免疫细胞浸润越多
# TumorPurity越高表示肿瘤细胞比例越高

# 保存最终结果
write.table(result, file = "LIHC_fpkm_mRNA_01A_estimate_score.txt",sep = "\t",row.names = T,col.names = NA,quote = F) # 保存并覆盖得分

# 可以进行分组比较，可以评估单组肿瘤纯度 etc.