# TCGA数据的下载：
# xena官网：https://xenabrowser.net/datapages/

# 设置工作目录
original_dir <- "/Users/paperz/Desktop/ABenMao生信分析/analysis"
setwd(original_dir)  # 切换工作目录到目标文件夹

#install.packages("tidyverse")
library(tidyverse)

#### counts数据的格式化 ####
# 读取tsv文件
counts_data = read.table(file = '/Users/paperz/Desktop/ABenMao生信分析/data/TCGA-LUAD.star_counts.tsv.gz', sep = '\t', header = TRUE) 

# 提取所有以"_PAR_Y"结尾的行，这里代表的是Y染色体的伪染色体区，其数值一般为0
par_y_rows <- counts_data[grepl("_PAR_Y$", counts_data$Ensembl_ID), ]
print(paste("找到", nrow(par_y_rows), "个_PAR_Y结尾的基因"))
View(par_y_rows)

# 从原始数据中剔除这些行
counts_clean <- counts_data[!grepl("_PAR_Y$", counts_data$Ensembl_ID), ]

# 将第一列ID转换为行名
rownames(counts_clean) <- counts_clean[,1]
counts_clean = counts_clean[,-1]

# 查询数据中患者的类别（01-09是癌症，10-19是正常，20-29是癌旁）
# 即使是肿瘤组织，01-09意义各不相同，比如，01代表原发灶，02代表转移灶
# A、B、C则代表不同的样本处理手段，如：01A：癌症组织，01B：福尔马林浸泡样本
# 若有需要，可进一步查阅官网注释
table(substr(colnames(counts_clean),14,16))

# 取决于研究目的，这里选择01A和11A数据进行提取
counts_clean <- counts_clean[,substr(colnames(counts_clean),14,16)%in% c("01A","11A")]
table(substr(colnames(counts_clean),14,16)) # 查询分别的数目

# 调整行名格式，使其更方便后续分析
rownames(counts_clean) <- substr(rownames(counts_clean),1,15)

# 反转换counts数据
counts_clean <- ceiling(2^(counts_clean)-1)

# 注释counts文件
counts <- counts_clean

# 加载基因注释参考文件,根据研究目的筛选
Ginfo_0 <- read.table(
  "/Users/paperz/Desktop/ABenMao生信分析/reference_data/gene_length_Table.txt",
  sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Ginfo <- Ginfo_0[which(Ginfo_0$genetype == "protein_coding"),] # 只要编码RNA

# 数据预处理
comgene <- intersect(rownames(counts), rownames(Ginfo))  # 取基因交集
counts <- counts[comgene, ]  # 保留共有的基因
Ginfo <- Ginfo[comgene, ]    # 同步注释信息

# 添加基因名称列并去重
counts$Gene <- as.character(Ginfo$genename)  # 新增Gene列（基因符号）
counts <- counts[!duplicated(counts$Gene), ] # 按基因符号去重
rownames(counts) <- counts$Gene              # 将基因符号设为行名
counts <- counts[, -ncol(counts)]            # 删除多余的Gene列

# 文件的输出
# 创建并切换到新工作目录
output_dir <- "formatted_data"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
setwd(output_dir)  # 切换工作目录到目标文件夹

# 输出为文本到formatted_data文件夹
write.table(counts_clean, file = "counts_TCGA.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

# 输出为CSV到formatted_data文件夹
# write.csv(counts_clean, file = "counts_TCGA.csv")

# 保存处理后的完整数据
write.table(counts, file = "TCGA_counts_mRNA_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# 保存癌症患者的counts
tumor <- colnames(counts)[substr(colnames(counts),14,16) == "01A"]
counts_01A <- counts[,tumor]
write.table(counts_01A, file = "TCGA_counts_mRNA_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# 清空环境
rm(list = ls()) # 移除所有对象
gc() # 清理并释放内存

