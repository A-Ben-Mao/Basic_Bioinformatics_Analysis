# 设置工作目录
original_dir <- "文件目录"
setwd(original_dir)

# 加载R包
# install.packages("tidyverse")
# install.packages("BiocManager")
# library("BiocManager")
# BiocManager::install('clusterProfiler')
# BiocManager::install('org.Hs.eg.db')
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)

# 读取差异基因结果文件，这里不需要筛选差异基因
load("TCGA_DEG.rda文件目录")
DEG <- as.data.frame(res)

# 读取参考基因集，根据研究目的确定
# 一般选择c2.all.v7.0.entrez.gmt (gsea)
# 或 c5.all.v7.0.entrez.gmt (GO)
reference_gmt <- read.gmt("reference_data/msigdb_v7.0_GMTs/c5.all.v7.0.entrez.gmt文件目录")

# 基因ID转换
DEG <- DEG %>% rownames_to_column("Gene")
genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))

# 准备排序基因列表
geneList = DEG[,3]
names(geneList) = as.character(DEG[,'ENTREZID'])
geneList = sort(geneList, decreasing = TRUE)

# 创建并切换到新工作目录
output_dir <- "enrichment_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
setwd(output_dir)  # 切换工作目录到目标文件夹

# 执行GSEA分析
# set.seed(1)
gsea<-GSEA(geneList,TERM2GENE = reference_gmt) #GSEA分析

# 结果保存
gsea_result_df <- as.data.frame(gsea)
write.table(gsea_result_df,file="GSEA_MSigDb_C5_result.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
save(gsea,gsea_result_df,file = "GSEA_deg_SPP1.rda")

# 可视化
# 单通路绘制
library(enrichplot)
gseaplot2(gsea,1,color="red")
gseaplot2(gsea,3,color="red",pvalue_table = T)

# 多通路绘制（第一个参数为通路序号，第二个参数为图像包含部分）
gseaplot2(gsea, geneSetID = c(1,3), subplots = 1:3)
gseaplot2(gsea, geneSetID = 1:3, subplots = 1:2)
gseaplot2(gsea, geneSetID = 1:10, subplots = 1:3)

# 清空环境
dev.off()
rm(list = ls()) # 移除所有对象
gc() # 清理并释放内存
