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

# 读取差异基因结果文件
load("TCGA_DEG.rda文件目录")

# 差异基因筛选
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05)

# 基因ID转换
DEG <- DEG %>% rownames_to_column("Gene")
genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))

# 创建并切换到新工作目录
output_dir <- "enrichment_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
setwd(output_dir)  # 切换工作目录到目标文件夹

#### GO分析 ####
go <- enrichGO(
  gene = DEG$ENTREZID,          # 输入基因列表（ENTREZ ID）
  OrgDb = org.Hs.eg.db,         # 指定物种数据库
  ont = "all",                  # 同时分析BP/MF/CC三个类别
  pAdjustMethod = "BH",         # p值校正方法（Benjamini-Hochberg）
  minGSSize = 5,                # 最小基因集大小（通常设为5-10）
  maxGSSize = 500,              # 最大基因集大小
  pvalueCutoff = 0.05,          # p值阈值
  qvalueCutoff = 0.05,          # q值阈值（FDR校正后）
  readable = TRUE               # 结果中显示基因Symbol而非ENTREZ ID
)

# 提取结果
go_res <- go@result

# 保存GO分析结果
save(go,go_res,file = "GO_analysis.Rdata")

# 可视化（可调整Top数量）
# 柱状图
barplot(go, showCategory = 20,color = "pvalue")

# 气泡图
dotplot(go, showCategory = 20)

# 分类展示
barplot(go, drop = TRUE, showCategory =10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')
dotplot(go,showCategory = 10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')

#### KEGG分析 ####
kegg <- enrichKEGG(gene         = DEG$ENTREZID,
                   organism     = 'hsa',
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.1)

# 提取结果
kegg_res <- kegg@result

# 创建并切换到新工作目录
output_dir <- "enrichment_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
setwd(output_dir)  # 切换工作目录到目标文件夹

# 保存KEGG分析结果
save(kegg,kegg_res,file = "KEGG_analysis.Rdata")

# 可视化（可调整Top数量）
# 柱状图
barplot(kegg, showCategory = 20,color = "pvalue")
# 气泡图
dotplot(kegg, showCategory = 20)

# 清空环境
dev.off()
rm(list = ls()) # 移除所有对象
gc() # 清理并释放内存
