# 设置工作目录
original_dir <- "文件目录"
setwd(original_dir)  # 切换工作目录到目标文件夹

# 创建并切换到新工作目录
output_dir <- "GEO_data"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
setwd(output_dir)  # 切换工作目录到目标文件夹

# 加载R包
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)

# 读取差异基因文件
deg <- read.table("GEO_deg_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

# 筛选差异基因
logFC=1
P.Value = 0.05
k1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)
k2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)
deg <- deg %>% filter(change!="stable")

# 基因ID转换
DEG <- deg
DEG <- DEG %>% rownames_to_column("Gene")

genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))

# GO分析
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)
go_res <- ego@result # 提取GO分析结果

# KEGG
kk <- enrichKEGG(gene         = DEG$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.1,
                 qvalueCutoff =0.1)
kegg_res <- kk@result

# GSEA
msigdb_GMTs <- "msigdb_v7.0_GMTs"
msigdb <- "c5.all.v7.0.entrez.gmt"    #c2.all.v7.0.entrez.gmt 或 c5.all.v7.0.entrez.gmt
#读取上面指定的gmt文件
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))

geneList = DEG[,2]
names(geneList) = as.character(DEG[,'ENTREZID'])
head(geneList)
geneList = sort(geneList, decreasing = TRUE)

set.seed(1)
KEGG<-GSEA(geneList,TERM2GENE = kegmt) #GSEA分析

# 等等，与TCGA的分析大致相同
