# 设置工作目录
original_dir <- "文件目录"
setwd(original_dir)  # 切换工作目录到目标文件夹

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

#### GO分析 ####
go <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)
go_res <- go@result # 提取GO分析结果

# 保存GO分析结果
save(go,go_res,file = "GO_analysis.Rdata")

# 可视化（可调整Top数量）
# 柱状图
barplot(go, showCategory = 10,color = "pvalue")

# 气泡图
dotplot(go, showCategory = 10)

# 分类展示
barplot(go, drop = TRUE, showCategory =10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')
dotplot(go,showCategory = 10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')

#### KEGG ####
kegg <- enrichKEGG(gene         = DEG$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.1,
                 qvalueCutoff =0.1)
kegg_res <- kegg@result

# 保存KEGG分析结果
save(kegg,kegg_res,file = "KEGG_analysis.Rdata")

# 可视化（可调整Top数量）
# 柱状图
barplot(kegg, showCategory = 20,color = "pvalue")
# 气泡图
dotplot(kegg, showCategory = 20)

# 等等，与TCGA的分析大致相同
