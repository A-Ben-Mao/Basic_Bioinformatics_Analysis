# ssGSEA，全称为单样本基因集富集分析（single-sample GSEA），
# 本质上是将每一个样本的表达谱，转换为对预定义基因集的富集打分，
# 从而反映该样本中某种生物通路、免疫细胞的相对活性或丰度。
# 因此借助免疫细胞的marker基因集，也可实现评估免疫细胞浸润的目的
# 此分析方法可用于非肿瘤分析

# 设置工作目录
original_dir <- "文件目录"
setwd(original_dir)

# 加载R包
library(tidyverse)
library(data.table)
library(GSVA)

# 创建并切换到新工作目录
output_dir <- "Imm_ssGSEA"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
setwd(output_dir)  # 切换工作目录到目标文件夹

# 读取细胞标记基因文件
cellMarker <- data.table::fread("reference_data/cellMarker.csv文件目录",data.table = F)
colnames(cellMarker)[2] <- "celltype"

# 预处理细胞标记基因
type <- split(cellMarker, cellMarker$celltype)  # 按细胞类型拆分成列表
cellMarker <- lapply(type, function(x){
  dd = x$Metagene
  unique(dd)
})
save(cellMarker, file = "cellMarker_ssGSEA.Rdata") # 保存处理好的文件

# 处理基因表达矩阵
expr <- data.table::fread("TCGA_fpkm_mRNA_01A.txt文件目录",data.table = F)   #读取表达文件
rownames(expr) <- expr[,1]   #将第一列作为行名
expr <- expr[,-1]            #去除第一列
expr <- as.matrix(expr)      #将expr转换为矩阵格式

# ssGSEA量化免疫浸润
# gsva_data <- gsva(expr,cellMarker, method = "ssgsea")
gsvaPar <- ssgseaParam(exprData = expr, 
                       geneSets = cellMarker,
                       normalize = TRUE)
gsva_data <- gsva(gsvaPar, verbose = FALSE)

# 结果整理
a <- gsva_data %>% t() %>% as.data.frame()

# 添加分组信息（需要已经进行分组）
load("group_data.RData文件目录")
identical(rownames(a),rownames(group_df))
a$group <- group_df$group
a <- a %>% rownames_to_column("sample")

# 保存结果
write.table(a,"ssGSEA.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# 可视化
library(ggsci)
library(tidyr)
library(ggpubr)

# 转换数据为长格式
b <- gather(a,key=ssGSEA,value = Expression,-c(group,sample))

# 绘制分组箱线图
ggboxplot(b, x = "ssGSEA", y = "Expression",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 

dev.off()
