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
library(GEOquery)
library(limma)

# 这一步和之前一样，只是为了获得分组信息演示
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

# 读取经过格式化的表达数据
exp <- read.table("GEO_data/GEO_GSE85841.txt文件目录",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

# 差异分析
design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
write.table(deg, file = "GEO_deg_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# 标记上下调基因
logFC=1
P.Value = 0.05
k1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)
k2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable")) #增加新的列
table(deg$change)

# 差异基因热图
cg = rownames(deg)[deg$change !="stable"]
diff=exp[cg,]
library(pheatmap)
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(diff) 
pheatmap(diff,
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)
dev.off()

