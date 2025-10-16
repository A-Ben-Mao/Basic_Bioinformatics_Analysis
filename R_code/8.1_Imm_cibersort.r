# 设置工作目录
original_dir <- "/Users/paperz/Desktop/ABenMao生信分析/analysis"
setwd(original_dir)

# 创建并切换到新工作目录
output_dir <- "Imm_cibersort"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
setwd(output_dir)  # 切换工作目录到目标文件夹

# 加载R包
# install.packages('e1071')
# install.packages('parallel')
# #install.packages("BiocManager")
# BiocManager::install("preprocessCore")
library(e1071)
library(parallel)
library(preprocessCore)
library(CIBERSORT)

# 读取相关数据
sig_matrix <- "/Users/paperz/Desktop/ABenMao生信分析/reference_data/LM22.txt" # 22种免疫细胞的特征基因矩阵
mixture_file = '/Users/paperz/Desktop/ABenMao生信分析/analysis/formatted_data/TCGA_fpkm_mRNA_01A.txt'   # 肿瘤患者表达谱

# 运行CIBERSORT
res_cibersort <- cibersort(sig_matrix, mixture_file, perm=100, QN=TRUE)
save(res_cibersort,file = "res_cibersort.Rdata")   #保存中间文件
# 可以用来做很多图片，比如：https://mp.weixin.qq.com/s/vxar-e-JXpFQwutyj-jbyg

# 提取结果（排出得分，细胞比例作图用）
res_cibersort_filter <- res_cibersort[,1:22]   
ciber.res <- res_cibersort_filter[,colSums(res_cibersort_filter) > 0]   #去除丰度全为0的细胞

# 可视化细胞比例
mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7) #创建彩虹色板（带70%透明度）
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(ciber.res)),
        border = NA, # 柱子无边框写
        names.arg = rep("",nrow(ciber.res)), # 无横坐标样本名
        yaxt = "n", # 先不绘制y轴
        ylab = "Relative percentage", # 修改y轴名称
        col = mycol) # 采用彩虹色板
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), # 补齐y轴添加百分号
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-20, # 这里-20要根据实际出图的图例位置情况调整
       par("usr")[4], 
       legend = colnames(ciber.res), 
       xpd = T,
       fill = mycol,
       cex = 0.6, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()

#### 基因表达与免疫细胞比例的相关性 ####
library(corrplot)
library(tidyverse)

# 加载基因表达谱
expr <- read.table(mixture_file,sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

# 提取目标基因的相关数据
gene <- c("DKK1","VAX1","KRT6A","TRPA1",'GJB3',
          'RGS20','IGF2BP1','C1QTNF6','KCNV1',"IGFBP1")
exp <- expr[gene,]
exp <- exp %>% t() %>% as.data.frame()

# 加载cibersort的结果
# load("res_cibersort.Rdata")
ciber <- res_cibersort[,1:22]
ciber <- as.data.frame(ciber)
identical(rownames(ciber),rownames(exp)) # 检查样本一致性

# 算基因表达与免疫细胞比例的Spearman相关系数
cor<-sapply(ciber,function(x,y) cor(x,y,method="spearman"),exp)
rownames(cor)<-colnames(exp)

cor_res <- cor.mtest(cor,#计算p值
                     conf.level = 0.95)#置信区间
corrplot(cor,
         method = "color",#相关性矩阵展示的图形
         col=colorRampPalette(c("#01468b","white","#ee0000"))(100),
         addCoef.col = "black",#为相关系数添加颜色
         tl.col="black",#设置文本标签的颜色
         number.cex = 0.5,
         tl.cex = 0.7,
         cl.align = "l")
dev.off()

#### 分组对CIBERSORT结果进行可视化 ####
# 加载cibersort结果
# load("res_cibersort.Rdata")
a <- res_cibersort[,1:22]
a <- as.data.frame(a)

# 添加分组信息（需要已经进行分组）
load("/Users/paperz/Desktop/ABenMao生信分析/analysis/group_data.RData")
identical(rownames(a),rownames(group_df))
b <- group_df
class(b$group)
a$group <- b$group
a <- a %>% rownames_to_column("sample")

# 加载R包
library(ggsci)
library(tidyr)
library(ggpubr)

b <- gather(a,key=CIBERSORT,value = Proportion,-c(group,sample))

ggboxplot(b, x = "CIBERSORT", y = "Proportion",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 
dev.off()
