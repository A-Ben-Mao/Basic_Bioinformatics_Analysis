# TCGA生存数据的下载：
# xena官网：https://xenabrowser.net/datapages/

# 设置工作目录
original_dir <- "文件目录"
setwd(original_dir)

# 安装加载R包
# install.packages("survival")
# install.packages("forestplot")
library(survival)
library(forestplot)
library(tidyverse)

# 读取生存信息tsv文件
surv = read.table(file = 'TCGA-LUAD.survival.tsv.gz文件目录', sep = '\t', header = TRUE) 

# 整理生存信息数据
surv$sample <- gsub("-",".",surv$sample)
rownames(surv) <- surv$sample
surv <- surv[,-1]
surv <- surv[,-3]

# 读取表达数据
expr <- read.table("TCGA_fpkm_mRNA_all.txt文件目录",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

# 匹配样本数据及格式
comgene <- intersect(colnames(expr),rownames(surv))
table(substr(comgene,14,16)) # 查询对应类别数目
expr <- expr[,comgene]
surv <- surv[comgene,]

# 读取差异基因结果文件
load("TCGA_DEG.rda文件目录")

# 目的是挑选要进行单因素cox的基因，可根据目的修改
res_deseq2 <- as.data.frame(res)%>%
  arrange(padj) %>%
  dplyr::filter(abs(log2FoldChange) > 2, padj < 0.05)

# 整合
deg_expr <- expr[rownames(res_deseq2),] %>% t() %>% as.data.frame() #挑选expr中之前筛选出来的差异基因，随后转置
surv.expr <- cbind(surv,deg_expr) #合并生存时间的数据

# 单基因Cox分析
Coxoutput <- NULL 
for(i in 3:ncol(surv.expr)){
  g <- colnames(surv.expr)[i]
  cox <- coxph(Surv(OS.time,OS) ~ surv.expr[,i], data = surv.expr) # 单变量cox模型
  coxSummary = summary(cox)
  
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                           lower = as.numeric(coxSummary$conf.int[,3][1]),
                                           upper = as.numeric(coxSummary$conf.int[,4][1]),
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}

# 创建并切换到新工作目录
output_dir <- "single_gene_cox"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
setwd(output_dir)  # 切换工作目录到目标文件夹

# 保存结果文件
write.table(surv.expr, file = "surv.expr.txt",sep = "\t",row.names = T,col.names = T,quote = F)
write.table(Coxoutput, file = "single_cox_results.txt",sep = "\t",row.names = F,col.names = T,quote = F)

# 筛选显著基因
pcutoff <- 0.05

# 提取所有P值小于cutoff的基因，并按P值升序排序
# topgene <- Coxoutput %>%
#   filter(pvalue < pcutoff) %>%    # 筛选所有符合P值条件的行
#   arrange(pvalue)                 # 按P值从小到大排序（最显著的在前）

# 提取满足条件的Top基因，并按P值升序排序
topgene <- Coxoutput %>%
  filter(pvalue < pcutoff) %>%    # 筛选所有符合P值条件的行
  arrange(pvalue) %>%             # 按P值从小到大排序（最显著的在前）
  slice_head(n = 10)              # 取前10个最显著的基因


# 可视化
# 输入表格的制作
tabletext <- cbind(c("Gene",topgene$gene),
                   c("HR",format(round(as.numeric(topgene$HR),3),nsmall = 3)),
                   c("lower 95%CI",format(round(as.numeric(topgene$lower),3),nsmall = 3)),
                   c("upper 95%CI",format(round(as.numeric(topgene$upper),3),nsmall = 3)),
                   c("pvalue",format(round(as.numeric(topgene$p),3),nsmall = 3)))
# 绘制森林图
total_rows <- 1 + nrow(topgene)  # 计算总行数

hrzl_lines <- setNames(
  list(gpar(lwd=2, col="black"),
       gpar(lwd=1.5, col="black"),
       gpar(lwd=2, col="black")),
  c("1", "2", as.character(total_rows + 1)) # 动态制表
)

forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(topgene$HR)),
           lower=c(NA,as.numeric(topgene$lower)), 
           upper=c(NA,as.numeric(topgene$upper)),
           graph.pos=5,# 图在表中的列位置
           graphwidth = unit(.25,"npc"),# 图在表中的宽度比
           fn.ci_norm="fpDrawDiamondCI",# box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),# box颜色
           
           boxsize=0.4,# box大小固定
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=T,# 显示区间
           zero=1,# zero线横坐标
           lwd.zero=1.5,# zero线宽
           xticks = c(0.5,1,1.5),# 横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab="Hazard ratios",
           txt_gp=fpTxtGp(label=gpar(cex=1.2),# 各种字体大小设置
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           hrzl_lines = hrzl_lines, # 在最后一行上画黑色实线
           lineheight = unit(.75,"cm"),# 固定行高
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
# 清空环境
dev.off()
rm(list = ls()) # 移除所有对象
gc() # 清理并释放内存

