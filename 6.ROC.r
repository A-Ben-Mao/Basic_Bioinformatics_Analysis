# 设置工作目录
original_dir <- "/Users/paperz/Desktop/ABenMao生信分析/analysis"
setwd(original_dir)

# 加载R包
library(tidyverse)

# 读取生存信息
surv = read.table(file = '/Users/paperz/Desktop/ABenMao生信分析/data/TCGA-LUAD.survival.tsv.gz', sep = '\t', header = TRUE) 

# 生存信息数据格式调整
surv$sample <- gsub("-",".",surv$sample)
rownames(surv) <- surv$sample
surv <- surv[,-1]
surv <- surv[,-3]

# 创建并切换到新工作目录
output_dir <- "ROC"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
setwd(output_dir)  # 切换工作目录到目标文件夹

# 保存整理好的生存信息
write.table(surv, file = "survival.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# 读取表达数据
expr <- read.table("/Users/paperz/Desktop/ABenMao生信分析/analysis/formatted_data/TCGA_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

# 取生存数据和表达数据的交集样本
comgene <- intersect(colnames(expr),rownames(surv))
table(substr(comgene,14,16))
expr <- expr[,comgene]
surv <- surv[comgene,]

# 选择目标基因，并提取数据
gene <- c("DKK1","VAX1","KRT6A","TRPA1",'GJB3',
          'RGS20','IGF2BP1','C1QTNF6','KCNV1',"IGFBP1")
exp10 <- expr[gene,] %>% t() %>% as.data.frame()

# 整合表达谱与生存信息
exp_sur <- cbind(exp10,surv)
write.table(exp_sur, file = "exp_sur_ROC.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#准备R包
# install.packages("ROCR")
# install.packages("rms")
library(ROCR)
library(rms)

#### 单基因ROC曲线 ####
# 构建ROC模型(这里可以修改为循环代码)
ROC1 <- prediction(exp_sur$DKK1,exp_sur$OS)   # 构建ROC预测模型，这里只选用了其中一个作为示例的单基因ROC分析
ROC2 <- performance(ROC1,"tpr","fpr")   # 计算TPR（真阳性率）和FPR（假阳性率）
AUC <- performance(ROC1, "auc")@y.values[[1]]  # 提取AUC值

#1.4 绘制ROC曲线
plot(ROC2,
     col="red",   #曲线的颜色
     xlab="False positive rate", ylab="True positive rate",   #x轴和y轴的名称
     lty=1,lwd=3,
     main=paste("AUC=",AUC))
abline(0, 1, lty=2, lwd=3)   #绘制对角线
dev.off()

#### timeROC ####
# 加载R包
library(timeROC)
library(survival)
library(tidyverse)

# 数据的整理
exp_sur$OS.time <- exp_sur$OS.time/365 # 生存时间单位从"天"转为"年"
exp_sur_01A <- exp_sur[substr(rownames(exp_sur),14,16) == "01A",] # 只保留肿瘤患者
write.table(exp_sur_01A, file = "exp_sur_ROC_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# 构建ROC曲线函数(注意修改基因)
ROC3 <- timeROC(T=exp_sur_01A$OS.time,   #结局时间
                delta=exp_sur_01A$OS,   #结局指标
                marker=exp_sur_01A$DKK1,   #预测变量
                cause=1,   #阳性结局指标数值
                weighting="marginal",   #计算方法，默认为marginal
                times=c(1, 3, 5),   #时间点，选取1年，3年和5年的生存率
                iid=TRUE)
ROC3   #查看模型变量信息

# 绘制多时间点ROC曲线
plot(ROC3,
     time=1, col="red")   #time是时间点，col是线条颜色
plot(ROC3,
     time=3, col="green", add=TRUE)   #add指是否添加在上一张图中
plot(ROC3,
     time=5, col="blue", add=TRUE)
legend("bottomright",
       c("Year-1", "Year-3", "Year-5"),
       col=c("red", "green", "blue"),
       lty=1, lwd=2)   #添加标签信息

dev.off()