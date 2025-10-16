# 设置工作目录
original_dir <- "/Users/paperz/Desktop/ABenMao生信分析/analysis"
setwd(original_dir)

# 加载R包
# install.packages("glmnet")
# install.packages("survival")
library("glmnet")
library("survival")
library("tidyverse")

# 读取带有生存时间的DEG数据
expr_for_cox=read.table("/Users/paperz/Desktop/ABenMao生信分析/analysis/single_gene_cox/surv.expr.txt",header=T,sep="\t",row.names=1)           
colnames(expr_for_cox)[2] <- 'fustat'
colnames(expr_for_cox)[1] <- 'futime'
rt <- expr_for_cox         
rt$futime=rt$futime/365   

# 读取单因素cox选出的基因
# 这里直接选用单因素cox显著的基因进行演示
# 根据实际情况选择基因
Coxoutput=read.table("/Users/paperz/Desktop/ABenMao生信分析/analysis/single_gene_cox/single_cox_results.txt",header=T,sep="\t",row.names=1)           

# 提取所有P值小于cutoff的基因
pcutoff <- 0.05
topgene <- Coxoutput %>%
  filter(pvalue < pcutoff) %>%    # 筛选所有符合P值条件的行
  arrange(pvalue)                 # 按P值从小到大排序（最显著的在前）

sig_genes <- rownames(topgene) # 提取需要使用的基因

# 检查并保留rt中存在的基因列,并保留前两列
rt <- rt[, c("futime", "fustat", 
             intersect(sig_genes, colnames(rt)))]

# lasso-cox回归分析
# 设置随机种子，保证结果可复现
set.seed(3)

# 构建矩阵
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))

# 拟合模型
fit=glmnet(x, y, family = "cox", maxit = 1000)
cvfit = cv.glmnet(x, y, family="cox", maxit = 1000)
plot(fit, xvar = "lambda", label = TRUE)
plot(cvfit)

# 输出预测模型的相关系数与riskScore
# 输出相关系数
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
geneCoef   #查看模型的相关系数

# 计算riskScore
FinalGeneExp = rt[,lassoGene]
myFun = function(x){crossprod(as.numeric(x),actCoef)}
riskScore = apply(FinalGeneExp,1,myFun)
outCol = c("futime", "fustat", lassoGene)
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
dat = cbind(rt[,outCol], riskScore=as.vector(riskScore), risk)

# # 创建并切换到新工作目录
# output_dir <- "lasso_cox"
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir)
# }
# setwd(output_dir)  # 切换工作目录到目标文件夹
# 
# # 保存分组文件
# save(dat,file = "lasso_roup.Rdata")

# 风险相关结果可视化
# 风险评分箱线图
library(ggpubr)
p <- ggboxplot(dat, x = "fustat", y = "riskScore",
               color = "fustat", palette = "jco",
               add = "jitter")
p <- p + stat_compare_means()   #  Add p-value
p   #得出预测结果

# ggboxplot(dat, x = "fustat", y = "riskScore", 
#           color = "fustat", add = "jitter") +
#   stat_compare_means()

# ROC曲线
library(ROCR)   #使用ROCR包绘制预测模型的ROC曲线
library(caret)
pred <- prediction(dat$riskScore, dat$fustat)
perf <- performance(pred,"tpr","fpr")
AUC <- performance(pred,"auc")   #计算AUC
plot(perf,colorize=FALSE, col="red", print.auc =TRUE) #绘制ROC曲线
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
dev.off()

# 生存时间分布图
rt <- dat
color=as.vector(rt$fustat)
color[color==1]="indianred1"
color[color==0]="lightseagreen"
plot(rt$futime, pch=19,
     xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("indianred1","lightseagreen"),cex=1.2)
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
abline(v=lowLength,lty=2)
dev.off()

# 风险评分分布图
rt <- rt[order(rt$riskScore),] # 排序
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10
plot(line, type="p", pch=20,
     xlab="Patients (increasing risk socre)", ylab="Risk score",
     col=c(rep("lightseagreen",lowLength),rep("indianred1",highLength)) )
abline(h=median(rt$riskScore),v=lowLength,lty=2)
legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("indianred1","lightseagreen"),cex=1.2)
dev.off()