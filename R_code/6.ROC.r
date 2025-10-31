# 设置工作目录
original_dir <- "工作目录"
setwd(original_dir)

# 加载R包
# install.packages("ROCR")
# install.packages("rms")
library(tidyverse)
library(ROCR)
library(rms)
library(timeROC)
library(survival)

# 创建并切换到新工作目录
output_dir <- "ROC"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
setwd(output_dir)  # 切换工作目录到目标文件夹

# 读取“单基因cox分析”代码中已经保存好的生存数据与表达基因文件
surv_expr <- read.table("single_gene_cox\\surv.expr.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

#### 单基因ROC曲线 ####
# 构建ROC模型
aim_gene <- "DKK1" # 设置你想分析的基因名

# 构建ROC对象并计算性能
pred <- prediction(surv_expr[[aim_gene]], surv_expr$OS)
perf <- performance(pred, "tpr", "fpr") # 计算TPR（真阳性率）和FPR（假阳性率）
auc  <- performance(pred, "auc")@y.values[[1]]

# 绘制ROC曲线
plot(perf,
     col = "red",
     lwd = 2,
     main = paste0("ROC Curve for ", aim_gene, " (AUC = ", round(auc, 3), ")"),
     xlab = "False Positive Rate",
     ylab = "True Positive Rate")
abline(0, 1, lty = 2, col = "gray") # 绘制对角线

#### timeROC ####
# 数据的整理
surv_expr$OS.time <- surv_expr$OS.time/365 # 生存时间单位从"天"转为"年"

# 构建时间依赖 ROC 模型
ROC_time <- timeROC(
  T = surv_expr$OS.time,       # 结局时间
  delta = surv_expr$OS,        # 生存状态 (0=生存, 1=死亡)
  marker = surv_expr[[aim_gene]], # 基因表达量
  cause = 1,
  weighting = "marginal",
  times = c(1, 3, 5),          # 1年、3年、5年
  iid = TRUE
)

# 绘制多时间点 ROC 曲线
plot(ROC_time$FP[,1], ROC_time$TP[,1], type = "l", col = "red", lwd = 2,
     xlab = "False Positive Rate", ylab = "True Positive Rate",
     main = paste0("Time-dependent ROC for ", aim_gene))
lines(ROC_time$FP[,2], ROC_time$TP[,2], col = "green", lwd = 2)
lines(ROC_time$FP[,3], ROC_time$TP[,3], col = "blue", lwd = 2)
abline(0, 1, lty = 2, col = "gray")

# 添加图例
legend("bottomright",
       legend = c(
         paste0("1-year (AUC=", round(ROC_time$AUC[1], 3), ")"),
         paste0("3-year (AUC=", round(ROC_time$AUC[2], 3), ")"),
         paste0("5-year (AUC=", round(ROC_time$AUC[3], 3), ")")
       ),
       col = c("red", "green", "blue"), lty = 1, lwd = 2, bty = "n")

#### 多基因ROC曲线 ####
# 设定你要比较的基因集合
genes <- c("DKK1", "PYCR1", "OTUD1")

# 创建空图
plot(0, 0, type = "n",
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = "False Positive Rate",
     ylab = "True Positive Rate",
     main = "ROC Curves for Multiple Genes")

# 设置颜色
colors <- rainbow(length(genes))

# 存放AUC结果
AUCs <- numeric(length(genes))

# 循环绘制每个基因的ROC
for (i in seq_along(genes)) {
  gene <- genes[i]
  
  # 构建预测对象
  pred <- prediction(surv_expr[[gene]], surv_expr$OS)
  perf <- performance(pred, "tpr", "fpr")
  
  # 计算AUC
  auc <- performance(pred, "auc")@y.values[[1]]
  AUCs[i] <- auc
  
  # 绘制曲线
  lines(perf@x.values[[1]], perf@y.values[[1]],
        col = colors[i], lwd = 2)
}

# 添加对角线
abline(0, 1, lty = 2, col = "gray")

# 添加图例
legend("bottomright",
       legend = paste0(genes, " (AUC=", round(AUCs, 3), ")"),
       col = colors[seq_along(genes)], lty = 1, lwd = 2, bty = "n")

#### 筛选基因表现良好的基因，用于随后的建模 ####
# 获取基因名（排除结局列和生存时间列）
genes <- setdiff(colnames(surv_expr), c("OS", "OS.time"))

# 初始化空向量存储 AUC
AUCs <- numeric(length(genes))

# 循环计算每个基因的 AUC
for (i in seq_along(genes)) {
  gene <- genes[i]
  
  # 构建 ROCR 对象
  pred <- prediction(surv_expr[[gene]], surv_expr$OS)
  auc <- performance(pred, "auc")@y.values[[1]]
  
  AUCs[i] <- auc
}

# 整理成数据框
roc_results <- data.frame(
  Gene = genes,
  AUC = round(AUCs, 3)
)

# 按照 AUC 值降序排序（AUC 越大预测能力越好）
roc_results <- roc_results %>% arrange(desc(AUC))

# 保存到文件
write.table(roc_results, file = "ROC_AUC_results.txt", sep = "\t", row.names = F, quote = F)
