# GEO数据下载
# 官网：https://www.ncbi.nlm.nih.gov/geo/

# 本示例使用GSE31210
# 多阅读文献，参考科研人员使用过的GEO数据集

# 设置工作目录
original_dir <- "文件目录"
setwd(original_dir)  # 切换工作目录到目标文件夹

# 加载R包
library(tidyverse)
library(GEOquery)
library(limma) 

#### 文件获取及格式化 ####
# 下载数据，如果文件夹中有会直接读入
# chooseBioCmirror()
gset = getGEO('GSE31210', destdir=".", AnnotGPL = F, getGPL = F) # 建议直接官网下载
gset[[1]] # 基因集信息
pdata <- pData(gset[[1]]) #查看数据集基本信息

# 数据矫正
exp <- exprs(gset[[1]])
exp=normalizeBetweenArrays(exp)
range(exp)
exp <- log2(exp+1) # 进行log转换，这一步根据原始数据来判断是否需要对数转换
range(exp)
exp <- as.data.frame(exp)

#### 基因ID手动注释 ####
index = gset[[1]]@annotation #检查测序平台
# 读取GPL注释文件
GPL <- read.delim("注释文件目录", row.names=1)

# 一些特殊情况
GPL$Gene.Symbol <- sub(" ///.*", "", GPL$Gene.Symbol)  # 删除"///"后的别名

# 数据交集处理
comname <- intersect(rownames(exp),rownames(GPL))
exp <- exp[comname,]
GPL <- GPL[comname,]

# 合并表达数据与注释信息
# 合并表达数据与注释信息
exp1 <- as.data.frame(exp)
exp1 <- cbind(exp,GPL)

# 去除 NA 或空的基因名（注意修改列名）
exp1 <- exp1[!is.na(exp1$Gene.Symbol) & exp1$Gene.Symbol != "", ]

# 统计重复基因数量并打印
dup_genes <- sum(duplicated(exp1$Gene.Symbol))
cat("共有", dup_genes, "个重复的基因名（probe->gene 多对一）。\n")

# 明确样本列（即只选择数值列）
sample_cols <- colnames(exp)# 原始 exp 的样本列名

# 按 Gene.Symbol 对样本列取均值
agg_df <- aggregate(exp1[, sample_cols], 
                    by = list(Gene.Symbol = exp1$Gene.Symbol), 
                    FUN = mean, na.rm = TRUE)
rownames(agg_df) <- agg_df$Gene.Symbol
agg_df$Gene.Symbol <- NULL

# 保存表达矩阵
write.table(agg_df, file = "GEO_GSE31210.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#### 整理生存信息 ####
# 加载R包
library(survival)
library(readr)

# 提取临床信息，先前已经进行过
# pdata <- pData(gset[[1]])

# 保存临床信息后手动整理生存信息
write.table(pdata,"GSE31210_clinic.txt" ,quote=FALSE,col.name=NA,sep="\t")

# 随后再读取已经筛选好的数据
surv <- read_tsv("GSE31210_OS.txt", 
                        locale = locale(encoding = "UTF-16LE"))
surv <- as.data.frame(surv)
rownames(surv) <- surv$Sample # 调整行名
surv <- surv[,-1]
surv$OS.time <- surv$OS.time/365*12 # 将生存时间转换为月份

# 整合基因表达和生存数据
exp2 <- agg_df %>% t() %>% as.data.frame()
common_samples <- intersect(rownames(surv), rownames(exp2)) # 找出两个数据集共有的样本
surv <- surv[common_samples, , drop = FALSE]
exp2 <- exp2[common_samples, , drop = FALSE]
exp2 <- exp2[rownames(surv), , drop = FALSE]  # 确保行名顺序一致
surv_exp <- cbind(surv, exp2) # 合并数据

#### 生存分析参数设置 ####
aim_gene <- "DKK1"   # 想分析的基因名
k <- 2               # 设置分组数，比如 2=二分法，3=三分位，4=四分位

# 取基因表达向量
gene_expr <- surv_exp[[aim_gene]]

# 计算分位点
probs <- seq(0, 1, length.out = k + 1)
cut_points <- quantile(gene_expr, probs = probs, na.rm = TRUE)

# 先分成 Q1~Qk
surv_exp$quantile_group <- cut(gene_expr,
                               breaks = cut_points,
                               include.lowest = TRUE,
                               labels = paste0("Q", 1:k))

# 只比较最高组(Qk) vs 其他
top_label <- paste0("Q", k)
surv_exp$group <- ifelse(surv_exp$quantile_group == top_label, "High", "Low")
surv_exp$group <- factor(surv_exp$group, levels = c("Low", "High"))

# 检查结果
print(cut_points)
table(surv_exp$group)

# 保存分组（如果需要），一般用于后续分析
group_df <- surv_exp
save(group_df, file = "surv_exp.expr.group.RData")

# 生存差异检验（对数秩检验）
fitd <- surv_expdiff(surv_exp(OS.time, OS) ~ group,
                     data      = surv_exp,
                     na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

# 拟合生存曲线(Kaplan-Meier)
fit <- surv_expfit(surv_exp(OS.time, OS)~ group, data = surv_exp)
summary(fit)

# 绘制生存曲线
#### 绘制生存曲线 ####
# 基础图像
plot(fit, conf.int = T,
     col = c("blue", "red"),
     lwd = 2,
     xlab = "Time(Months)",
     ylab = "Survival probablity"
)

#添加标签
legend("topright",
       title = "Group",
       c("Low", "High"),
       lwd = 2, lty = 1,
       col = c("blue", "red"))

#添加P值
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))
text(25, 0.2, p.lab)

# 图像2
cols <- c("#1F77B4", "#D62728")  # 蓝/红

# 计算最大随访时间和最小生存率，用于自动放置p值
x_pos <- max(surv$OS.time, na.rm = TRUE) * 0.7
y_pos <- min(fit$surv, na.rm = TRUE) + 0.05

# 绘图
plot(fit,
     conf.int = FALSE,          # 不画置信区间，图会更干净
     col = cols,
     lwd = 2.2,
     mark.time = TRUE,          # 显示删失点
     xlab = "Time (Months)",
     ylab = "Survival probability (%)",
     xlim = c(0, max(surv$OS.time, na.rm = TRUE)),
     ylim = c(0, 1),
     axes = FALSE,
     main = paste0("Kaplan-Meier survival: ", aim_gene)
)

# 添加坐标轴（让纵轴显示百分比）
axis(1)
axis(2, at = seq(0, 1, 0.2), labels = paste0(seq(0, 100, 20), "%"), las = 1)
box()

# 添加图例
legend("topright",
       legend = levels(surv$group),
       col = cols,
       lwd = 2,
       bty = "n",
       title = "Group",
       cex = 0.9)

# 添加P值
y_pos <- max(0, min(fit$surv, na.rm = TRUE) - 0.05)
x_pos <- max(surv$OS.time, na.rm = TRUE) * 0.75

p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001",
                            paste0(" = ", round(pValue, 3))))
text(x_pos, y_pos, p.lab, cex = 1.1, font = 2)

# 图像3
library(survminer)
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间
           risk.table = TRUE, # 显示风险表
           risk.table.col = "strata",
           palette = "jco", # 配色采用jco
           legend.labs = c("Low", "High"), # 图例
           size = 1,
           xlim = c(0,120), # x轴长度，一般为0-10年
           break.time.by = 20, # x轴步长为20个月
           legend.title = "",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Months)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
