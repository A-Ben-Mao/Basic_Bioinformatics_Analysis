# 设置工作目录
original_dir <- "文件目录"
setwd(original_dir)

# 加载R包
library(survival)

# 使用时间依赖ROC同款整理文件
surv <- read.table("\single_gene_cox\surv.expr.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
surv$OS.time <- surv$OS.time/365*12 # 将生存时间转换为月份

#### 参数设置 ####
aim_gene <- "DKK1"   # 想分析的基因名
k <- 2               # 设置分组数，比如 2=二分法，3=三分位，4=四分位

# 取基因表达向量
gene_expr <- surv[[aim_gene]]

# 计算分位点
probs <- seq(0, 1, length.out = k + 1)
cut_points <- quantile(gene_expr, probs = probs, na.rm = TRUE)

# 先分成 Q1~Qk
surv$quantile_group <- cut(gene_expr,
                           breaks = cut_points,
                           include.lowest = TRUE,
                           labels = paste0("Q", 1:k))

# 只比较最高组(Qk) vs 其他
top_label <- paste0("Q", k)
surv$group <- ifelse(surv$quantile_group == top_label, "High", "Low")
surv$group <- factor(surv$group, levels = c("Low", "High"))

# 检查结果
print(cut_points)
table(surv$group)

# 生存差异检验（对数秩检验）
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

# 拟合生存曲线(Kaplan-Meier)
fit <- survfit(Surv(OS.time, OS)~ group, data = surv)
summary(fit)

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

