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

# 创建并切换到新工作目录
output_dir <- "GEO_data"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
setwd(output_dir)  # 切换工作目录到目标文件夹

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
exp1 <- as.data.frame(exp)
exp1 <- cbind(exp,GPL)
exp1 <- exp1[!duplicated(exp1$Gene.Symbol),] # 注意修改列名
rownames(exp1) <- exp1$Gene.Symbol # 注意修改列名
exp1 <- exp1[,(1:246)] # 筛选患者列，注意修改数据
write.table(exp1, file = "GEO_GSE31210.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#### 生存分析 ####
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
exp2 <- exp1 %>% t() %>% as.data.frame()
common_samples <- intersect(rownames(surv), rownames(exp2)) # 找出两个数据集共有的样本
surv <- surv[common_samples, , drop = FALSE]
exp2 <- exp2[common_samples, , drop = FALSE]
exp2 <- exp2[rownames(surv), , drop = FALSE]  # 确保行名顺序一致
surv_exp <- cbind(surv, exp2) # 合并数据

# 计算中位数并分组
median(surv_exp$DKK1) # 注意修改基因
surv_exp$group <- ifelse(surv_exp$DKK1 > median(surv_exp$DKK1),"High","Low")
surv_exp$group <- factor(surv_exp$group, levels = c("Low","High")) 
class(surv_exp$group)
table(surv_exp$group)

# 生存差异检验（对数秩检验）
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = surv_exp,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

# 拟合生存曲线(Kaplan-Meier)
fit <- survfit(Surv(OS.time, OS)~ group, data = surv_exp)
summary(fit)

# 绘制生存曲线
# 基础图像(可选)
###3. 设置颜色，坐标
plot(fit, conf.int = T,
     col = c("blue", "red"),
     lwd = 2,
     xlab = "Time(Months)",
     ylab = "Survival probablity(%)"
)
###添加标签
legend("topright",
       title = "Group",
       c("Low", "High"),
       lwd = 2, lty = 1,
       col = c("blue", "red"))
###添加P值
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))
text(25, 0.2, p.lab)
dev.off()


# 更好看的图像(可选)
library(survminer)
ggsurvplot(fit,
           data = surv_exp,
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
dev.off()
