# 设置工作目录
original_dir <- "文件目录"
setwd(original_dir)  # 切换工作目录到目标文件夹

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

# 差异基因火山图
library(ggplot2)
library(ggrepel)

# 构造用于绘图的数据框
volcano_df <- deg
volcano_df$gene <- rownames(volcano_df)
volcano_df$negLog10P <- -log10(volcano_df$P.Value + 1e-300)  # 防止P.Value==0导致Inf

# 分别挑选上调和下调基因
top_n <- 10  # 上调和下调各取前10个
label_up <- volcano_df %>%
  filter(change == "up") %>%
  arrange(P.Value) %>%
  slice(1:top_n) %>%
  pull(gene)

label_down <- volcano_df %>%
  filter(change == "down") %>%
  arrange(P.Value) %>%
  slice(1:top_n) %>%
  pull(gene)

label_genes <- c(label_up, label_down)

# 绘图并保存
outfile <- "GEO_volcano.png"
png(filename = outfile, width = 2000, height = 1600, res = 300)
p <- ggplot(volcano_df, aes(x = logFC, y = negLog10P)) +
  geom_point(aes(color = change), alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "stable" = "grey70")) +
  geom_vline(xintercept = c(-logFC, logFC), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(P.Value), linetype = "dashed", color = "black") +
  labs(title = "Volcano plot",
       subtitle = paste0("Thresholds: |logFC| >", logFC, ", P.Value <", P.Value),
       x = "log2 Fold Change",
       y = "-log10(P.Value)",
       color = "") +
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

p <- p + geom_text_repel(data = subset(volcano_df, gene %in% label_genes),
                         aes(label = gene),
                         size = 3,
                         max.overlaps = 30)

print(p)
dev.off()
cat("火山图已保存到：", outfile, "\n")
