# 🧬 基础生信分析代码  
# 🧬 Basic Bioinformatics Analysis Pipeline Scripts

---

## 📖 项目简介 | Project Overview

本项目整理了一个常规的较为基础的 **医学生信分析文章的流程**，涵盖从测序数据格式化，到差异表达分析、功能富集分析、生存分析及免疫浸润分析的常见步骤。  
The repository provides a basic ** Bioinformatics analysis pipeline**, including data preprocessing, differential expression, enrichment analysis, survival modeling, and immune infiltration estimation.

本项目的脚本大部分来源于：
- 开源生信教程与公共代码；
- 本人学习与理解后修改完善；
- 借助 AI（如 ChatGPT）优化代码结构与注释。  

Most scripts are based on:
- Open-source bioinformatics tutorials and public repositories;  
- My own learning, understanding, and improvements;  
- Assistance from AI tools (e.g., ChatGPT) for optimization and documentation.  

🙏 **感谢所有开源项目的贡献者们！**  
希望本项目能帮助生信初学者理解标准分析流程，并继续传承开源精神。  
🙏 **Thanks to all open-source contributors!**  
This project aims to help bioinformatics beginners learn the standard workflow and continue the spirit of open sharing.

---

## 🧩 文件结构与功能说明 | File Structure & Descriptions

### 📁 主体分析脚本 | Main Analysis Scripts

| 文件名 / Script | 功能说明 | Description |
|------------------|------------------------|------------------------|
| **1.1_TCGA_counts_data_Formatting.r** | TCGA counts 数据格式化 | Formatting TCGA counts data |
| **1.2_TCGA_fpkm_data_Formatting.r** | TCGA FPKM 数据格式化 | Formatting TCGA FPKM data |
| **2.1_counts_Differential_gene_analysis.r** | 基于 counts 的差异基因分析 | Differential gene expression (DEG) using counts data |
| **2.2_image_of_DEG.r** | 差异基因可视化（火山图、热图） | DEG visualization (volcano & heatmap) |
| **3.1_Enrichment_analysis.r** | GO/KEGG 富集分析 | Functional enrichment (GO, KEGG) |
| **3.2_GSEA_enrichment_analysis.r** | 基因集富集分析 | GSEA-based enrichment |
| **4.Single_gene_cox.r** | 单基因 Cox 生存分析 | Single-gene Cox regression |
| **5.lasso_cox.r** | LASSO Cox 模型分析 | LASSO Cox regression modeling |
| **6.ROC.r** | 绘制 ROC 曲线评估模型性能 | Plot ROC curves for model evaluation |
| **7.Survival_analysis.r** | 生存曲线与分组分析 | Kaplan-Meier survival analysis |
| **8.0_group_example.r** | 分组示例脚本 | Grouping example (for immune analysis) |
| **8.1_imm_cibersort.r** | 免疫浸润分析：CIBERSORT | Immune infiltration via CIBERSORT |
| **8.2_imm_estimate.r** | 免疫浸润分析：ESTIMATE | Immune scoring via ESTIMATE |
| **8.3_imm_ssGSEA.r** | 免疫浸润分析：ssGSEA | Immune infiltration via ssGSEA |

---

### 📚 reference_data 参考数据文件 | Reference Data Files

除了分析脚本外，项目还包含一组参考数据文件，用于支持部分分析模块的运行。  
In addition to analysis scripts, this repository contains a **`reference_data`** folder with auxiliary files required by certain analysis steps.

| 文件名 / File | 说明 / Description | 说明 / Description |
|----------------|--------------------|--------------------|
| **gene_length_Table.txt** | 基因注释及长度参考文件，用于数据格式化 | Gene annotation and length table used in data formatting |
| **cellMarker.csv** | 免疫细胞 marker 基因集，用于 ssGSEA 分析 | Immune cell marker gene sets for ssGSEA |
| **LM22.txt** | CIBERSORT 的 signature matrix 文件 | Signature matrix file for CIBERSORT |
| **msigdb_v7.0_GMTs/** | 来自 MSigDB v7.0 的多种基因集（GO、KEGG、Reactome、Hallmark 等），用于 GSEA 分析 | Gene sets from MSigDB v7.0 (GO, KEGG, Reactome, Hallmark, etc.) for GSEA analysis |

---

## 🧠 分析流程概览 | Workflow Overview

1. **数据准备与格式化 / Data Preparation**

2. **差异基因分析 / Differential Expression Analysis**

3. **功能富集分析 / Functional Enrichment**

4. **生存与预后分析 / Survival & Prognostic Analysis**

5. **免疫浸润分析 / Immune Infiltration**

---

## 🧩 未来计划 | Future Updates

📅 **即将更新内容 / Upcoming Updates**

- **GEO 数据处理与分析流程 / GEO Data Processing & Analysis Pipeline**  
  将在后续版本中加入基于 GEO 公共数据集的下载、预处理与差异分析模块。  
  Future updates will include scripts for processing GEO datasets, including data download, normalization, and DEG analysis.

---

## 💻 环境依赖 | Environment Requirements

建议使用 R 4.2 或更高版本。  
Recommended: R ≥ 4.2.  

---

🌟 如果觉得有帮助 | If You Found This Useful

欢迎点一个 ⭐️ 支持开源！
Give this repo a ⭐️ if you find it helpful!
