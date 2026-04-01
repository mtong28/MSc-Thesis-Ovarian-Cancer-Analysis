# MSc-Thesis-Ovarian-Cancer-Analysis

# Transcriptomic Analysis of Ovarian Cancer (MSc Thesis)

This repository contains core scripts used for my MSc thesis at Queen Mary University of London, focusing on the transcriptomic analysis of ovarian cancer patients.

## 🔬 Project Overview
- **Objective:** To identify non-coding RNA isoforms associated with REV3L and analyze their clinical relevance.
- **Data Source:** TCGA datasets (Ovarian, Colon, Breast, Lung) and clinical metadata.

## 💻 Workflow & Tools
- **Read Alignment:** `Hisat2`
- **Transcript Assembly:** `StringTie`
- **Quantification:** `HTSeq`
- **Differential Expression:** `DESeq2` (R)
- **Downstream Analysis:** Survival analysis, correlation analysis, and data visualization in R.

## 📂 Repository Structure
- `DE_analysis.R`: R script for DESeq2 and survival plots.
- `alignment_pipeline.sh`: Shell script for read alignment and assembly.

*Note: Raw sequencing data are not included due to privacy and file size constraints.*
