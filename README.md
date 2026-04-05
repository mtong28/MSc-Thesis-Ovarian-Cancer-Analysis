# Transcriptomic Analysis Pipeline & Multi-Cancer Validation

**MSc Thesis Project: Identification of disease-associated lncRNAs in HGSOC and cross-cancer validation.**

---

## 🔬 Project Overview
This repository (originally developed for my MSc thesis: *Ovarian Cancer Analysis*) contains a **robust bioinformatics pipeline** designed for the end-to-end processing of transcriptomic data. 

While the primary application was **High-Grade Serous Ovarian Cancer (HGSOC)**, the workflow is architected for **pan-cancer extensibility**. I have specifically implemented validation logic to leverage large-scale datasets from **TCGA**, including **Breast Cancer (BRCA)**, to identify shared regulatory axes and prognostic signatures.

---

## 📊 Datasets Handled
- **Discovery Dataset:** Whole-transcriptome RNA-Seq from omental biopsies of 35 metastasised HGSOC patients (Canbuild study).
- **Validation Datasets:** 
  - TCGA Ovarian Cancer Dataset (`TCGA-OV`, n=374).
  - Pan-Cancer TCGA Datasets including Colon (`TCGA-COAD`), Bladder (`TCGA-BLCA`), Lung (`TCGA-LUAD`/`LUSC`), Breast (`TCGA-BRCA`), and others.
- **Experimental Controls:** Primary adipose tissue and HGSOC samples from Professor Balkwill’s laboratory.

---

## 📂 Repository Structure & Workflow

### 1. Upstream Processing & HPC Operations
- **`commands.sh`**
  - Logs of high-performance computing operations on the cluster.
  - Demonstrates competency in handling restricted root permissions (manual local compilation of `HTSeq v0.6.1p1`), environment module swapping (Python 3.8.3 and R 4.0.2), and secure data transfer via SSH port forwarding from a firewalled cluster.
- **`htseq_counts25.py`**
  - Custom Python script used to automate the iterative parsing of individual sample count files, merging them into a unified raw count matrix (`htseq_all_count.csv`).

### 2. Downstream Analysis & Clinical Validation
- **`01_normalization_and_correlation.R`**
  - **Data Wrangling:** Resolved critical sample-name mismatch issues between TCGA clinical metadata and HTSeq output matrices (handling `X` prefixes and `.` vs `-` delimiters).
  - **Normalisation:** Applied `edgeR` and `limma-voom` workflows. Filtered out low-abundance genes present in less than 33% of the cohort.
  - **Differential Expression & Hypothesis Testing:** Implemented linear modeling (`lmFit`) to identify candidate lncRNAs. Conducted Shapiro-Wilk tests for normality, followed by independent T-tests and Pearson’s correlation to validate the relationship between `l_LNC14112` and its host gene `REV3L`.

- **`02_survival_analysis.R`**
  - **Algorithm:** Implements a data-driven approach to identify the optimal prognostic threshold. A custom function `FindBestCutoff_FromVectors` iterates through the 20th to 80th percentiles of the risk scores to fit Cox proportional hazards models and select the cut-off that minimizes the Wald test P-value.
  - **Visualization:** Generates high-quality Kaplan-Meier curves with risk tables using `survminer` and `ggplot2` for `l_LNC14112`, `REV3L`, and their ratio.

---

*Note: Raw sequencing BAM/FASTQ files and controlled-access TCGA datasets are not included in this repository due to size and privacy constraints.*

