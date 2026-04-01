# MSc-Thesis-Ovarian-Cancer-Analysis (Transcriptomic Analysis)

This repository contains the core workflow and R scripts used for my MSc thesis at Queen Mary University of London, focusing on the identification and validation of disease-associated lncRNAs in High-Grade Serous Ovarian Cancer (HGSOC) and other cancer types.

## 🔬 Project Overview
The project involves processing whole-transcriptome bulk RNA-Seq data from clinical biopsies (Canbuild study) and public databases (TCGA) to discover lncRNAs correlated with disease severity, tumor purity, and patient survival.

---

## 📊 Datasets Handled
- **Discovery Dataset:** Whole-transcriptome RNA-Seq from omental biopsies of 35 metastasised HGSOC patients (Canbuild study).
- **Validation Datasets:** 
  - TCGA Ovarian Cancer Dataset (`TCGA-OV`, n=374).
  - Pan-Cancer TCGA Datasets including Colon (`TCGA-COAD`), Bladder (`TCGA-BLCA`), Lung (`TCGA-LUAD`/`LUSC`), Breast (`TCGA-BRCA`), and others.
- **Experimental Controls:** Primary adipose tissue and HGSOC samples from Professor Balkwill’s laboratory.

---

## 💻 Bioinformatics Workflow & Tools

### 1. Upstream Processing (Alignment & Assembly)
- **Read Alignment:** `Hisat2 (v2.1.0)` mapped to Gencode v29 / GRCh38.p12 reference genome.
- **De novo Transcript Assembly:** `StringTie (v1.3.5)` used to merge and assemble transcripts into a master GTF.

### 2. LncRNA Filtering Criteria
To identify high-confidence long non-coding RNAs from de novo assembly, transcripts were filtered based on:
- Removal of mono-exonic genes and genes < 200nt.
- Exclusion of Ensembl coding genes and pseudogenes.
- Coding potential assessment using `CPAT (v1.2.4)`.
- Specific biotype filtering (lincRNA, antisense, processed_transcript, etc.).

### 3. Quantification & Downstream Analysis
- **Quantification:** `HTSeq (v0.11.1)` to generate raw count matrices.
- **Normalisation:** `limma-voom` transformation via the `limma` R-package (filtering out genes with < 1 count in 1/3 of samples).
- **Tumor Purity Estimation:** Computed via `ESTIMATE` immune and stromal signatures and cross-referenced with disease-centric ESTIMATE results.
- **Statistical Analysis:** 
  - Pearson’s correlation coefficient for expression analyses.
  - Kaplan-Meier Survival Analysis using `survminer` and `ggplot2` R-packages.

---

## 📂 Repository Structure
*(Tip: Update this section based on the exact files you upload!)*
- `survival_analysis.R`: R scripts for Kaplan-Meier plotting and determining best cut-offs.
- `deg_and_correlation.R`: R scripts for expression normalisation (limma-voom) and Pearson correlation.
- `lncRNA_filtering.py` / `commands.sh`: Scripts or command logs used for filtering and upstream processing.

*Note: Raw sequencing BAM/FASTQ files and controlled-access TCGA datasets are not included in this repository due to size and privacy constraints.*

