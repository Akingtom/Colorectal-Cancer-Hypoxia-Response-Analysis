# DESeq2 Differential Gene Expression Analysis

This repository contains a complete differential gene expression analysis pipeline using DESeq2 in R.

## Overview

This analysis identifies differentially expressed genes between treated and untreated samples while controlling for sequencing batch effects. The pipeline includes data preprocessing, statistical testing, and comprehensive visualization.

## Files Structure

```
├── README.md
├── deseq2_analysis.R          # Main analysis script
├── data/
│   ├── Count_matrix.csv       # Raw count matrix
│   └── design.csv            # Sample metadata
├── results/
│   ├── de_results.all.csv     # All DESeq2 results
│   ├── de_results.filtered.csv # Filtered significant results
│   └── normalized_counts.csv  # DESeq2 normalized counts
└── plots/
    ├── dispersion_plot.png
    ├── pca_plot.png
    ├── sample_heatmap.png
    ├── top_genes_heatmap.png
    ├── zscore_heatmap.png
    ├── ma_plot.png
    └── volcano_plot.png
```

## Requirements

### R packages
```r
install.packages(c("dplyr", "ggplot2", "ggrepel", "RColorBrewer", "pheatmap"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
```

## Input Data

### Count Matrix (`Count_matrix.csv`)
- Rows: genes (with gene IDs as row names)
- Columns: samples
- Values: raw read counts

### Design Matrix (`design.csv`)
- Rows: samples (matching count matrix columns)
- Columns: `Treatment` (treated/untreated), `Sequencing` (batch information)

## Analysis Pipeline

### 1. Data Loading and Preprocessing
- Load count matrix and sample metadata
- Set factor levels for experimental design
- Create DESeq2 dataset object with design formula: `~ Sequencing + Treatment`

### 2. Quality Control and Filtering
- Filter low-count genes (>10 reads in minimum group size)
- Set reference level for treatment comparison

### 3. Differential Expression Analysis
- Perform DESeq2 statistical testing
- Extract results with significance testing
- Apply multiple testing correction (Benjamini-Hochberg)

### 4. Results Processing
- Filter significant genes (padj < 0.05, |log2FC| > 1)
- Export all results and filtered results
- Generate normalized count matrix

### 5. Visualization

#### Quality Assessment
- **Dispersion Plot**: Shows gene-wise dispersion estimates
- **PCA Plot**: Sample clustering by treatment and sequencing batch
- **Sample Distance Heatmap**: Hierarchical clustering of samples

#### Results Visualization
- **Top Genes Heatmap**: Expression patterns of most significant genes
- **Z-score Heatmap**: Standardized expression of top genes
- **MA Plot**: Log2 fold change vs. mean expression
- **Volcano Plot**: Statistical significance vs. biological significance

## Key Results

The analysis identifies genes with:
- Adjusted p-value < 0.05
- Absolute log2 fold change > 1
- Controlled for sequencing batch effects

## Usage

1. **Clone the repository:**
```bash
git clone https://github.com/yourusername/deseq2-analysis.git
cd deseq2-analysis
```

2. **Prepare your data:**
   - Place `Count_matrix.csv` and `design.csv` in the `data/` directory
   - Ensure sample names match between files

3. **Run the analysis:**
```r
source("deseq2_analysis.R")
```

4. **Check results:**
   - Statistical results in `results/` directory
   - Visualization plots in `plots/` directory

## Output Files

- `de_results.all.csv`: Complete DESeq2 results for all genes
- `de_results.filtered.csv`: Significantly differentially expressed genes only
- `normalized_counts.csv`: DESeq2 size-factor normalized counts
- Various plots for quality control and results visualization

## Statistical Model

The analysis uses a generalized linear model accounting for:
- **Treatment effect**: Primary comparison of interest
- **Sequencing batch**: Technical covariate to control batch effects

Model formula: `~ Sequencing + Treatment`

## Interpretation

- **Positive log2FoldChange**: Higher expression in treated samples
- **Negative log2FoldChange**: Higher expression in untreated samples
- **padj**: Benjamini-Hochberg adjusted p-values for multiple testing

## Citation

If you use this analysis pipeline, please cite:
- Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 550 (2014).

## Contact

For questions about this analysis, please open an issue or contact [your-email@domain.com].
