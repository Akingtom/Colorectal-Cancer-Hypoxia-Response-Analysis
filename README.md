# Colorectal Cancer Hypoxia Response Analysis

## Overview

This project performs comprehensive differential expression and pathway enrichment analysis on colorectal cancer RNA-seq data (GEO accession: GSE197576) to investigate the transcriptional response to hypoxia and the effects of RELB and ITPR3 knockdown under both normoxic and hypoxic conditions.

## Dataset

**GEO Accession**: [GSE197576](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197576)

The dataset contains RNA-seq gene counts from SW480 colorectal cancer cells under various conditions:
- Control (sgCTRL) - Normoxia and Hypoxia
- RELB knockdown (sgRELB) - Normoxia and Hypoxia  
- ITPR3 knockdown (sgITPR3) - Normoxia and Hypoxia

## Analysis Pipeline

### 1. Data Acquisition
- Downloads raw gene count matrix from GEO
- Processes count data for downstream analysis

### 2. Differential Expression Analysis (DESeq2)
The analysis performs multiple comparisons:
- Hypoxia vs Normoxia (sgCTRL)
- sgRELB_Norm vs sgCTRL_Norm
- sgITPR3_Norm vs sgCTRL_Norm
- sgRELB_Hyp vs sgCTRL_Hyp
- sgITPR3_Hyp vs sgCTRL_Hyp

**Significance thresholds**: 
- Adjusted p-value ≤ 0.05
- |log2FoldChange| ≥ 2

### 3. Pathway Enrichment Analysis

#### GO Term Enrichment
- Biological Process (BP) ontology
- Background: all detected genes in the experiment
- Multiple testing correction: Benjamini-Hochberg

#### Hallmark Gene Set Enrichment
- Over-representation analysis (ORA) using MSigDB Hallmark gene sets
- Gene Set Enrichment Analysis (GSEA) with ranked gene lists

### 4. Visualization
- Volcano plots for all comparisons
- Barplots and dotplots for enriched pathways
- GSEA enrichment plots for key pathways (G2M checkpoint, Hypoxia response)

## Key Findings

The hypoxia vs normoxia comparison identified **444 significant genes** with strong enrichment in:
- **Hypoxia response** (HALLMARK_HYPOXIA)
- **Epithelial-mesenchymal transition** 
- **KRAS signaling**
- **TNFα signaling via NF-κB**
- **Inflammatory response**
- **Angiogenesis** pathways

## Requirements

### R Packages
```r
# Bioconductor packages
BiocManager::install(c("GEOquery", "DESeq2", "clusterProfiler", 
                       "org.Hs.eg.db", "enrichplot"))

# CRAN packages
install.packages(c("dplyr", "readr", "here", "purrr", "tibble", 
                   "stringr", "ggplot2", "msigdbr"))
```

## Usage

1. Clone this repository
2. Open the R script in RStudio
3. Run the analysis sequentially
4. Results will be saved to your specified output directory

## Project Structure

```
colorectal/
├── GSE197576_raw_gene_counts_matrix.tsv.gz  # Downloaded count matrix
├── analysis.R                                # Main analysis script
└── results/
    ├── sgCTRL_Hyp_vs_sgCTRL_Norm.pdf        # Pathway enrichment plots
    ├── sgRELB_Norm_vs_sgCTRL_Norm.pdf
    ├── sgITPR3_Norm_vs_sgCTRL_Norm.pdf
    ├── sgRELB_Hyp_vs_sgCTRL_Hyp.pdf
    └── sgITPR3_Hyp_vs_sgCTRL_Hyp.pdf
```

## Methodology Highlights

- **Functional programming approach**: Uses `purrr::map()` and list-columns for batch processing multiple comparisons
- **Robust GSEA scoring**: Handles infinite p-values and missing gene mappings
- **Comprehensive QC**: Pre-filtering of low-count genes (>10 reads in ≥2 samples)

## Author

Jesutofunmi

## References

- Love MI, Huber W, Anders S (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." *Genome Biology*
- Yu G et al. (2024). "Thirteen years of clusterProfiler." *The Innovation*
- Korotkevich G et al. (2019). "Fast gene set enrichment analysis." *bioRxiv*

## License

This project is available for academic and research use.
