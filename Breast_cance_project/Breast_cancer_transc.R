
setwd("C:/Users/JESUTOFUNMI/OneDrive/Desktop/Chuggy/Breast_cance_project")

#load libraries
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

#read count data
R.utils::gunzip("GSE272101_gene_count_ADR.txt.gz")
count_data <- read.delim("GSE272101_gene_count_ADR.txt", row.names = 1)
head(count_data)
dim(count_data)


#read sample_info
sample_info <- read.csv("Design.csv")
sample_info

#Separate counts and annotations
count_matrix <- as.matrix(count_data[, 1:6])  # First 6 columns are counts
gene_annot <- count_data[, 7:ncol(count_data)]  # Remaining columns are annotations

#set factor levels
sample_info$Group <- factor(sample_info$Group)
levels(sample_info$Group)

sample_info$Group


# Trim whitespace from all factor levels
levels(sample_info$Group) <- trimws(levels(sample_info$Group))
levels(sample_info$Group)

#create Deseq objecr
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,  # Must be numeric
  colData = sample_info,     # Must match count columns
  design = ~ Group           # Compare groups
)

# Add gene annotations
mcols(dds) <- cbind(mcols(dds), gene_annot[rownames(dds),])

dds


levels(dds$Group) <- trimws(levels(dds$Group))
levels(dds$Group)

#set the refrerence for the group level
dds$Group <- relevel(dds$Group, ref="Control")
dds$Group

#filter genes
#keep gene with at least N counts >=10, where N= size of smallest group
keep <- rowSums(counts(dds) >=10) >= min(table(sample_info$Group)) 
dds <- dds[keep,]
dds

#perform statistical test to identify differentially expressed genes
dds <- DESeq(dds,test = "Wald",sfType = "poscounts")
deseq_results <- results(dds)
deseq_results

deseq_results <- as.data.frame(deseq_results)
class(deseq_results)

head(deseq_results)

names(deseq_results)

# Creating a new column 'GeneID' from the row names for easier export and use.
deseq_results$GeneID <- rownames(deseq_results)

#Merge with gene annotations
deseq_results <- merge(deseq_results, gene_annot, by.x = "GeneID", by.y = "row.names")

#Save results
write.table(deseq_results, "deseq_results.tsv", sep = "\t", row.names = FALSE)


#Identify DEGs 
deg <- subset(deseq_results, padj < 0.05 & abs(log2FoldChange) >= 1)
deg <- deg[order(deg$padj),]
write.table(deg, "deseq_deg.tsv", sep = "\t", row.names = FALSE)


# MA Plot
  plotMA(res, ylim = c(-2, 2))
  
#plot dispersion estimates
plotDispEsts(dds,main="GSE272101 Dispersion Estimates")
plotDispEsts

#histogram plot of p-values
hist(deseq_results$padj, 
     breaks = seq(0, 1, length = 21), 
     col = "grey", 
     border = "white",
     xlab = "",
     ylab = "",
     ylim = c(0,8000),
     main = "GSE272101 Frequencies of padjvalue")

#plot volcano

#Get the group levels
groups <- levels(dds$Group)

#Create the title
volcano_title <- paste(groups[1], "vs", groups[2])

#set colours
old.pal <- palette(c("#00BFFF", "#FF3030"))

#Set margin size
par(mar=c(4,4,2,1), cex.main=1.5)

#set title
title = paste(groups[1], "vs", groups[2])

#plot values
plot(deseq_results$log2FoldChange, -log10(deseq_results$padj),
     main = title,
     xlab = "log2FC",
     ylab = "-log10(padj)",
     pch = 20,
     cex = 0.5)


with( subset(deseq_results, padj <0.05 & abs(log2FoldChange) >=1),
      points(log2FoldChange, -log10(padj), pch=20, col=(sign(log2FoldChange) +3)/2, cex=1))

legend("bottomleft", title=paste("padj<", 0.05, sep=""),
       legend=c("down", "up"), pch=20, col=1:2)


# Heatmap of top genes
vsd <- vst(dds, blind = FALSE)

#plot PCA

plotPCA(vsd, intgroup=c("Group"))

#heatmap of log transformed normalized counts.
#top 10 genes

normalized_counts <- counts(dds, normalized=T)
head(normalized_counts)

transformed_counts <- log2(normalized_counts+1)
head(transformed_counts)

top_hits <- row.names(deg[1:10, ])
top_hits

colnames(deg)
top_hits <- deg$GeneID[1:10]
top_hits <- transformed_counts[top_hits, ]
head(top_hits)

pheatmap(top_hits,cluster_rows = FALSE,cluster_cols = FALSE)

# 8. Functional Enrichment 
# Convert gene IDs to symbols if needed
gene_list <- deg$gene_name  # Using gene_name column from your data
gene_list

library(gprofiler2)

# GO enrichment
go_results <- gost(
  query = gene_list,
  organism = "hsapiens",
  sources = c("GO", "KEGG"),
  significant = TRUE
)

# View the full results table
(go_results$result)

keggResults <- gost(query = gene_list,
                    organism = 'hsapiens',
                    source = 'KEGG',
                    user_threshold = 0.05,
                    correction_method = 'g_SCS')
# View the full KEGG results table
keggResults$result

# Visualize top GO terms
gostplot(go_results, capped = FALSE, interactive = TRUE)

# Visualize KEGG pathways
gostplot(keggResults, capped = FALSE, interactive = TRUE)


