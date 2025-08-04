#load libraries
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

setwd("C:/Users/JESUTOFUNMI/OneDrive/Desktop/Chuggy/Project")

#read count data
count_table <- read.csv("raw_counts.tsv", sep = "\t", row.names=1)
head(count_table)
dim(count_table)

#read sample info.
sample_info <- read.csv("design.tsv", sep = "\t", row.names=1)
sample_info
dim(sample_info)

#set factors_level
factors <- factor(sample_info$Group)
groups <- unique(sample_info$Group)
groups

groups <- rev(groups)
groups

sample_info$Group <- factors

sample_info$Group

#create Deseq object
dds <- DESeqDataSetFromMatrix(
  countData = count_table,
  colData = sample_info,
  design = ~ Group
)


#set the refrerence for the group level
dds$Group <- relevel(dds$Group, ref="control")



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

deseq_results$Genename <- row.names(deseq_results)
head(deseq_results)

deseq_results <- subset(deseq_results,
                        select=c("Genename","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean"))

names(deseq_results)

write.table(deseq_results, file="deseq_results.tsv", row.names=F, sep="\t")




#extract results for top genes
#extract genes with padj < 0.05 and log2foldchange <=-1 or >=1

deg <- subset(deseq_results,padj<0.05 & abs(log2FoldChange) >=1)
head(deg)
dim(deg)


deg <- deg[order(deg$padj),]
head(deg)

write.table(deg,file="deseq_deg.tsv",row.names=F, sep="\t")


#Gene expression data visualization

#plot dispersion estimates
plotDispEsts(dds,main="GSE203159 Dispersion Estimates")

#histogram plot of p-values
hist(deseq_results$padj, 
     breaks = seq(0, 1, length = 21), 
     col = "grey", 
     border = "white",
     xlab = "",
     ylab = "",
     ylim = c(0,8000),
     main = "GSE203159 Frequencies of padjvalue"
     

#plot volcano
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


#variance stabilizing transformation
vsd <- vst(dds, blind = F)

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

top_hits <- transformed_counts[top_hits, ]

head(top_hits)

pheatmap(top_hits,cluster_rows = FALSE,cluster_cols = FALSE)


#plot MA
plotMA(dds, ylim=c(-2,2))

resultsNames(dds)

#remove noise
# Assuming the name is "Group_tgf.beta_vs_control"
resLFC <- lfcShrink(dds, coef = "Group_tgf.beta_vs_control", type = "apeglm")

plotMA(resLFC, ylim=c(-2,2))

library(DESeq2)
library(knitr)
library(gProfileR)

install.packages("gProfileR")
install.packages("gprofiler2")
library(gprofiler2)

#calculate enriched GO terms

# First, filter for significant genes
sig_genes <- subset(deseq_results, padj < 0.05 & abs(log2FoldChange) >= 1)

# Now, create a vector of the gene names
gene_list <- sig_genes$Genename

goResults <- gost(query = gene_list,
                  organism = 'hsapiens',
                  source = 'GO',
                  user_threshold = 0.05,
                  correction_method = 'g_SCS')
# View the full results table
(goResults$result)

#create KEGG Enrichment
# A list of significant genes is needed for the query
# Assuming you have already created a 'gene_list' variable
# containing the names of your differentially expressed genes.

keggResults <- gost(query = gene_list,
                    organism = 'hsapiens',
                    source = 'KEGG',
                    user_threshold = 0.05,
                    correction_method = 'g_SCS')
# View the full KEGG results table
keggResults$result
