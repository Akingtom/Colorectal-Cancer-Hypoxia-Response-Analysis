count_data <- read.csv("Count_matrix.csv", header = TRUE, row.names = 1)
colnames(count_data)
head(count_data)
sample_info <- read.csv("design.csv", header = TRUE, row.names = 1)
colnames(sample_info)
head(sample_info)

#set factor levels
sample_info$Treatment <- factor(sample_info$Treatment)
sample_info$Sequencing <- factor(sample_info$Sequencing)


library(DESeq2)

#create a deseq object and import the count_data and sample_info
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = sample_info,
  design = ~ Sequencing + Treatment
)

#set the reference for the treatment factor
dds$Treatment <- factor(dds$Treatment, levels=c("untreated", "treated"))

#filter genes
keep <- rowSums(counts(dds)>10) >= min(table(sample_info$Treatment)) 
dds <- dds[keep,]

# Check the factor levels of your design variables after filtering
table(colData(dds)$Treatment)
table(colData(dds)$Sequencing)


#perform statistical test to identify differentially expressed genes
dds <- DESeq(dds)
deseq_results <- results(dds)
deseq_results


#change Deseq objects to object(dataframe)
deseq_results <- as.data.frame(deseq_results)
class(deseq_results)
head(deseq_results)

#order the result table by increasing p value
deseq_results_ordered <- deseq_results[order(deseq_results$pvalue),]
head(deseq_results_ordered)


deseq_results["FBgn0003360", ]

deseq_results["FBgn0261552",]

library(dplyr)

filtered <- deseq_results %>% filter(deseq_results$padj < 0.05)

filtered <- filtered %>% filter(abs(filtered$log2FoldChange) > 1)

dim(deseq_results)
dim(filtered)

write.csv(deseq_results, 'de_results.all.csv')
write.csv(filtered, 'de_results.filtered.csv')

normalized_counts <- counts(dds,normalized=TRUE)
head(normalized_counts)

write.csv(normalized_counts, 'normalized_counts.csv')

#Visualization

#plot Dispersion Estimates
plotDispEsts(dds)

vsd <- vst(dds, blind=FALSE)

#plot PCA
plotPCA(vsd, intgroup=c('Sequencing', 'Treatment'))

#generate distance matrix
SampleDists <- dist(t(assay(vsd)))

SampleDistMatrix <- as.matrix(SampleDists)
colnames(SampleDistMatrix)

library(RColorBrewer)

#set a color scheme
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

library(pheatmap)

#generate the heatmap for distance matrix
pheatmap(SampleDistMatrix, clustering_distance_rows=SampleDists, 
             clustering_distance_cols=SampleDists,col=colors)
#top genes
top_hits <- deseq_results[order(deseq_results$padj), ][1:10, ]
top_hits <- row.names(top_hits)
top_hits

rld <- rlog(dds,blind=FALSE)


pheatmap(assay(rld)[top_hits,],cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE)
pheatmap(assay(rld)[top_hits,])


annot_info <- as.data.frame(colData(dds)[,c('Sequencing','Treatment')])
pheatmap(assay(rld)[top_hits,],cluster_rows=FALSE, show_rownames=TRUE, 
         cluster_cols=FALSE, annotation_col =annot_info)

#Heat map of Z score
cal_z_score <- function(x)((x-mean(x)) / sd(x))

zscore_all <- t(apply(normalized_counts, 1, cal_z_score))
zscore_subset <- zscore_all[top_hits,]
pheatmap(zscore_subset)


#MA plot
plotMA(dds, ylim=c(-2,2))

#remove noise
resLFC <- lfcShrink(dds,coef = "Treatment_treated_vs_untreated", type="apeglm")

plotMA(resLFC, ylim=c(-2,2))

resLFC <- as.data.frame(resLFC)

#label the genes

resLFC$diffexpressed <- "No"
resLFC$diffexpressed[resLFC$log2FoldChange>0.1 & resLFC$padj<0.05] <- "UP"
resLFC$diffexpressed[resLFC$log2FoldChange<0.1 & resLFC$padj<0.05] <- "DOWN"

resLFC$delabel<-NA


library(ggplot2)
library(ggrepel)

volcano_plot <- ggplot(data=resLFC, aes(x = log2FoldChange, y = -log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c('blue','black','red')) +
  theme(text=element_text(size=20))
print(volcano_plot)
