# RNA-Seq Differential Gene Expression Analysis
# Tool: DESeq2 in R

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# 1. Load Data
counts <- read.table("../counts/gene_counts.txt", header=TRUE, row.names=1, sep="\t")
col_data <- data.frame(
    row.names = colnames(counts),
    condition = factor(c("Healthy", "Healthy", "Healthy", "Disease", "Disease", "Disease"))
)

# 2. Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = col_data,
                              design = ~ condition)

# 3. Run DESeq
dds <- DESeq(dds)
res <- results(dds)

# 4. Visualization - PCA Plot
vsd <- vst(dds, blind=FALSE)
pca_data <- plotPCA(vsd, intgroup="condition")
pdf("PCA_Plot.pdf")
print(ggplot(pca_data, aes(PC1, PC2, color=condition)) + 
      geom_point(size=3) + 
      ggtitle("PCA Plot: Healthy vs Disease") +
      theme_minimal())
dev.off()

# 5. Visualization - Volcano Plot
res_df <- as.data.frame(res)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Yes", "No")

pdf("Volcano_Plot.pdf")
print(ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
      geom_point(alpha=0.5) +
      scale_color_manual(values=c("No"="grey", "Yes"="red")) +
      geom_vline(xintercept=c(-1, 1), col="blue") +
      geom_hline(yintercept=-log10(0.05), col="blue") +
      labs(title="Volcano Plot", x="Log2 Fold Change", y="-log10(adj P-value)") +
      theme_minimal())
dev.off()

# 6. Visualization - Heatmap
top_genes <- head(orderBy(res, order=FALSE, decrease=TRUE), 20)
mat <- assay(vsd)[rownames(top_genes), ]
mat <- mat - rowMeans(mat) # Center the data

pdf("Heatmap.pdf")
pheatmap(mat, 
         annotation_col=col_data, 
         main="Top 20 Differentially Expressed Genes",
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

# 7. Save Results
write.csv(res, "DGE_Results.csv")
print("Analysis Complete. Plots and results saved.")
