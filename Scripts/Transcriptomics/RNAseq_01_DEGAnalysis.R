###############################################################################
# Transcriptomic Analysis - Differential expression Analysis          
#                                                                     
# Description: Differential gene expression analysis of Arabidopsis thaliana   
# seedlings, comparing the 2cpab mutant against wild type. Includes quality 
# control (correlation between replicates), PCA, hierarchical clustering and 
# two independent differential expression approaches: DESeq2 and limma-voom.                                  
#
# Organism: Arabidopsis thaliana                                      
# Replicates: 3 WT   (WT_21-1,    WT_21-2,    WT_21-3)                
#             3 2cpab (2cpab_21-4, 2cpab_21-5, 2cpab_21-6)             
# Input: WT_21-1.tsv, WT_21-2.tsv, WT_21-3.tsv                        
#        2cpab_21-4.tsv, 2cpab_21-5.tsv, 2cpab_21-6.tsv               
#        gene_count_matrix.csv                                        
# Output: 
#   tables/activated_genes_deseq2.txt   : list of activated DEGs (DESeq2)
#   tables/repressed_genes_deseq2.txt   : list of repressed DEGs (DESeq2)
#   tables/activated_genes_voom.txt     : list of activated DEGs (limma-voom)
#   tables/repressed_genes_voom.txt     : list of repressed DEGs (limma-voom)
#   images/PCA_RNAseq.png               : PCA plot of normalized expression
#   images/Volcano_DESeq2.png           : volcano plot (DESeq2)
#   images/Venn_activated.png           : DESeq2 vs limma-voom overlap (up)
#   images/Venn_repressed.png           : DESeq2 vs limma-voom overlap (down)
#   images/Heatmap_top_DEGs_RNAseq.jpg  : heatmap of top 35 DEGs
#
# Author: Nardy Celeste Ayala Pérez                                   
# Date: 2026                                                          
###############################################################################        

#-------------------------------------------------------------------------------
# PACKAGE LOADING
#-------------------------------------------------------------------------------
library(edgeR)    # DGEList, voom pipeline
library(limma)    # Linear models, eBayes, topTable     
library(DESeq2)   # Primary differential expression analysis
library(ggplot2)  # Plotting
library(pheatmap) # Heatmaps
library(VennDiagram) # Venn diagrams
library(FactoMineR) # PCA, clustering
library(factoextra) # PCA/dendrogram plots
library(NormalyzerDE) # Normalization

# Create output directories if they do not already exist
dir.create("tables", showWarnings = FALSE)
dir.create("images", showWarnings = FALSE)

# Color palette used in the volcano plot to indicate the direction of
# differential expression

mycolors        <- c("dodgerblue", "firebrick2", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

## -----------------------------------------------------------------------------
## READ DATA
## -----------------------------------------------------------------------------

# StringTie output tables (one per sample), containing TPM values per gene.

WT_1   <- read.table(file = "WT_21-1.tsv",    header = T, sep = "\t")
WT_2   <- read.table(file = "WT_21-2.tsv",    header = T, sep = "\t")
WT_3   <- read.table(file = "WT_21-3.tsv",    header = T, sep = "\t")
cpab_1 <- read.table(file = "2cpab_21-4.tsv", header = T, sep = "\t")
cpab_2 <- read.table(file = "2cpab_21-5.tsv", header = T, sep = "\t")
cpab_3 <- read.table(file = "2cpab_21-6.tsv", header = T, sep = "\t")

# Extract TPM values as named vectors (gene ID -> TPM)

gene_expression_WT_1   <- setNames(WT_1$TPM,   WT_1$Gene.ID)
gene_expression_WT_2   <- setNames(WT_2$TPM,   WT_2$Gene.ID)
gene_expression_WT_3   <- setNames(WT_3$TPM,   WT_3$Gene.ID)
gene_expression_cpab_1 <- setNames(cpab_1$TPM, cpab_1$Gene.ID)
gene_expression_cpab_2 <- setNames(cpab_2$TPM, cpab_2$Gene.ID)
gene_expression_cpab_3 <- setNames(cpab_3$TPM, cpab_3$Gene.ID)

# Use the gene order of the first sample as the reference for all others
gene.ids.tpm <- names(gene_expression_WT_1)

# Combine all six samples into a single expression matrix (genes x samples)
gene.expression <- matrix(
  data = c(gene_expression_WT_1[gene.ids.tpm],
           gene_expression_WT_2[gene.ids.tpm],
           gene_expression_WT_3[gene.ids.tpm],
           gene_expression_cpab_1[gene.ids.tpm],
           gene_expression_cpab_2[gene.ids.tpm],
           gene_expression_cpab_3[gene.ids.tpm]),
  nrow = length(gene.ids.tpm),
  ncol = 6
)
rownames(gene.expression) <- gene.ids.tpm
colnames(gene.expression) <- c("WT_1","WT_2","WT_3","cpab_1","cpab_2","cpab_3")

# Experimental design table, used throughout the script for group assignment
design.tpm <- data.frame(
  sample = colnames(gene.expression),
  group  = c(rep("WT", 3), rep("cpab", 3))
)
design.tpm

dim(gene.expression)
head(gene.expression)

## -----------------------------------------------------------------------------
## NORMALIZATION WITH NormalyzerDE
## -----------------------------------------------------------------------------
## Several normalization methods are evaluated automatically. Quantile
## normalization was selected as the best-performing method after visual
## inspection of the diagnostic PDF report generated by this function

# An offset of +1 is added before normalization to avoid zero values, which
# would cause errors during log transformation and division-based normalization

gene.expression.1 <- gene.expression + 1
write.table(x     = gene.expression.1,
            file  = "at_gene_expression.tsv",
            quote = FALSE, row.names = TRUE, col.names = NA, sep = "\t")

write.table(x         = design.tpm,
            file      = "normalyzer_design.tsv",
            quote     = FALSE, row.names = FALSE, sep = "\t")

normalyzer(jobName    = "rna_seq_WT_pcab",
           designPath = "normalyzer_design.tsv",
           dataPath   = "at_gene_expression.tsv",
           outputDir  = ".")

# Read back the Quantile-normalized matrix generated by NormalyzerDE
normalized.gene.expression <- read.table(
  file = "rna_seq_WT_pcab/Quantile-normalized.txt", header = TRUE)
rownames(normalized.gene.expression) <- gene.ids.tpm

# Rename columns from "pcab" to "cpab" for naming consistency 
colnames(normalized.gene.expression) <- gsub("pcab_", "cpab_", colnames(normalized.gene.expression))
head(normalized.gene.expression)

## -----------------------------------------------------------------------------
## CORRELATION BETWEEN BIOLOGICAL REPLICATES
## -----------------------------------------------------------------------------

## Pairwise scatter plots with Pearson correlation coefficients are used
## to assess reproducibility within each genotype group

png("images/Replicate_Correlation_RNAseq.png", width = 1200, 
    height = 800, res = 150)

par(mfrow = c(2, 3))

# WT replicate 1 vs replicate 2
plot(normalized.gene.expression[,"WT_1"],
     normalized.gene.expression[,"WT_2"],
     pch = 19, col = "grey", xlab = "WT_1", ylab = "WT_2",
     cex = 0.5, main = "Similarity WT 1-2")
text(x = 4, y = 12, labels = paste0("cor = ",
                      round(100 * cor(normalized.gene.expression[,"WT_1"],
                      normalized.gene.expression[,"WT_2"],
                      use = "complete.obs"), 2), "%"))

# WT replicate 1 vs replicate 3
plot(normalized.gene.expression[,"WT_1"],
     normalized.gene.expression[,"WT_3"],
     pch = 19, col = "grey", xlab = "WT_1", ylab = "WT_3",
     cex = 0.5, main = "Similarity WT 1-3")
text(x = 4, y = 12, labels = paste0("cor = ",
                      round(100 * cor(normalized.gene.expression[,"WT_1"],
                      normalized.gene.expression[,"WT_3"],
                      use = "complete.obs"), 2), "%"))

# WT replicate 2 vs replicate 3
plot(normalized.gene.expression[,"WT_2"],
     normalized.gene.expression[,"WT_3"],
     pch = 19, col = "grey", xlab = "WT_2", ylab = "WT_3",
     cex = 0.5, main = "Similarity WT 2-3")
text(x = 4, y = 12, labels = paste0("cor = ",
                     round(100 * cor(normalized.gene.expression[,"WT_2"],
                     normalized.gene.expression[,"WT_3"],
                     use = "complete.obs"), 2), "%"))

# 2cpab replicate 1 vs replicate 2
plot(normalized.gene.expression[,"cpab_1"],
     normalized.gene.expression[,"cpab_2"],
     pch = 19, col = "grey", xlab = "2cpab_1", ylab = "2cpab_2",
     cex = 0.5, main = "Similarity 2cpab 1-2")
text(x = 4, y = 12, labels = paste0("cor = ",
                     round(100 * cor(normalized.gene.expression[,"cpab_1"],
                     normalized.gene.expression[,"cpab_2"],
                     use = "complete.obs"), 2), "%"))

# 2cpab replicate 1 vs replicate 3
plot(normalized.gene.expression[,"cpab_1"],
     normalized.gene.expression[,"cpab_3"],
     pch = 19, col = "grey", xlab = "2cpab_1", ylab = "2cpab_3",
     cex = 0.5, main = "Similarity 2cpab 1-3")
text(x = 4, y = 12, labels = paste0("cor = ",
                      round(100 * cor(normalized.gene.expression[,"cpab_1"],
                      normalized.gene.expression[,"cpab_3"],
                      use = "complete.obs"), 2), "%"))

# 2cpab replicate 2 vs replicate 3
plot(normalized.gene.expression[,"cpab_2"],
     normalized.gene.expression[,"cpab_3"],
     pch = 19, col = "grey", xlab = "2cpab_2", ylab = "2cpab_3",
     cex = 0.5, main = "Similarity 2cpab 2-3")
text(x = 4, y = 12, labels = paste0("cor = ",
                      round(100 * cor(normalized.gene.expression[,"cpab_2"],
                      normalized.gene.expression[,"cpab_3"],
                      use = "complete.obs"), 2), "%"))

par(mfrow = c(1, 1))

## -----------------------------------------------------------------------------
## PRINCIPAL COMPONENT ANALYSIS (PCA) AND HIERARCHICAL CLUSTERING
## -----------------------------------------------------------------------------

## PCA on quantile-normalized expression to assess overall sample separation by 
## genotype. HCPC is used as a complementary clustering approach to confirm the 
## grouping observed in the PCA.

# Transpose expression matrix so that samples are rows and genes are columns,
# as required by FactoMineR's PCA() function

pca.gene.expression <- data.frame(colnames(normalized.gene.expression), 
                                  t(normalized.gene.expression))
colnames(pca.gene.expression)[1] <- "Sample"

# Remove any columns (genes) containing NA values, which would otherwise
# cause PCA() to fail

pca.gene.expression.clean <- pca.gene.expression[,
                             colSums(is.na(pca.gene.expression)) == 0]

res.pca <- PCA(pca.gene.expression.clean, graph = FALSE, scale.unit = TRUE,
               quali.sup  = 1)

# scale.unit = TRUE standardizes each gene to unit variance before PCA,
# preventing highly expressed genes from dominating the analysis.
# quali.sup = 1 excludes the "Sample" column from the PCA calculation,
# treating it only as a qualitative label

pca.plot <- fviz_pca_ind(res.pca, col.ind      = design.tpm$group,
                         pointsize    = 2,
                         pointshape   = 21,
                         fill         = "black",
                         repel        = TRUE,
                         addEllipses  = TRUE,
                         ellipse.type = "confidence",
                         legend.title = "Conditions",
                         title        = "RNA-seq Gene Expression")
ggsave("images/PCA_RNAseq.png", pca.plot, width = 10, height = 8, dpi = 600)
pca.plot

# Hierarchical Clustering

res.hcpc <- HCPC(res.pca, graph = FALSE, nb.clust = 2)

dend.plot <- fviz_dend(res.hcpc,
                       k                   = 2,
                       cex                 = 0.75,
                       palette             = "jco",
                       rect                = TRUE,
                       rect_fill           = TRUE,
                       rect_border         = "jco",
                       type                = "rectangle",
                       labels_track_height = 900)

ggsave("images/HCPC_dendrogram_RNAseq.png", dend.plot, width = 8, 
       height = 6, dpi = 300)
dend.plot

## -----------------------------------------------------------------------------
## READ RAW COUNT MATRIX FOR DIFFERENTIAL EXPRESSION
## -----------------------------------------------------------------------------

## DESeq2 and limma-voom require raw integer read counts (not TPM), as their
## internal normalization methods model count data directly

gene.count.matrix.raw <- read.table( file   = "gene_count_matrix.csv", 
                                     header = TRUE, sep    = ",")
head(gene.count.matrix.raw)

## Gene IDs in the original file follow the format AT1G01010|NAC001;
## only the TAIR locus identifier (before "|") is retained.

# Keep only the TAIR locus ID before "|"
gene.count.matrix.raw$gene_id <- sub("\\|.*", "", gene.count.matrix.raw$gene_id)
gene.ids <- gene.count.matrix.raw$gene_id
head(gene.ids)

# Remove the gene_id column (already saved separately as gene.ids) and use
# it as row names, producing a numeric matrix required by DESeq2

gene.count.matrix <- gene.count.matrix.raw[, -1]
rownames(gene.count.matrix) <- gene.ids

# Original column order in file is cpab_4/5/6, WT_1/2/3, rename and reorder
colnames(gene.count.matrix) <- c("cpab_1","cpab_2","cpab_3",
                                 "WT_1","WT_2","WT_3")

gene.count.matrix <- gene.count.matrix[, c("WT_1","WT_2","WT_3",
                                           "cpab_1","cpab_2","cpab_3")]
head(gene.count.matrix)
dim(gene.count.matrix)

## -----------------------------------------------------------------------------
## SAMPLE INFORMATION
## Shared experimental design table used by both DESeq2 and limma-voom.
## -----------------------------------------------------------------------------

condition <- c(rep("WT", 3), rep("cpab", 3))
type      <- rep("single-end", 6)
coldata   <- data.frame(condition, type)
rownames(coldata) <- colnames(gene.count.matrix)
coldata

## -----------------------------------------------------------------------------
## DESeq2 — PRIMARY DIFFERENTIAL EXPRESSION ANALYSIS
## -----------------------------------------------------------------------------

# Build the DESeq2 dataset object from the raw count matrix, sample metadata,
# and the experimental design formula (condition: WT vs 2cpab)

dds <- DESeqDataSetFromMatrix(countData = gene.count.matrix, 
                              colData   = coldata, design    = ~ condition)
dds

dds <- DESeq(dds)

# Extract results for the contrast 2cpab vs WT
# 2cpab = numerator (mutant), WT = denominator (control)
# Positive log2FoldChange = higher expression in 2cpab, negative = higher in WT

res.deseq2    <- results(dds, contrast = c("condition", "cpab", "WT"))
res.deseq2.df <- as.data.frame(res.deseq2)
head(res.deseq2.df)

# Extract vectors
gene.ids.deseq2 <- rownames(res.deseq2.df)
log.fc.deseq2   <- res.deseq2.df$log2FoldChange
adj.p.deseq2    <- res.deseq2.df$padj
names(log.fc.deseq2) <- gene.ids.deseq2
names(adj.p.deseq2)  <- gene.ids.deseq2

# Differentially expressed genes (DEGs) are defined using a fold-change
# threshold of 1.5 (|log2FC| > log2(1.5)) combined with statistical
# significance (BH-adjusted p-value < 0.05)

activated.genes.deseq2 <- gene.ids.deseq2[
  !is.na(log.fc.deseq2) & !is.na(adj.p.deseq2) &
    log.fc.deseq2 > log2(1.5) & adj.p.deseq2 < 0.05]

repressed.genes.deseq2 <- gene.ids.deseq2[
  !is.na(log.fc.deseq2) & !is.na(adj.p.deseq2) &
    log.fc.deseq2 < -log2(1.5) & adj.p.deseq2 < 0.05]

cat("Activated genes (DESeq2):", length(activated.genes.deseq2))
cat("Repressed genes (DESeq2):", length(repressed.genes.deseq2))

# Save the final DEG lists (one gene ID per line) for downstream functional
# enrichment analyses (GO, KEGG)

write.table(activated.genes.deseq2,
            file = "tables/activated_genes_deseq2.txt",
            sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(repressed.genes.deseq2,
            file = "tables/repressed_genes_deseq2.txt",
            sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(res.deseq2.df,
            file = "tables/full_results_deseq2.tsv",
            sep = "\t", quote = FALSE)

## -----------------------------------------------------------------------------
## VOLCANO PLOT — DESeq2 
## -----------------------------------------------------------------------------

# Classify each gene as UP, DOWN or NO based on the same thresholds used
# to define the DEG lists above (|log2FC| > log2(1.5), adjusted p < 0.05)

res.deseq2.df$diffexpressed <- "NO"
res.deseq2.df$diffexpressed[
  !is.na(res.deseq2.df$log2FoldChange) & !is.na(res.deseq2.df$padj) &
    res.deseq2.df$log2FoldChange > log2(1.5) &
    res.deseq2.df$padj < 0.05] <- "UP"
res.deseq2.df$diffexpressed[
  !is.na(res.deseq2.df$log2FoldChange) & !is.na(res.deseq2.df$padj) &
    res.deseq2.df$log2FoldChange < -log2(1.5) &
    res.deseq2.df$padj < 0.05] <- "DOWN"

# Volcano plot: log2 fold change vs statistical significance, with dashed
# red lines marking the fold-change and significance thresholds

volcano.deseq2 <- ggplot(data = res.deseq2.df,
                         aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point(size = 0.8) +
  theme_minimal() +
  scale_color_manual(values = mycolors) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), col = "red", lty = 2) +
  geom_hline(yintercept = -log10(0.05), col = "red", lty = 2) +
  ggtitle("Volcano Plot: 2cpab vs WT (DESeq2)") +
  xlab("log2(Fold Change)") + ylab("-log10(adj. p-value)")

ggsave("images/Volcano_DESeq2.png", volcano.deseq2, width = 8, height = 6, dpi = 300)
volcano.deseq2

## -----------------------------------------------------------------------------
## limma-voom — DIFFERENTIAL EXPRESSION (COMPARISON)
## Used as a sensitivity comparison against the primary DESeq2 analysis.
## -----------------------------------------------------------------------------

# Build a DGEList object from the raw count matrix and apply TMM
# (Trimmed Mean of M-values) normalization to account for differences
# in library size and composition between samples

d0      <- DGEList(gene.count.matrix)
d0.norm <- calcNormFactors(d0, method = "TMM")

# Filter low-expressed genes: keep genes with CPM > 2 in at least one sample
# reducing noise from genes with unreliable low counts

cutoff <- 2
drop   <- which(apply(cpm(d0.norm), 1, max) < cutoff)
d      <- d0.norm[-drop, ]
cat("Genes before filtering:", nrow(d0.norm))
cat("Genes after  filtering:", nrow(d))

group <- as.factor(condition)

# Model matrix with one coefficient per group, no intercept
mm <- model.matrix(~ 0 + group)
colnames(mm) <- levels(group)   # "WT" and "cpab"

# voom transforms raw counts to log2-CPM and estimates the mean-variance
# relationship to assign precision weights to each observation, which are
# then used in the linear model fit

y <- voom(d, mm, plot = FALSE)

## -----------------------------------------------------------------------------
## LINEAR MODEL FIT,  CONTRASTS AND DEGs IDENTIFICATION(limma-voom)
## -----------------------------------------------------------------------------

# Fit a linear model to the voom-transformed expression values, then define
# and apply the contrast of interest (2cpab vs WT)

fit                 <- lmFit(y, mm)
contrast.matrix     <- makeContrasts(cpab - WT, levels = colnames(coef(fit)))
contrast.linear.fit <- contrasts.fit(fit, contrast.matrix)

# eBayes applies empirical Bayes shrinkage to the variance estimates,
# improving statistical power for genes with low replicate numbers

contrast.results    <- eBayes(contrast.linear.fit)

cpab.vs.wt.voom <- topTable(contrast.results,
                            number  = nrow(contrast.results),
                            coef    = 1,
                            sort.by = "logFC")
head(cpab.vs.wt.voom)

# Extract log2 fold change and adjusted p-value as named vectors

fold.change.voom <- cpab.vs.wt.voom$logFC
adj.pval.voom    <- cpab.vs.wt.voom$adj.P.Val
gene.ids.voom    <- rownames(cpab.vs.wt.voom)
names(fold.change.voom) <- gene.ids.voom
names(adj.pval.voom)    <- gene.ids.voom

# Same DEG thresholds as DESeq2 (|log2FC| > log2(1.5), adjusted p < 0.05),
# applied here for direct comparability between the two methods

activated.genes.voom <- gene.ids.voom[
  fold.change.voom > log2(1.5) & adj.pval.voom < 0.05]
repressed.genes.voom <- gene.ids.voom[
  fold.change.voom < -log2(1.5) & adj.pval.voom < 0.05]

cat("Activated genes (voom):", length(activated.genes.voom))
cat("Repressed genes (voom):", length(repressed.genes.voom))

# Save DEG lists for the Venn diagram comparison against DESeq2 results

write.table(activated.genes.voom,
            file = "tables/activated_genes_voom.txt",
            sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(repressed.genes.voom,
            file = "tables/repressed_genes_voom.txt",
            sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

## -----------------------------------------------------------------------------
## VENN DIAGRAMS — DESeq2 vs limma-voom DEG OVERLAP
## -----------------------------------------------------------------------------

## In this dataset, DESeq2 identified more DEGs than limma-voom, with voom
## DEGs largely contained within the DESeq2 DEG set, supporting the use of
## DESeq2 for downstream functional analyses. This reflects the specific
## results observed here rather than a general property of either method
# Overlap between activated DEGs identified by each method

png("images/Venn_activated.png", width = 800, height = 600, res = 150)
grid.newpage()
draw.pairwise.venn(
  area1      = length(activated.genes.deseq2),
  area2      = length(activated.genes.voom),
  cross.area = length(intersect(activated.genes.deseq2, activated.genes.voom)),
  lwd        = 3,
  category   = c("DESeq2", "limma-voom"),
  euler.d    = TRUE,
  col        = c("dodgerblue", "firebrick2"),
  fill       = c("dodgerblue", "firebrick2"),
  alpha      = 0.3,
  cex        = 1.5,
  cat.cex    = 1.5,
  main       = "Activated genes: DESeq2 vs limma-voom"
)

dev.off()

# Overlap between repressed DEGs identified by each method

png("images/Venn_repressed.png", width = 800, height = 600, res = 150)
grid.newpage()

draw.pairwise.venn(
  area1      = length(repressed.genes.deseq2),
  area2      = length(repressed.genes.voom),
  cross.area = length(intersect(repressed.genes.deseq2, repressed.genes.voom)),
  lwd        = 3,
  category   = c("DESeq2", "limma-voom"),
  euler.d    = TRUE,
  col        = c("dodgerblue", "firebrick2"),
  fill       = c("dodgerblue", "firebrick2"),
  alpha      = 0.3,
  cex        = 1.5,
  cat.cex    = 1.5,
  main       = "Repressed genes: DESeq2 vs limma-voom"
)
dev.off()
## -----------------------------------------------------------------------------
## HEATMAP: TOP 35 DEGs BY ADJUSTED P-VALUE
## -----------------------------------------------------------------------------

# Combine activated and repressed DEGs, then select the top 35 by adjusted
# p-value

degs.deseq2         <- unique(c(activated.genes.deseq2, repressed.genes.deseq2))
top35.df    <- res.deseq2.df[degs.deseq2, ]
top35.df    <- top35.df[order(top35.df$padj), ]
top35       <- head(rownames(top35.df), 35)

# Use quantile-normalized expression values for visualization

heatmap_data <- as.matrix(normalized.gene.expression[top35, ])
heatmap_data <- heatmap_data[complete.cases(heatmap_data), ]
heatmap_data <- heatmap_data[apply(heatmap_data, 1, var, na.rm = TRUE) > 0, ]

# Display genotype labels as "WT" / "2cpab" in column names and annotation
colnames(heatmap_data) <- gsub("cpab_", "2cpab_", colnames(heatmap_data))

annotation_col_top <- data.frame(
  Condition = factor(c(rep("WT", 3), rep("2cpab", 3)), levels = c("WT", "2cpab"))
)
rownames(annotation_col_top) <- colnames(heatmap_data)
ann_colors_top <- list(Condition = c(WT = "#36454F", `2cpab` = "#2E9B9B"))

# Row-wise scaling (Z-score) highlights relative expression patterns across
# samples for each gene, independent of absolute expression level

pheatmap(heatmap_data, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE,
         show_rownames = TRUE, fontsize_row = 7, fontsize_col = 13,
         annotation_col = annotation_col_top, annotation_colors = ann_colors_top,
         color = colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(50),
         main = "Top DEGs \u2014 2cpab vs WT (RNA-seq)",
         border_color = "grey60",
         filename = "images/Heatmap_top_DEGs_RNAseq.jpg",
         width = 8, height = 10)

# End of script
