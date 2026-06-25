###############################################################################
# Proteomics Analysis - Differential Protein Abundance Analysis (SWATH-MS)
#
# Description: Differential protein abundance analysis of Arabidopsis
# thaliana seedlings, comparing the 2cpab mutant against wild-type (WT),
# based on SWATH-MS quantitative proteomics data. Includes quality control
# (correlation between replicates, PCA) and differential abundance testing
# using limma, followed by volcano plot and heatmap visualization of the
# most significant differentially abundant proteins (DEPs).
#
# Organism: Arabidopsis thaliana
# Replicates: 3 WT (wt1, wt2, wt3), 3 2cpab (ab1, ab2, ab3)
#
# Input:
#   SWATH WT-ab(Normalizadas).csv : normalized SWATH-MS protein intensities,
#                                   semicolon-separated, European decimal format
#
# Output:
#   tables/full_results_proteomics.csv   : complete limma results table
#   tables/DEPs_supplementary_table.xlsx : full DEP table with UniProt ID,
#                                          gene name, and subcellular
#                                          localization, split by Up/Down
#   images/Replicate_Correlation_Proteomics.png : correlation between replicates
#   images/PCA_proteomics.png               : PCA plot of normalized abundance
#   images/Volcano_proteomics.png           : volcano plot (limma)
#   images/Heatmap_top_DEPs.jpg             : heatmap of top DEPs
#
# Author: Nardy Celeste Ayala Pérez
# Date: 2026
###############################################################################

#-------------------------------------------------------------------------------
# PACKAGE LOADING
#-------------------------------------------------------------------------------
library(limma)            # Linear models, eBayes, topTable
library(ggplot2)          # Plotting
library(tidyverse)        # Data manipulation
library(pheatmap)         # Heatmap visualization
library(FactoMineR)       # PCA
library(factoextra)       # PCA visualization
library(clusterProfiler)  # bitr(), ID conversion
library(org.At.tair.db)   # Arabidopsis thaliana genome annotation
library(AnnotationDbi)    # select(), annotation queries
library(GO.db)            # GO term descriptions
library(writexl)          # Exporting results to Excel

dir.create("tables", showWarnings = FALSE)
dir.create("images", showWarnings = FALSE)


#-------------------------------------------------------------------------------
# LOAD SWATH-MS PROTEIN ABUNDANCE DATA
#-------------------------------------------------------------------------------

data <- read.csv("SWATH WT-ab(Normalizadas).csv", sep = ";",
                 header = TRUE, check.names = FALSE)

# Confirm that the expected sample columns are present
stopifnot(all(c("ab1", "ab2", "ab3", "wt1", "wt2", "wt3") %in% colnames(data)))


#-------------------------------------------------------------------------------
# LOG2 TRANSFORMATION
#-------------------------------------------------------------------------------

# Sample columns use the European decimal comma format and must be
# converted to numeric. 

sample.cols <- c("ab1", "ab2", "ab3", "wt1", "wt2", "wt3")

num.expr <- as.data.frame(
  lapply(data[, sample.cols], function(x) as.numeric(gsub(",", ".", x)))
)
rownames(num.expr) <- data[["Peak Name"]]

# log2(x + 1) is applied to stabilize variance, avoiding errors for proteins 
# with zero intensity.

log.expr <- log2(num.expr + 1)

#-------------------------------------------------------------------------------
# CORRELATION BETWEEN BIOLOGICAL REPLICATES
#-------------------------------------------------------------------------------

# Pairwise scatter plots with Pearson correlation coefficients are used
# to assess reproducibility within each genotype group

png("images/Replicate_Correlation_Proteomics.png", width = 1200, 
    height = 800, res = 150)

par(mfrow = c(2, 3))

# WT replicate 1 vs replicate 2
plot(log.expr[, "wt1"], log.expr[, "wt2"],
     pch = 19, col = "grey", xlab = "WT_1", ylab = "WT_2",
     cex = 0.5, main = "WT replicate 1 vs 2")
text(x = min(log.expr[, "wt1"]), y = max(log.expr[, "wt2"]),
     labels = paste0("cor = ", round(100 * cor(log.expr[, "wt1"], log.expr[, "wt2"],
                                               use = "complete.obs"), 2), "%"))

# WT replicate 1 vs replicate 3
plot(log.expr[, "wt1"], log.expr[, "wt3"],
     pch = 19, col = "grey", xlab = "WT_1", ylab = "WT_3",
     cex = 0.5, main = "WT replicate 1 vs 3")
text(x = min(log.expr[, "wt1"]), y = max(log.expr[, "wt3"]),
     labels = paste0("cor = ", round(100 * cor(log.expr[, "wt1"], log.expr[, "wt3"],
                                               use = "complete.obs"), 2), "%"))

# WT replicate 2 vs replicate 3
plot(log.expr[, "wt2"], log.expr[, "wt3"],
     pch = 19, col = "grey", xlab = "WT_2", ylab = "WT_3",
     cex = 0.5, main = "WT replicate 2 vs 3")
text(x = min(log.expr[, "wt2"]), y = max(log.expr[, "wt3"]),
     labels = paste0("cor = ", round(100 * cor(log.expr[, "wt2"], log.expr[, "wt3"],
                                               use = "complete.obs"), 2), "%"))

# 2cpab replicate 1 vs replicate 2
plot(log.expr[, "ab1"], log.expr[, "ab2"],
     pch = 19, col = "grey", xlab = "2cpab_1", ylab = "2cpab_2",
     cex = 0.5, main = "2cpab replicate 1 vs 2")
text(x = min(log.expr[, "ab1"]), y = max(log.expr[, "ab2"]),
     labels = paste0("cor = ", round(100 * cor(log.expr[, "ab1"], log.expr[, "ab2"],
                                               use = "complete.obs"), 2), "%"))

# 2cpab replicate 1 vs replicate 3
plot(log.expr[, "ab1"], log.expr[, "ab3"],
     pch = 19, col = "grey", xlab = "2cpab_1", ylab = "2cpab_3",
     cex = 0.5, main = "2cpab replicate 1 vs 3")
text(x = min(log.expr[, "ab1"]), y = max(log.expr[, "ab3"]),
     labels = paste0("cor = ", round(100 * cor(log.expr[, "ab1"], log.expr[, "ab3"],
                                               use = "complete.obs"), 2), "%"))

# 2cpab replicate 2 vs replicate 3
plot(log.expr[, "ab2"], log.expr[, "ab3"],
     pch = 19, col = "grey", xlab = "2cpab_2", ylab = "2cpab_3",
     cex = 0.5, main = "2cpab replicate 2 vs 3")
text(x = min(log.expr[, "ab2"]), y = max(log.expr[, "ab3"]),
     labels = paste0("cor = ", round(100 * cor(log.expr[, "ab2"], log.expr[, "ab3"],
                                               use = "complete.obs"), 2), "%"))

par(mfrow = c(1, 1))
dev.off()


#-------------------------------------------------------------------------------
# PRINCIPAL COMPONENT ANALYSIS (PCA)
#-------------------------------------------------------------------------------

condition.prot <- factor(c(rep("2cpab", 3), rep("WT", 3)))

res.pca.prot <- PCA(t(log.expr), graph = FALSE, scale.unit = TRUE)

pca.plot.prot <- fviz_pca_ind(res.pca.prot,
                              col.ind      = condition.prot,
                              pointsize    = 2,
                              pointshape   = 21,
                              fill         = "black",
                              repel        = TRUE,
                              addEllipses  = TRUE,
                              ellipse.type = "confidence",
                              legend.title = "Condition",
                              title        = "PCA \u2014 SWATH Proteomics")

ggsave("images/PCA_proteomics.png", pca.plot.prot, width = 8, height = 6, dpi = 300)

#-------------------------------------------------------------------------------
# EXPERIMENTAL DESIGN AND LINEAR MODEL (limma)
#-------------------------------------------------------------------------------

# Sample order: ab1, ab2, ab3, wt1, wt2, wt3 -> group 1 (2cpab), group 2 (WT)

experimental.design <- model.matrix(~ -1 + factor(c(1, 1, 1, 2, 2, 2)))
colnames(experimental.design) <- c("cpab", "wt")

linear.fit <- lmFit(log.expr, experimental.design)

# Contrast: 2cpab vs wt (positive logFC = higher abundance in 2cpab)

contrast.matrix <- makeContrasts(cpab - wt, levels = c("cpab", "wt"))

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results    <- eBayes(contrast.linear.fit)

results_table <- topTable(contrast.results, number  = nrow(log.expr),
                          coef    = 1, sort.by = "logFC")

write.csv(results_table, file = "tables/full_results_proteomics.csv")


#-------------------------------------------------------------------------------
# DIFFERENTIALLY ABUNDANT PROTEINS (DEPs)
#-------------------------------------------------------------------------------

# DEPs are defined as |log2FC| > log2(1.5) and unadjusted P.Value < 0.05.
# A stricter subset (adj.P.Val < 0.05) is identified separately below for
# the supplementary table.

FC_threshold <- log2(1.5)
pval_thresh  <- 0.05

fold.change <- results_table$logFC
prot.ids    <- rownames(results_table)

up.proteins   <- prot.ids[fold.change > FC_threshold & results_table$P.Value < pval_thresh]
down.proteins <- prot.ids[fold.change < -FC_threshold & results_table$P.Value < pval_thresh]
DEPs          <- c(up.proteins, down.proteins)

cat("Up-regulated proteins:", length(up.proteins))
cat("Down-regulated proteins:", length(down.proteins))
cat("Total DEPs:", length(DEPs))

# Proteins remaining significant after multiple-testing correction
up.proteins.sig   <- rownames(results_table[results_table$logFC > FC_threshold &
                                              results_table$adj.P.Val < 0.05, ])
down.proteins.sig <- rownames(results_table[results_table$logFC < -FC_threshold &
                                              results_table$adj.P.Val < 0.05, ])

cat("Up-regulated proteins (adj.P.Val < 0.05):", length(up.proteins.sig))
cat("Down-regulated proteins (adj.P.Val < 0.05):", length(down.proteins.sig))


#-------------------------------------------------------------------------------
# VOLCANO PLOT
#-------------------------------------------------------------------------------

volcano_data <- results_table %>% mutate(diffexpressed = case_when(
  logFC >=  FC_threshold & P.Value <= pval_thresh ~ "UP",
  logFC <= -FC_threshold & P.Value <= pval_thresh ~ "DOWN",
  TRUE ~ "NO"))

mycolors <- c("UP" = "firebrick3", "NO" = "grey50", "DOWN" = "dodgerblue3")

volcano_plot <- ggplot(volcano_data, aes(logFC, -log10(P.Value), col = diffexpressed)) +
  geom_point(size = 0.8) +
  theme_minimal() +
  scale_color_manual(values = mycolors) +
  geom_vline(xintercept = c(-FC_threshold, FC_threshold), col = "red", lty = 2) +
  geom_hline(yintercept = -log10(pval_thresh), col = "red", lty = 2) +
  ggtitle("Volcano Plot: 2cpab vs WT (SWATH proteomics)") +
  xlab("log2(Fold Change)") + ylab("-log10(P-Value)")

ggsave("images/Volcano_proteomics.png", volcano_plot, width = 8, height = 6, dpi = 300)


#-------------------------------------------------------------------------------
# HEATMAP OF TOP DEPs
#-------------------------------------------------------------------------------

if (length(DEPs) >= 2) {
  
  top_DEPs     <- head(DEPs, min(50, length(DEPs)))
  heatmap_data <- as.matrix(log.expr[top_DEPs, ])
  
  # Clean UniProt-style row names (e.g., sp|Q9C5R8|BAS1B_ARATH -> BAS1B)
  rownames(heatmap_data) <- gsub(".*\\|(.+)_ARATH", "\\1", rownames(heatmap_data))
  
  heatmap_data <- heatmap_data[, c("wt1", "wt2", "wt3", "ab1", "ab2", "ab3")]
  colnames(heatmap_data) <- c("WT_1", "WT_2", "WT_3", "2cpab_1", "2cpab_2", "2cpab_3")
  
  annotation_col <- data.frame(
    Condition = factor(c(rep("WT", 3), rep("2cpab", 3)), levels = c("WT", "2cpab"))
  )
  rownames(annotation_col) <- colnames(heatmap_data)
  
  ann_colors <- list(
    Condition = c(WT = "#36454F", `2cpab` = "#2E9B9B")
  )
  
  pheatmap(
    heatmap_data,
    scale             = "row",
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,
    show_rownames     = TRUE,
    fontsize_row      = 5,
    annotation_col    = annotation_col,
    annotation_colors = ann_colors,
    color             = colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(50),
    main              = "Top DEPs \u2014 2cpab vs WT (proteomics)",
    border_color      = "grey60",
    filename          = "images/Heatmap_top_DEPs.jpg",
    width             = 10,
    height            = 12
  )
}


#-------------------------------------------------------------------------------
# SUPPLEMENTARY TABLE: DEPs WITH GENE NAME AND SUBCELLULAR LOCALIZATION
#-------------------------------------------------------------------------------

# Extract gene symbols from UniProt-style protein IDs and convert to TAIR
dep.names <- gsub(".*\\|(.+)_ARATH", "\\1", DEPs)

tair.mapping <- bitr(dep.names,
                     fromType = "SYMBOL",
                     toType   = "TAIR",
                     OrgDb    = org.At.tair.db)

gene.info.prot <- AnnotationDbi::select(org.At.tair.db,
                                        keys    = tair.mapping$TAIR,
                                        columns = c("GENENAME", "SYMBOL"),
                                        keytype = "TAIR")
gene.info.prot <- gene.info.prot[!duplicated(gene.info.prot$TAIR), ]

# Subcellular localization (GO Cellular Component)
go.cc.prot <- AnnotationDbi::select(org.At.tair.db,
                                    keys    = tair.mapping$TAIR,
                                    columns = c("GO", "ONTOLOGY"),
                                    keytype = "TAIR")
go.cc.prot <- go.cc.prot[!is.na(go.cc.prot$ONTOLOGY) & go.cc.prot$ONTOLOGY == "CC", ]

go.terms.prot <- AnnotationDbi::select(GO.db,
                                       keys    = unique(go.cc.prot$GO),
                                       columns = "TERM",
                                       keytype = "GOID")
go.cc.prot <- merge(go.cc.prot, go.terms.prot, by.x = "GO", by.y = "GOID", all.x = TRUE)

go.cc.collapsed.prot <- aggregate(TERM ~ TAIR, data = go.cc.prot,
                                  FUN = function(x) paste(unique(x), collapse = "; "))
colnames(go.cc.collapsed.prot) <- c("TAIR", "Subcellular_localization")

# Build the supplementary table for a given direction (Up/Down)
build.prot.table <- function(proteins, direction) {
  prot.names  <- gsub(".*\\|(.+)_ARATH", "\\1", proteins)
  uniprot.ids <- gsub("^[a-z]+\\|(.+)\\|.+", "\\1", proteins)
  
  df <- data.frame(
    Protein_ID   = proteins,
    UniProt_ID   = uniprot.ids,
    Protein_name = prot.names,
    Direction    = direction,
    log2FC       = round(results_table[proteins, "logFC"], 3),
    P.Value      = round(results_table[proteins, "P.Value"], 6),
    adj.P.Val    = round(results_table[proteins, "adj.P.Val"], 6),
    stringsAsFactors = FALSE
  )
  
  df <- merge(df, tair.mapping, by.x = "Protein_name", by.y = "SYMBOL", all.x = TRUE)
  df <- merge(df, gene.info.prot[, c("TAIR", "GENENAME")], by = "TAIR", all.x = TRUE)
  df <- merge(df, go.cc.collapsed.prot, by = "TAIR", all.x = TRUE)
  df <- df[order(abs(df$log2FC), decreasing = TRUE), ]
  df
}

supp.up.prot   <- build.prot.table(up.proteins,   "Up")
supp.down.prot <- build.prot.table(down.proteins, "Down")

write_xlsx(list("Up_DEPs"   = supp.up.prot, "Down_DEPs" = supp.down.prot),
           path = "tables/DEPs_supplementary_table.xlsx")

# End of script
