###############################################################################
# Proteomics Analysis - Functional Enrichment Analysis (GO)
#
# Description: Gene Ontology (Biological Process) enrichment analysis of
# differentially abundant proteins (DEPs) identified in the 2cpab vs WT
# comparison (SWATH-MS proteomics), performed separately for up and
# downregulated proteins. GO terms and their adjusted p-values are exported
# in the format required by REVIGO, used to group enriched terms by semantic
# similarity and generate a more structured visualization.
#
# Organism: Arabidopsis thaliana
#
# Input:
#   tables/full_results_proteomics.csv : complete limma results table from
#                                        Proteomics_01_DEPAnalysis.R
#
# Output:
#   tables/GO_proteomics_up.tsv     : full GO enrichment table, up-regulated
#   tables/GO_proteomics_down.tsv   : full GO enrichment table, down-regulated
#   tables/revigo_prot_up.txt       : GO IDs + p.adjust, formatted for REVIGO
#   tables/revigo_prot_down.txt     : GO IDs + p.adjust, formatted for REVIGO
#
# Author: Nardy Celeste Ayala Pérez
# Date: 2026
###############################################################################

#-------------------------------------------------------------------------------
# PACKAGE LOADING
#-------------------------------------------------------------------------------
library(clusterProfiler)  # enrichGO, bitr
library(org.At.tair.db)   # Arabidopsis thaliana genome annotation

dir.create("tables", showWarnings = FALSE)

#-------------------------------------------------------------------------------
# LOAD RESULTS TABLE FROM Proteomics_01_DEPAnalysis.R OUTPUT AND DEFINE
# THE SIGNIFICANT DEP SETS (adj.P.Val < 0.05)
#-------------------------------------------------------------------------------

results_table <- read.csv("tables/full_results_proteomics.csv", row.names = 1)

FC_threshold <- log2(1.5)

up.proteins.sig   <- rownames(results_table[results_table$logFC > FC_threshold &
                                              results_table$adj.P.Val < 0.05, ])
down.proteins.sig <- rownames(results_table[results_table$logFC < -FC_threshold &
                                              results_table$adj.P.Val < 0.05, ])

cat("Up-regulated proteins (adj.P.Val < 0.05):", length(up.proteins.sig))
cat("Down-regulated proteins (adj.P.Val < 0.05):", length(down.proteins.sig))

# Extract gene symbols from UniProt-style protein IDs
# (e.g., sp|Q9C5R8|BAS1B_ARATH -> BAS1B)
up.names   <- gsub(".*\\|(.+)_ARATH", "\\1", up.proteins.sig)
down.names <- gsub(".*\\|(.+)_ARATH", "\\1", down.proteins.sig)

# BAS1B is not recognized as a valid SYMBOL by org.At.tair.db and must be
# renamed to BAS1 for successful TAIR mapping
down.names <- gsub("BAS1B", "BAS1", down.names)

#-------------------------------------------------------------------------------
# GO ENRICHMENT - UPREGULATED PROTEINS (adj.P.Val < 0.05)
#-------------------------------------------------------------------------------

up.tair <- bitr(up.names,
                fromType = "SYMBOL",
                toType   = "TAIR",
                OrgDb    = org.At.tair.db)

ego.prot.up <- enrichGO(gene          = up.tair$TAIR,
                        OrgDb         = org.At.tair.db,
                        keyType       = "TAIR",
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05)

cat("GO terms enriched (up-regulated proteins):", nrow(as.data.frame(ego.prot.up)))

write.table(as.data.frame(ego.prot.up), file = "tables/GO_proteomics_up.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

#-------------------------------------------------------------------------------
# GO ENRICHMENT - DOWNREGULATED PROTEINS (adj.P.Val < 0.05)
#-------------------------------------------------------------------------------

down.tair <- bitr(down.names,
                  fromType = "SYMBOL",
                  toType   = "TAIR",
                  OrgDb    = org.At.tair.db)

ego.prot.down <- enrichGO(gene          = down.tair$TAIR,
                          OrgDb         = org.At.tair.db,
                          keyType       = "TAIR",
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05)

cat("GO terms enriched (down-regulated proteins):", 
    nrow(as.data.frame(ego.prot.down)))

write.table(as.data.frame(ego.prot.down), 
            file = "tables/GO_proteomics_down.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

#-------------------------------------------------------------------------------
# EXPORT GO TERMS FOR REVIGO
#-------------------------------------------------------------------------------

# GO IDs and adjusted p-values are exported in the tab-separated format
# required by REVIGO (http://revigo.irb.hr/), used to group enriched
# terms by semantic similarity and generate a more structured
# visualization than a standard dotplot.

revigo.up <- data.frame(
  GO_ID  = ego.prot.up@result$ID,
  pvalue = ego.prot.up@result$p.adjust
)
revigo.up <- revigo.up[revigo.up$pvalue < 0.05, ]
write.table(revigo.up, "tables/revigo_prot_up.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")

revigo.down <- data.frame(
  GO_ID  = ego.prot.down@result$ID,
  pvalue = ego.prot.down@result$p.adjust
)
revigo.down <- revigo.down[revigo.down$pvalue < 0.05, ]
write.table(revigo.down, "tables/revigo_prot_down.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")

# End of script