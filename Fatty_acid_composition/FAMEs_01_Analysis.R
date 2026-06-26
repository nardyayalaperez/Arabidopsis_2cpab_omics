###############################################################################
# Fatty Acid Methyl Esters (FAMEs) Analysis - WT vs 2cpab
#
# Description: Analysis of fatty acid composition in seven-day-old Arabidopsis
# thaliana seedlings, comparing the 2cpab mutant against wild-type (WT).
# Relative abundance (%) of each fatty acid species is compared between
# genotypes using Welch's two-sample t-test, with Benjamini-Hochberg (BH)
# correction for multiple comparisons. Results are visualized as a bar plot
# showing the relative abundance of each fatty acid across samples.
#
# Organism: Arabidopsis thaliana
# Replicates: 5 WT and 5 2cpab (10 samples total)
#
# Input:
#   fames.csv : raw FAMEs data, semicolon-separated, containing a "%" section
#               with the relative abundance of each fatty acid per sample
#
# Output:
#   tables/FAMEs_summary.csv      : mean and SD of each fatty acid per genotype
#   tables/FAMEs_statistics.csv   : summary table with Welch's t-test p-values
#                                   and BH-adjusted p-values
#   images/FAMEs_barplot.png      : bar plot of fatty acid composition (%)
#
# Author: Nardy Celeste Ayala Pérez
# Date: 2026
###############################################################################

#-------------------------------------------------------------------------------
# PACKAGE LOADING
#-------------------------------------------------------------------------------
library(car)  # leveneTest()
library(ggplot2)    # Bar plot and boxplot visualization
library(tidyverse)   # Data manipulation
library(reshape2)    # Data reshaping

dir.create("tables", showWarnings = FALSE)
dir.create("images", showWarnings = FALSE)

#-------------------------------------------------------------------------------
# LOAD RAW FAMEs DATA
#-------------------------------------------------------------------------------

# The raw file contains multiple data sections (absolute amounts,
# percentages, etc.). Only the "%" section, giving the relative abundance
# of each fatty acid, is used for this analysis.

data <- read.csv("fames.csv",
                 sep = ";", header = FALSE,
                 stringsAsFactors = FALSE,
                 check.names = FALSE)


#-------------------------------------------------------------------------------
# EXTRACT THE PERCENTAGE (%) DATA SECTION
#-------------------------------------------------------------------------------

# Locate the row marking the start of the "%" section
pct.start <- which(data[, 1] == "%")

# Extract the percentage matrix (rows: fatty acid species, columns: 10 samples)
pct.data <- data[(pct.start + 1):(pct.start + 14), 1:11]
colnames(pct.data) <- c("FA", paste0("S", 1:10))
rownames(pct.data) <- pct.data$FA
pct.data <- pct.data[, -1]

# Convert values to numeric, replacing the European decimal comma with a point
pct.matrix <- apply(pct.data, 2, function(x) as.numeric(gsub(",", ".", x)))
rownames(pct.matrix) <- rownames(pct.data)

# Remove the "Total" row
pct.matrix <- pct.matrix[rownames(pct.matrix) != "Total", ]

cat("Fatty acids detected:", nrow(pct.matrix))
cat("Samples:", ncol(pct.matrix))


#-------------------------------------------------------------------------------
# DEFINE EXPERIMENTAL GROUPS
#-------------------------------------------------------------------------------

condition <- c(rep("WT", 5), rep("2cpab", 5))
colnames(pct.matrix) <- paste0(condition, "_", 1:10)

wt.cols   <- grep("WT",    colnames(pct.matrix))
cpab.cols <- grep("2cpab", colnames(pct.matrix))


#-------------------------------------------------------------------------------
# CALCULATE MEAN AND SD PER GENOTYPE FOR EACH FATTY ACID
#-------------------------------------------------------------------------------

wt.mean   <- rowMeans(pct.matrix[, wt.cols],   na.rm = TRUE)
cpab.mean <- rowMeans(pct.matrix[, cpab.cols], na.rm = TRUE)
wt.sd     <- apply(pct.matrix[, wt.cols],   1, sd, na.rm = TRUE)
cpab.sd   <- apply(pct.matrix[, cpab.cols], 1, sd, na.rm = TRUE)

summary.df <- data.frame(
  FA        = rownames(pct.matrix),
  WT.mean   = round(wt.mean,   2),
  WT.sd     = round(wt.sd,     2),
  cpab.mean = round(cpab.mean, 2),
  cpab.sd   = round(cpab.sd,   2)
)

print(summary.df)
write.csv(summary.df, "tables/FAMEs_summary.csv", row.names = FALSE)

#-------------------------------------------------------------------------------
# NORMALITY AND HOMOGENEITY OF VARIANCE ASSESSMENT
#-------------------------------------------------------------------------------

# Shapiro-Wilk test is used to assess normality within each genotype group,
# and Levene's test is used to assess homogeneity of variance between
# groups, for each fatty acid species.

shapiro.wt   <- sapply(rownames(pct.matrix), function(fa)
  shapiro.test(pct.matrix[fa, wt.cols])$p.value)
shapiro.cpab <- sapply(rownames(pct.matrix), function(fa)
  shapiro.test(pct.matrix[fa, cpab.cols])$p.value)

levene.p <- sapply(rownames(pct.matrix), function(fa) {
  values <- c(pct.matrix[fa, wt.cols], pct.matrix[fa, cpab.cols])
  groups <- factor(condition)
  leveneTest(values ~ groups)$`Pr(>F)`[1]
})

summary.df$Shapiro.WT    <- round(shapiro.wt, 4)
summary.df$Shapiro.2cpab <- round(shapiro.cpab, 4)
summary.df$Levene.p      <- round(levene.p, 4)

summary.df

# Once normality and variance homogeneity are assessed, Welch's two-sample
# t-test is applied to each fatty acid species, as it does not require
# equal variances between groups and is robust to minor deviations from
# normality.
             
#-------------------------------------------------------------------------------
# STATISTICAL TESTING
#-------------------------------------------------------------------------------

# WELCH'S T-TEST PER FATTY ACID
# t.test() uses var.equal = FALSE by default, performing Welch's
# two-sample t-test, which does not assume equal variances between
# genotypes. P-values are corrected for multiple comparisons using the
# Benjamini-Hochberg (BH) method.

pvalues <- sapply(rownames(pct.matrix), function(fa) {
  t.test(pct.matrix[fa, wt.cols],
         pct.matrix[fa, cpab.cols])$p.value
})

summary.df$pvalue <- round(pvalues, 4)
summary.df$padj   <- round(p.adjust(pvalues, method = "BH"), 4)
summary.df$sig    <- ifelse(summary.df$padj < 0.05, "*", "")

cat("STATISTICAL RESULTS")
print(summary.df[, c("FA", "WT.mean", "cpab.mean", "pvalue", "padj", "sig")])

write.csv(summary.df, "tables/FAMEs_statistics.csv", row.names = FALSE)

#-------------------------------------------------------------------------------
# BAR PLOT: FATTY ACID COMPOSITION (%)
#-------------------------------------------------------------------------------

plot.df <- data.frame(
  FA        = rep(rownames(pct.matrix), 2),
  Mean      = c(wt.mean, cpab.mean),
  SD        = c(wt.sd,   cpab.sd),
  Condition = rep(c("WT", "2cpab"), each = nrow(pct.matrix))
)

# Order fatty acids by their original row order (chain length / saturation)
fa.order <- rownames(pct.matrix)
plot.df$FA <- factor(plot.df$FA, levels = fa.order)
plot.df$Condition <- factor(plot.df$Condition, levels = c("WT", "2cpab"))

barplot.fames <- ggplot(plot.df, aes(x = FA, y = Mean, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                position = position_dodge(width = 0.8),
                width = 0.25, linewidth = 0.5) +
  scale_fill_manual(values = c(WT = "dodgerblue3", "2cpab" = "firebrick3")) +
  theme_classic(base_size = 12) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
        plot.title   = element_text(face = "bold", hjust = 0.5),
        legend.title = element_text(size = 11)) +
  labs(title = "Fatty Acid Composition - WT vs 2cpab",
       x     = "Fatty Acid",
       y     = "Relative abundance (%)",
       fill  = "Genotype")

ggsave("images/FAMEs_barplot.png", barplot.fames, width = 10, height = 6, dpi = 300)

# End of script
