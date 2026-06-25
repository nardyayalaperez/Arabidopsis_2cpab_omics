# Arabidopsis_2cpab_omics

This repository contains the analysis scripts used in the Master's Thesis (TFM) "Redox control of lipid catabolism at early stages of plant development", investigating the impact of 2-Cys peroxiredoxin (2-Cys Prx) deficiency on lipid catabolism in Arabidopsis thaliana seedlings, using complementary transcriptomic, proteomic and fatty acid composition analyses.

Authors: Nardy Celeste Ayala-Pérez, Fernando Rodríguez Marín, Francisco Javier Cejudo, Juan Manuel Pérez-Ruiz, María Luisa Hernández.

**Overview**

The *2cpab* mutant, lacking both 2-Cys peroxiredoxin isoforms, was compared against wild-type (WT) *Arabidopsis thaliana* seedlings across three complementary datasets:

- Transcriptomics (RNA-seq): differential gene expression analysis using DESeq2 (primary) and limma-voom (comparison), followed by GO and KEGG functional enrichment and cross-referencing with a curated lipid metabolism gene list.
- Proteomics (SWATH-MS): differential protein abundance analysis using limma, followed by GO functional enrichment.
- Fatty acid composition (FAMEs): relative abundance of fatty acid species compared between genotypes using Welch's t-test.

## Repository structure and analytical pipeline
 
The scripts used for the different analyses are organized into three main directories: "Transcriptomics", "Proteomics", and "FAMEs". The complete analytical pipeline developed for this Master's Thesis can be explored in the following sequential steps:

1. [RNA-seq data processing pipeline](Transcriptomics/Bash_code.md)
2. [RNA-seq analysis: Differential gene expression (DESeq2 + limma-voom)](Transcriptomics/RNAseq_01_DEGAnalysis.R)
3. [RNA-seq analysis: GO and KEGG functional enrichment](Transcriptomics/RNAseq_02_FunctionalEnrichment.R)
4. [RNA-seq analysis: Cross-reference with lipid metabolism genes](Transcriptomics/RNAseq_03_LipidGeneCross.R)
5. [Proteomics: Differential protein abundance analysis (limma)](Proteomics/Proteomics_01_DEPAnalysis.R)
6. [Proteomics: GO functional enrichment](Proteomics/Proteomics_02_FunctionalEnrichment.R)
7. [Fatty acid composition: Welch's t-test between genotypes](Fatty_acid_composition/FAMEs_01_Analysis.R)

All R scripts expect to be run from a working directory containing the corresponding raw data files (see the header of each script for the exact
input file names). Each script automatically creates `tables/` and `images/` subdirectories to store its outputs.
