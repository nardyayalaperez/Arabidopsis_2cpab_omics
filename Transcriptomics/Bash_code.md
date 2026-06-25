# **RNA-seq Preprocessing: Reference Preparation, Alignment and Quantification**

This document describes the bash/cluster pipeline used to process raw RNA-seq FASTQ files into transcript-level TPM quantification tables, prior to the downstream differential expression analysis performed in R (RNAseq_01_DEGAnalysis.R).

## **1. Reference Data and Preparation**

We download and decompress the Arabidopsis thaliana (TAIR10) reference genome and structural annotation files from Ensembl Plants.

```
wget -O arabidopsis_thaliana.fa.gz https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-62/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

wget -O arabidopsis_thaliana.gtf.gz https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-62/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gtf.gz

gunzip arabidopsis_thaliana.fa.gz
gunzip arabidopsis_thaliana.gtf.gz

````
## 2. Merging Multi-Lane FASTQ Files

Each sample was sequenced across 4 lanes (L001-L004). Reads from all lanes belonging to the same biological replicate are concatenated into a single FASTQ file per sample.

```
cat GENext102-21-1_S1_L001_R1_001.fastq.gz \
    GENext102-21-1_S1_L002_R1_001.fastq.gz \
    GENext102-21-1_S1_L003_R1_001.fastq.gz \
    GENext102-21-1_S1_L004_R1_001.fastq.gz > WT_21-1.fastq.gz

cat GENext102-21-2_S2_L001_R1_001.fastq.gz \
    GENext102-21-2_S2_L002_R1_001.fastq.gz \
    GENext102-21-2_S2_L003_R1_001.fastq.gz \
    GENext102-21-2_S2_L004_R1_001.fastq.gz > WT_21-2.fastq.gz

cat GENext102-21-3_S3_L001_R1_001.fastq.gz \
    GENext102-21-3_S3_L002_R1_001.fastq.gz \
    GENext102-21-3_S3_L003_R1_001.fastq.gz \
    GENext102-21-3_S3_L004_R1_001.fastq.gz > WT_21-3.fastq.gz

cat GENext102-21-4_S4_L001_R1_001.fastq.gz \
    GENext102-21-4_S4_L002_R1_001.fastq.gz \
    GENext102-21-4_S4_L003_R1_001.fastq.gz \
    GENext102-21-4_S4_L004_R1_001.fastq.gz > 2cpab_21-4.fastq.gz

cat GENext102-21-5_S5_L001_R1_001.fastq.gz \
    GENext102-21-5_S5_L002_R1_001.fastq.gz \
    GENext102-21-5_S5_L003_R1_001.fastq.gz \
    GENext102-21-5_S5_L004_R1_001.fastq.gz > 2cpab_21-5.fastq.gz

cat GENext102-21-6_S6_L001_R1_001.fastq.gz \
    GENext102-21-6_S6_L002_R1_001.fastq.gz \
    GENext102-21-6_S6_L003_R1_001.fastq.gz \
    GENext102-21-6_S6_L004_R1_001.fastq.gz > 2cpab_21-6.fastq.gz
````
## 3. Genome Index Generation

We build a STAR genome index using the TAIR10 reference genome and GTF file.

```
STAR --runMode genomeGenerate \
     --genomeDir genome/index \
     --genomeFastaFiles arabidopsis_thaliana.fa \
     --sjdbGTFfile arabidopsis_thaliana.gtf \
     --genomeSAindexNbases 10
```
## 4. Sample Processing Script

We define a bash script, process_sample.sh, that performs quality control (FastQC), read alignment (STAR) and transcript assembly and quantification (StringTie) for a single sample.

```
#!/bin/bash
FOLDER=$1
NAME=$2

## Access the sample folder
cd $FOLDER

## Sample quality control and read mapping to reference genome
module load FastQC/0.12.1-Java-11
fastqc ${NAME}.fastq.gz

## Read mapping
module load STAR/2.7.11b-GCC-12.3.0
STAR --genomeDir ../genome/index/ \
     --readFilesIn ${NAME}.fastq.gz \
     --readFilesCommand "gunzip -c" --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif \
     --outFilterIntronMotifs RemoveNoncanonical --alignIntronMax 300000 \
     --outFileNamePrefix $NAME

## Transcript assembly and Quantification
source /lustre/software/easybuild/common/software/Miniconda3/4.9.2/etc/profile.d/conda.sh
conda activate stringtie-3.0.3

stringtie -G ../genome/arabidopsis_thaliana.gtf -o $NAME.gtf \
          ${NAME}Aligned.sortedByCoord.out.bam

stringtie -e -G ../genome/arabidopsis_thaliana.gtf \
          -o $NAME.gtf -A $NAME.tsv \
          ${NAME}Aligned.sortedByCoord.out.bam

conda deactivate
```
## 5. Job Submission

Each sample is processed by submitting process_sample.sh as a cluster job, passing the sample folder and sample name as arguments.

```
sbatch --job-name=WT_21-1 --output=WT_21-1.out process_sample.sh samples/WT_21-1/ WT_21-1
sbatch --job-name=WT_21-1 --output=WT_21-1.out process_sample.sh samples/WT_21-3/ WT_21-1
sbatch --job-name=WT_21-1 --output=WT_21-1.out process_sample.sh samples/WT_21-3/ WT_21-1
sbatch --job-name=2pcab_21-4 --output=2pcab_21-4.out process_sample.sh samples/2pcab_21-4/ 2cpab_21-4
sbatch --job-name=2pcab_21-4 --output=2pcab_21-4.out process_sample.sh samples/2pcab_21-4/ 2cpab_21-5
sbatch --job-name=2pcab_21-4 --output=2pcab_21-4.out process_sample.sh samples/2pcab_21-4/ 2cpab_21-6

````

## 6. Gene-Level Count Matrix Generation

We generate a raw gene-level count matrix (gene_count_matrix.csv) for discrete expression analysis using the prepDE.py script distributed with StringTie.

```
#!/bin/bash
source /lustre/software/easybuild/common/software/Miniconda3/4.9.2/etc/profile.d/conda.sh
conda activate stringtie-3.0.3

cd /home/nayala/rnaseq
prepDE.py -i samples/

conda deactivate
```
Submitted as a cluster job:

```
sbatch --job-name=prepDE --output=prepDE.out --time=06:00:00 scripts/prep.sh
```




