#!/bin/bash

# Define directories
RAW_DATA_DIR="/path/to/raw/data"
QC_DIR="/path/to/quality_control"
ASSEMBLY_DIR="/path/to/assembly"
GENE_PREDICTION_DIR="/path/to/gene_prediction"
ANNOTATION_DIR="/path/to/annotation"


mkdir -p ${QC_DIR}
mkdir -p ${ASSEMBLY_DIR}
mkdir -p ${GENE_PREDICTION_DIR}
mkdir -p ${ANNOTATION_DIR}

## 0. Pre-processing of metatranscriptomic reads
# Cleaning and quality control of reads using fastp
for f in ${RAW_DATA_DIR}/*_R1.fastq.gz; do
    base=$(basename ${f} _R1.fastq.gz)
    fastp -i ${RAW_DATA_DIR}/${base}_R1.fastq.gz -I ${RAW_DATA_DIR}/${base}_R2.fastq.gz \
          -o ${QC_DIR}/${base}_R1_trimmed.fastq.gz -O ${QC_DIR}/${base}_R2_trimmed.fastq.gz \
          --cut_by_quality3 30 --cut_by_quality5 30 --thread 40 \
          --html ${QC_DIR}/${base}_fastp_report.html --json ${QC_DIR}/${base}_fastp_report.json
done

## 1. Assembly using Megahit
# Assemble the cleaned and quality-controlled reads into contigs
megahit -1 ${QC_DIR}/*_R1_trimmed.fastq.gz -2 ${QC_DIR}/*_R2_trimmed.fastq.gz \
        -o ${ASSEMBLY_DIR} \
        --k-min 31 --k-max 91 --k-step 10 --min-contig-len 300 -t 40

## 2. Gene prediction with Prodigal
# Predict genes on assembled contigs
for contig in ${ASSEMBLY_DIR}/final.contigs.fa; do
    prodigal -i $contig -a ${GENE_PREDICTION_DIR}/$(basename $contig .fa)_proteins.faa \
             -d ${GENE_PREDICTION_DIR}/$(basename $contig .fa)_genes.fna -p meta \
             -o ${GENE_PREDICTION_DIR}/$(basename $contig .fa).gbk
done

## 3. Annotation in accordance with the Joint Genome Institute’s protocols
# Mapping to KEGG orthologs for functional insights
for protein in ${GENE_PREDICTION_DIR}/*.faa; do
    blastp -query $protein \
           -db /path/to/blast/db/kegg_proteins \
           -out ${ANNOTATION_DIR}/$(basename $protein .faa)_annotated.txt \
           -outfmt 6 \
           -evalue 1e-5 \
           -num_threads 4 \
           -max_target_seqs 1
done

