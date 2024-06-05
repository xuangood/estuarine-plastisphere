#!/bin/bash

## 0. PRE-PROCESSING OF METATRANSCRIPTOMIC READS
# For each fastq file, the following commands were run to trim the reads then de-replicate them if needed
for f in /path/to/raw/data/*_R1.fastq.gz; do
    base=$(basename $f _R1.fastq.gz)
    fastp -i /path/to/raw/data/${base}_R1.fastq.gz -I /path/to/raw/data/${base}_R2.fastq.gz -o /path/to/qc/${base}_R1_trimmed.fastq.gz -O /path/to/qc/${base}_R2_trimmed.fastq.gz --cut_by_quality3 30 --cut_by_quality5 30 --thread 40
    usearch -fastq_mergepairs /path/to/qc/${base}_R1_trimmed.fastq.gz -relabel @ -fastqout /path/to/qc/${base}_merged.fastq -fastq_minlen 130 -threads 40
    usearch -fastx_uniques /path/to/qc/${base}_merged.fastq -fastaout /path/to/qc/${base}_unique.fastq -sizeout -threads 40
done

## 1. ASSEMBLY USING MEGAH. IT
mkdir -p /path/to/assembly
megahit -r /path/to/qc/*_unique.fastq -o /path/to/assembly --k-min 31 --k-max 91 --k-step 10 --min-contig-len 300 -t 40

## 2. GENE PREDICTION WITH PRODIGAL
mkdir -p /path/to/gene_prediction
for contig in /path/to/assembly/*.fa; do
    prodigal -i $contig -a /path/to/gene_prediction/$(basename $contig .fa)_proteins.faa -d /path/to/gene_prediction/$(basename $contig .fa)_genes.fna -p meta -o /path/to/gene_prediction/$(basename $contig .fa).gbk
done

## 3. ANNOTATION USING JGI'S PROTOCOL
mkdir -p /path/to/annotation
for protein in /path/to/gene_prediction/*.faa; do
    run_annotation --input $protein --output /path/to/annotation/$(basename $protein .faa)_annotated --db KEGG
done

echo "Metatranscriptomic analysis pipeline completed successfully!"
