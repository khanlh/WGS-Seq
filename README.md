## Download and Convert SRA Files to Compressed FASTQ

> **Note:**  
> This script demonstrates how to mine raw sequencing data from the NCBI Sequence Read Archive (SRA) by downloading `.sra` files, converting them to FASTQ format, and compressing the results for downstream analysis.

### Script Overview

```bash
#!/usr/bin/env bash

# Create output directory
mkdir -p fastq

# Function to process a single SRR ID
process_srr() {
    i=$1
    echo "Downloading SRR${i}.sra..."
    wget -q https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR${i}/SRR${i} -O SRR${i}.sra

    echo "Converting SRR${i}.sra to FASTQ..."
    fastq-dump --split-files --outdir ./fastq SRR${i}.sra

    echo "Compressing FASTQ files..."
    gzip fastq/SRR${i}_*.fastq

    echo "Cleaning up SRR${i}.sra..."
    rm SRR${i}.sra

    echo "Done SRR${i}"
    echo "-------------------------"
}

export -f process_srr

# Run up to 4 processes in parallel for SRR12287823 to SRR12287836
seq 12287823 12287836 | xargs -n 1 -P 4 -I {} bash -c 'process_srr "$@"' _ {}
```


# WGS-Seq

```bash
#!/bin/bash

# ------------------ CONFIGURATION ------------------

# Reference genome FASTA file (used for indexing and variant calling)
REF_GENOME="/media/khanle/Untitled/Monodon/P.mon_ref/GCF_015228065.2_NSTDA_Pmon_1_genomic.fna"

# Data directories
INPUT_DIR="/media/khanle/Untitled/Monodon/DNA_seq/PRJNA1056855/fastq"
TRIM_DIR="$INPUT_DIR/../trimmed"
MAP_DIR="$INPUT_DIR/../mapped"
VCF_DIR="$INPUT_DIR/../vcf"

# Activate Conda environment containing bwa, samtools, gatk
source /home/khanle/miniconda3/bin/activate
conda activate dna_env

# Create directories if they do not exist
mkdir -p "$TRIM_DIR" "$MAP_DIR" "$VCF_DIR"

# ------------------ BWA: Index the reference genome (run once) ------------------

if [[ ! -f "${REF_GENOME}.bwt" ]]; then
    echo "Indexing the reference genome with BWA..."
    bwa index "$REF_GENOME"
fi

# ------------------ BWA MEM: Mapping paired-end reads ------------------

for file in "$TRIM_DIR"/*_1_paired.fastq.gz; do
    base=$(basename "$file" _1_paired.fastq.gz)
    R1="$TRIM_DIR/${base}_1_paired.fastq.gz"
    R2="$TRIM_DIR/${base}_2_paired.fastq.gz"
    OUT_BAM="$MAP_DIR/${base}_sorted.bam"

    # Check if R2 file exists
    if [[ ! -f "$R2" ]]; then
        echo "Warning: Missing R2 file for sample $base, skipping..."
        continue
    fi

    echo "Mapping sample: $base"

    # Align reads with BWA MEM, convert to BAM, and sort
    bwa mem -t 4 "$REF_GENOME" "$R1" "$R2" | \
        samtools view -bS - | \
        samtools sort -o "$OUT_BAM"

    # Index the sorted BAM file
    samtools index "$OUT_BAM"
done

# ------------------ GATK: Variant Calling ------------------

# Create FASTA index file if not present
if [[ ! -f "${REF_GENOME}.fai" ]]; then
    echo "Creating FASTA index with samtools..."
    samtools faidx "$REF_GENOME"
fi

# Create sequence dictionary for GATK if not present
DICT="${REF_GENOME%.*}.dict"
if [[ ! -f "$DICT" ]]; then
    echo "Creating sequence dictionary with GATK..."
    gatk CreateSequenceDictionary -R "$REF_GENOME" -O "$DICT"
fi

# Call variants for each BAM file using GATK HaplotypeCaller
for BAM_FILE in "$MAP_DIR"/*.bam; do
    SAMPLE=$(basename "$BAM_FILE" _sorted.bam)
    echo "Calling variants for sample: $SAMPLE"

    gatk HaplotypeCaller \
        -R "$REF_GENOME" \
        -I "$BAM_FILE" \
        -O "$VCF_DIR/${SAMPLE}.vcf.gz"
done

echo "Pipeline completed successfully."
```

