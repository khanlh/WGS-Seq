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
