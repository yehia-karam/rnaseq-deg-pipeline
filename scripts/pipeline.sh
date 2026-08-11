#!/bin/bash

# NGS Analysis Pipeline
# This script performs Alignment and Gene Counting

# --- Configuration ---
REF_GENOME="reference_genome.fa"
REF_ANNOTATION="annotation.gtf"
INDEX_PREFIX="genome_index"
SAMPLES=(
    "HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22"
    "HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22"
    "HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22"
    "UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22"
    "UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22"
    "UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22"
)

# Create directories
mkdir -p alignment counts

echo "Step 1: Building HISAT2 Index..."
hisat2-build -p 4 $REF_GENOME $INDEX_PREFIX

echo "Step 2: Aligning reads..."
for sample in "${SAMPLES[@]}"; do
    echo "Aligning $sample..."
    # Assume trimmed files are in ../trimmed_data/
    R1="../trimmed_data/${sample}_R1_trimmed.fastq.gz"
    R2="../trimmed_data/${sample}_R2_trimmed.fastq.gz"
    
    hisat2 -p 4 -x $INDEX_PREFIX -1 $R1 -2 $R2 -S alignment/${sample}.sam
    
    echo "Converting to BAM and sorting..."
    samtools view -bS alignment/${sample}.sam | samtools sort -o alignment/${sample}_sorted.bam
    samtools index alignment/${sample}_sorted.bam
    rm alignment/${sample}.sam
done

echo "Step 3: Generating gene counts..."
featureCounts -p -a $REF_ANNOTATION -o counts/gene_counts.txt alignment/*_sorted.bam

echo "Pipeline Complete!"
