#!/bin/bash

# User-defined variables
# Paths and directories should be set by the user or passed as arguments
RAW_DIR=${RAW_DIR:-"/path/to/raw_data"}        # Raw data directory (set as a variable)
OUT_DIR=${OUT_DIR:-"/path/to/trimmed_data"}    # Output directory for trimmed files
REPORTS_DIR=${REPORTS_DIR:-"/path/to/reports"} # Directory for storing reports
REF=${REF:-"/path/to/reference/bowtie2_index"} # Reference genome (Bowtie2 index)
PICARD_PATH=${PICARD_PATH:-"/path/to/picard.jar"}  # Path to Picard tools
BAM_DIR=${BAM_DIR:-"/path/to/bam_files"}      # Directory for storing BAM files
NORMALIZED_DIR=${NORMALIZED_DIR:-"/path/to/normalized_coverages"} # Directory for normalized coverages
MACS2_PEAKS_DIR=${MACS2_PEAKS_DIR:-"/path/to/macs2_peaks"} # Directory for storing peak calling results

# Load required modules
module add trim_galore-0.6.2_py3
module add fastQC-0.11.5
module add python27-modules-gcc
module add bowtie2-2.4.2
module add samtools-1.11
module add picard-2.9.0
module add bamCompare # (assuming this is a separate module)

# Step 1: Trim raw data
echo "Starting trimming process..."
for R1_FASTQ in $RAW_DIR/*_R1_001.fastq.gz; do
    BASE_NAME=$(basename $R1_FASTQ "_R1.fastq.gz")
    R2_FASTQ="${RAW_DIR}/${BASE_NAME}_R2.fastq.gz"
    
    # Run Trim Galore with FastQC
    trim_galore --paired -q 25 -j 4 --fastqc $R1_FASTQ $R2_FASTQ -o $OUT_DIR
done

# Step 2: Map trimmed datasets using Bowtie2 and process with Samtools
echo "Starting mapping and sorting with Bowtie2..."
for TRIMMED_R1 in $OUT_DIR/*_R1_val_1.fq.gz; do
    BASE_NAME=$(basename $TRIMMED_R1 "_R1_val_1.fq.gz")
    TRIMMED_R2="${OUT_DIR}/${BASE_NAME}_R2_val_2.fq.gz"
    
    # Run Bowtie2 mapping
    bowtie2 -p 20 --no-mixed --no-discordant --very-sensitive-local -X 1000 -x $REF -1 $TRIMMED_R1 -2 $TRIMMED_R2 2>${REPORTS_DIR}/${BASE_NAME}_bowtie2_report.txt | \
    samtools view -uh -f3 -F 1024 -@20 - | \
    samtools sort -@ 20 - -o $BAM_DIR/${BASE_NAME}_sorted.bam
done

# Step 3: Index BAM files
echo "Indexing BAM files..."
for BAM_FILE in $BAM_DIR/*.bam; do
    samtools index -c $BAM_FILE
done

# Step 4: Remove duplicates using Picard
echo "Removing duplicates using Picard..."
for BAM_FILE in $BAM_DIR/*.bam; do
    java -jar $PICARD_PATH MarkDuplicates \
        I=$BAM_FILE \
        O=${BAM_FILE%%.bam}_noDups.bam \
        M=${BAM_FILE%%.bam}_noDups_stats.txt \
        REMOVE_DUPLICATES=true
done

# Step 5: Generate input-normalized coverages using bamCompare where -b1 is a ChIP-sample and -b2 is the corresponding Input
echo "Generating normalized coverages..."
bamCompare -b1 $BAM_DIR/*_noDups.bam \
           -b2 $BAM_DIR/*input_noDups.bam \
           --centerReads -e --effectiveGenomeSize 5300000000 --binSize 30 \
           --outFileName $NORMALIZED_DIR/coverage.bw \
           --outFileFormat bigwig -p 20

# Step 6: Peak calling using MACS2 where -t is a ChIP-sample and -c is the corresponding Input
echo "Calling peaks using MACS2..."
macs2 callpeak -t $BAM_DIR/*_noDups.bam \
               -c $BAM_DIR/*input_noDups.bam \
               --broad -g 5.30e+09 --broad-cutoff 0.1 \
               -n sample_name_peaks --outdir $MACS2_PEAKS_DIR

echo "Pipeline completed successfully!"

