#!/bin/bash

# Bismark BS-seq Analysis Pipeline
# Usage: ./bismark_pipeline.sh <trimmed_reads_dir> <genome_folder> <output_base_dir> <sample_prefix> [threads]

# Required inputs
TRIMMED_DIR=$1           # Directory with trimmed *.fq.gz files
GENOME_DIR=$2            # Path to genome directory prepared for Bismark
OUTPUT_BASE=$3           # Base directory for all outputs
SAMPLE_PREFIX=$4         # Prefix identifying sample (e.g., SRR5124893)
THREADS=${5:-8}          # Optional: default to 8 threads if not specified

# Load module
module load bismark/0.23.1

# Prepare genome (once per genome)
echo "Preparing Bismark genome..."
bismark_genome_preparation --bowtie2 --verbose "$GENOME_DIR" 2> "$OUTPUT_BASE/genome_prep_summary.txt"

# Mapping reads
echo "Mapping reads..."
for fq1 in "$TRIMMED_DIR"/*_1_val_1.fq.gz; do
    base=$(basename "$fq1" "_1_val_1.fq.gz")
    fq2="$TRIMMED_DIR/${base}_2_val_2.fq.gz"
    OUTMAP="$OUTPUT_BASE/mapped/${base}_run"
    mkdir -p "$OUTMAP"
    bismark -o "$OUTMAP" --parallel "$THREADS" --genome "$GENOME_DIR" -1 "$fq1" -2 "$fq2" \
        2> "$OUTPUT_BASE/mapped/${base}_mapping_summary.txt"
done

# Deduplicate
echo "Deduplicating..."
DEDUP_IN="$OUTPUT_BASE/mapped/${SAMPLE_PREFIX}_run/${SAMPLE_PREFIX}_1_val_1_bismark_bt2_pe.bam"
DEDUP_OUT="$OUTPUT_BASE/deduplicated"
mkdir -p "$DEDUP_OUT"
deduplicate_bismark "$DEDUP_IN" --output_dir "$DEDUP_OUT"

# Methylation extraction
echo "Extracting methylation..."
EXTRACTION_OUT="$OUTPUT_BASE/extraction"
mkdir -p "$EXTRACTION_OUT"
bismark_methylation_extractor \
    --multicore "$THREADS" \
    --bedGraph --CX --cytosine_report --comprehensive --gzip \
    "$DEDUP_OUT/${SAMPLE_PREFIX}_1_val_1_bismark_bt2_pe.deduplicated.bam" \
    -o "$EXTRACTION_OUT" \
    --genome_folder "$GENOME_DIR" \
    2> "$EXTRACTION_OUT/bismark_extraction_summary.txt"

# Cytosine report
echo "Generating cytosine content report..."
coverage2cytosine \
    "$EXTRACTION_OUT/${SAMPLE_PREFIX}_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz" \
    --CX --genome_folder "$GENOME_DIR" \
    -o "$EXTRACTION_OUT/cytosine_content.txt"

# Bismark reports
echo "Creating Bismark summary report..."
REPORT_DIR="$OUTPUT_BASE/report"
mkdir -p "$REPORT_DIR"
bismark2report --dir "$REPORT_DIR" \
    --alignment_report "$OUTPUT_BASE/mapped/${SAMPLE_PREFIX}_run/${SAMPLE_PREFIX}_1_val_1_bismark_bt2_PE_report.txt" \
    --dedup_report "$DEDUP_OUT/${SAMPLE_PREFIX}_1_val_1_bismark_bt2_pe.deduplication_report.txt" \
    --splitting_report "$EXTRACTION_OUT/${SAMPLE_PREFIX}_1_val_1_bismark_bt2_pe.deduplicated_splitting_report.txt" \
    --mbias_report "$EXTRACTION_OUT/${SAMPLE_PREFIX}_1_val_1_bismark_bt2_pe.deduplicated.M-bias.txt" \
    --nucleotide_report "none"

# Create bedGraph
echo "Generating bedGraph files..."
BEDGRAPH_DIR="$OUTPUT_BASE/bedgraphs"
mkdir -p "$BEDGRAPH_DIR"
cd "$EXTRACTION_OUT"
bismark2bedGraph --buffer 10G -dir "$BEDGRAPH_DIR" -o "${SAMPLE_PREFIX}_CpG_context.bedGraph.gz" CpG*
bismark2bedGraph --buffer 10G -dir "$BEDGRAPH_DIR" -o "${SAMPLE_PREFIX}_CHG_context.bedGraph.gz" CHG*
bismark2bedGraph --CX --ample_mem -dir "$BEDGRAPH_DIR" -o "${SAMPLE_PREFIX}_CHH_context.bedGraph.gz" CHH*

echo "Pipeline complete!"
