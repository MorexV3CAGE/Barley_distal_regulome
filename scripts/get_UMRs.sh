#!/bin/bash

# Load required modules
module load bedtools-2.26.0
module load python36-modules-gcc

# ==== Set file paths ====
GENOME_FILE="genome.txt"                       # Tab-separated: chrom <tab> length
METHYLATION_FILE="methylation_sorted.txt"      # Sorted BED-like file with methylation values in column 4
BINS_OUTPUT="genome_binned_100bp.bed"
METHYLATION_MAPPED="methylation_mapped_means.bed"
FILTERED_BINS="UMRs_filtered.bed"
FINAL_UMRs="UMR_final.bed"

# ==== Step 1: Create 100bp bins across the whole genome ====
bedtools makewindows -g "$GENOME_FILE" -w 100 > "$BINS_OUTPUT"

# ==== Step 2: Map methylation data onto the bins (mean methylation per bin) ====
bedtools map -a "$BINS_OUTPUT" -b "$METHYLATION_FILE" -c 4 -o mean > "$METHYLATION_MAPPED"

# ==== Step 3: Run R for filtering low-methylation bins (≤1%) ====
Rscript - <<EOF
library(dplyr)

# Load methylation-mapped bins
UMRs <- read.delim("$METHYLATION_MAPPED", header = FALSE, sep = "\t", quote = "\"", fill = TRUE)

# Step 3a: Filter for bins with ≤1% methylation and remove unplaced contigs (optional)
UMRs_true <- UMRs[UMRs\$V4 <= 1 & UMRs\$V1 != "chrUn", ]

# Step 3b: Remove isolated missing-data bins (".") unless flanked by real values
UMRs_filtered <- UMRs_true %>%
  filter(!(V4 == "." &
           (lead(V4, default = "") == "." |
            lag(V4, default = "" ) == "."))) %>%
  na.omit()

# Step 3c: Further filter bins that are not continuous with neighbors
UMRs_final <- UMRs_filtered %>%
  mutate(prev_end = lag(V3),
         next_start = lead(V2)) %>%
  filter(!(V4 == "." &
           (prev_end != V2 | next_start != V3)))

# Save intermediate result
write.table(UMRs_final, file="$FILTERED_BINS", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
EOF

# ==== Step 4: Merge adjacent bins into continuous UMRs ====
bedtools merge -i "$FILTERED_BINS" -c 4 -o mean > "$FINAL_UMRs"

# ==== Step 5: Keep only UMRs ≥300bp ====
awk '{if (($3 - $2) >= 300) print}' "$FINAL_UMRs" > tmp && mv tmp "$FINAL_UMRs"

echo "UMR identification complete. Final output: $FINAL_UMRs"