#!/bin/bash
# ========================================================
# Community Abundance Quantification Pipeline (CoverM)
# ========================================================
# Prerequisites: CoverM, Samtools, and GNU sed/awk installed

BASE_DIR="${BASE_DIR:-.}"
OUTPUT_DIR="${OUTPUT_DIR:-$BASE_DIR/abundance_results}"
MAGS_DIR="${MAGS_DIR:-$BASE_DIR/data/mags}"
VIRAL_FASTA="${VIRAL_FASTA:-$BASE_DIR/data/viral_catalog.fasta}"
READS_DIR="${READS_DIR:-$BASE_DIR/reads_trimmed}"

# 1. Create results directory
mkdir -p "$OUTPUT_DIR"

# 2. Correctly tag MAG contigs with their original MAG ID
echo "Tagging MAG contigs with their original MAG name..."
rm -f "$OUTPUT_DIR/labeled_mags.fasta"

if [ ! -d "$MAGS_DIR" ]; then
    echo "[ERROR] MAGs directory not found at: $MAGS_DIR"
    exit 1
fi

for mag in "$MAGS_DIR"/*.fa; do
    [ -f "$mag" ] || continue
    mag_name=$(basename "$mag")
    
    # Inject the file name and triple underscore right after the '>'
    sed "s/^>/>${mag_name}___/" "$mag" >> "$OUTPUT_DIR/labeled_mags.fasta"
done

# Verify viral catalog exists
if [ ! -f "$VIRAL_FASTA" ]; then
    echo "[ERROR] Viral catalog not found at: $VIRAL_FASTA"
    exit 1
fi

# 3. Concatenate labeled MAGs and vMAGs into a unified reference file
echo "Combining labeled MAGs and viral catalog into a single reference..."
cat "$OUTPUT_DIR/labeled_mags.fasta" "$VIRAL_FASTA" > "$OUTPUT_DIR/total_community_reference.fna"

# 4. Gather reads matching pattern
R1_FILES=("$READS_DIR"/*_R1.fastq.gz)
R2_FILES=("$READS_DIR"/*_R2.fastq.gz)

if [ ${#R1_FILES[@]} -eq 0 ] || [ ! -e "${R1_FILES[0]}" ]; then
    echo "[ERROR] No FASTQ read files found in $READS_DIR"
    exit 1
fi

R1_STR="${R1_FILES[*]}"
R2_STR="${R2_FILES[*]}"

# 5. Run CoverM once with the unified reference
echo "Running CoverM to normalize all samples together in the same analysis..."
coverm contig \
    --reference "$OUTPUT_DIR/total_community_reference.fna" \
    -1 ${R1_STR} \
    -2 ${R2_STR} \
    -m tpm mean \
    --output-file "$OUTPUT_DIR/total_community_abundance_tpm.tsv"

echo "========================================================"
echo "Finished!!"
echo "========================================================"
