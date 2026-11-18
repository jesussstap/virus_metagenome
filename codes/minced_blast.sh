#!/bin/bash
# ==============================================================================
# Integrated Host-Virus Interaction Pipeline (Provirus BLAST & CRISPR Spacers)
# ==============================================================================
# Prerequisites: BLAST+, MinCED, and python3 installed

# ==============================================================================
# PATH CONFIGURATION (Modify per your environment or use variables)
# ==============================================================================
BASE_DIR="${BASE_DIR:-.}"
OUTPUT_DIR="${OUTPUT_DIR:-$BASE_DIR/interaction_results}"
MAGS_FASTA="${MAGS_FASTA:-$BASE_DIR/data/mags/total_mags.fasta}"
VIRAL_FASTA="${VIRAL_FASTA:-$BASE_DIR/data/viral_catalog.fasta}"

# Create main results directory
mkdir -p "$OUTPUT_DIR"

echo "========================================================"
echo "Starting Integrated Host-Virus Interaction Analysis"
echo "Base directory: $BASE_DIR"
echo "========================================================"

# ==============================================================================
# PART 1: PROVIRUS DETECTION (Viral Catalog vs. Host MAGs via BLASTn)
# ==============================================================================
echo "-> [1/2] Running Provirus Integration Analysis..."

PROVIRUS_DIR="$OUTPUT_DIR/provirus_blast"
mkdir -p "$PROVIRUS_DIR"

echo "   + Creating BLAST database from host MAGs..."
makeblastdb -in "$MAGS_FASTA" -dbtype nucl -out "$PROVIRUS_DIR/db_mags" > /dev/null

echo "   + Executing BLASTn (Viral Catalog vs. Host MAGs)..."
blastn -query "$VIRAL_FASTA" -db "$PROVIRUS_DIR/db_mags" \
       -out "$PROVIRUS_DIR/provirus_results.tsv" \
       -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore qstart qend" \
       -perc_identity 90 -evalue 1e-5 -num_threads 4

echo "   + Formatting and filtering true proviruses..."
echo -e "qseqid\tsseqid\tpident\tlength\tqlen\tslen\tevalue\tbitscore\tqstart\tqend" > "$PROVIRUS_DIR/header.txt"
cat "$PROVIRUS_DIR/header.txt" "$PROVIRUS_DIR/provirus_results.tsv" > "$PROVIRUS_DIR/provirus_clean.tsv"

# Filter: alignment length > 5000 bp and alignment length >= 2 * query length (or proportional criteria)
awk 'NR==1 || ($4 > 5000 && $6 >= ($4 * 2))' "$PROVIRUS_DIR/provirus_clean.tsv" > "$PROVIRUS_DIR/true_proviruses.tsv"

# Clean up intermediate header file
rm "$PROVIRUS_DIR/header.txt"
echo "BLAST analysis completed. Results in: $PROVIRUS_DIR"

# ==============================================================================
# PART 2: CRISPR SPACER MATCHING (Host Spacers vs. Viral Catalog via MinCED)
# ==============================================================================
echo "------------------------------------------------------------------"
echo "-> [2/2] Running CRISPR Spacer vs. Virus Analysis..."

CRISPR_DIR="$OUTPUT_DIR/crispr_minced"
mkdir -p "$CRISPR_DIR"

echo "   + Extracting CRISPR spacers from host MAGs using MinCED..."
minced -spacers "$MAGS_FASTA" "$CRISPR_DIR/crispr_summary.txt" > /dev/null

echo "   + Creating BLAST database from viral catalog..."
makeblastdb -in "$VIRAL_FASTA" -dbtype nucl -out "$CRISPR_DIR/db_vmags" > /dev/null

echo "   + Crossing spacers against viral database using blastn-short..."
blastn -query "$CRISPR_DIR/crispr_summary_spacers.fa" \
       -db "$CRISPR_DIR/db_vmags" \
       -out "$CRISPR_DIR/crispr_temp.tsv" \
       -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" \
       -task blastn-short \
       -perc_identity 80 \
       -evalue 1e-5 \
       -num_threads 4

echo "   + Formatting CRISPR match results..."
echo -e "qseqid\tsseqid\tpident\tlength\tqlen\tslen\tevalue\tbitscore" > "$CRISPR_DIR/header.txt"
cat "$CRISPR_DIR/header.txt" "$CRISPR_DIR/crispr_temp.tsv" > "$CRISPR_DIR/crispr_matches_final.tsv"

# Clean up temporary files
rm -f "$CRISPR_DIR/header.txt" "$CRISPR_DIR/crispr_temp.tsv"

echo "CRISPR analysis completed. Results in: $CRISPR_DIR"

echo "========================================================"
echo "Finished!"
echo "========================================================"
