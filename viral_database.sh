#!/bin/bash
# ========================================================
# Viral Contigs Clustering & Dereplication Pipeline
# ========================================================
# Prerequisites: Biopython, BLAST+, CheckV, and python3 installed

echo "========================================================"
echo "Starting Smart Viral Contigs Integration (Global Catalog)"
echo "========================================================"

# ==========================================
# 1. PATH CONFIGURATION (Modify per your environment)
# ==========================================
BASE_DIR="${BASE_DIR:-.}"
ORIGINAL_CONTIGS="${ORIGINAL_CONTIGS:-$BASE_DIR/data/assemblies}"
VIBRANT_DIR="${VIBRANT_DIR:-$BASE_DIR/data/vibrant}"
VS2_DIR="${VS2_DIR:-$BASE_DIR/data/virsorter2}"
GENOMAD_DIR="${GENOMAD_DIR:-$BASE_DIR/data/genomad}"

OUT_DIR="${OUT_DIR:-$BASE_DIR/viral_database}"

mkdir -p "$OUT_DIR/clustering"
mkdir -p "$OUT_DIR/checkv_final"

RAW_FILE="$OUT_DIR/all_raw_viruses.fasta"
> "$RAW_FILE" # Clear file in case of re-run

# ==========================================
# 2. PYTHON SCRIPT: FASTA EXTRACTOR
# ==========================================
cat << 'EOF' > "$OUT_DIR/extractor.py"
import sys
from Bio import SeqIO
rep_file, fasta_in, fasta_out = sys.argv[1], sys.argv[2], sys.argv[3]
reps = set(line.strip() for line in open(rep_file))
records = [rec for rec in SeqIO.parse(fasta_in, "fasta") if rec.id.split()[0] in reps]
SeqIO.write(records, fasta_out, "fasta")
EOF

# ==========================================
# 3. MASSIVE COLLECTION & SAFE RENAMING
# ==========================================
echo "-> Collecting all predictions and protecting contig names..."

for contig_file in "$ORIGINAL_CONTIGS"/*.fa; do
    
    [ -f "$contig_file" ] || continue
    sample_full=$(basename "$contig_file" .fa)
    sample=$(echo "$sample_full" | sed 's/\.Assembly//g')
    
    echo "   Collecting: $sample"
    
    vib_fa=$(find "$VIBRANT_DIR/$sample" -name "*.phages_combined.fna" 2>/dev/null | head -n 1)
    vs2_fa="$VS2_DIR/${sample_full}.fa.fa/final-viral-combined.fa"
    genomad_fa="$GENOMAD_DIR/${sample_full}.fasta"
    
    # Extract and rename VIBRANT
    if [[ -n "$vib_fa" ]]; then
        awk -v prefix="${sample}_vib_" '/^>/{print ">"prefix substr($0,2); next} {print}' "$vib_fa" >> "$RAW_FILE"
    fi
    
    # Extract and rename VirSorter2
    if [[ -f "$vs2_fa" ]]; then
        awk -v prefix="${sample}_vs2_" '/^>/{print ">"prefix substr($0,2); next} {print}' "$vs2_fa" >> "$RAW_FILE"
    fi
    
    # Extract and rename geNomad
    if [[ -f "$genomad_fa" ]]; then
        awk -v prefix="${sample}_gn_" '/^>/{print ">"prefix substr($0,2); next} {print}' "$genomad_fa" >> "$RAW_FILE"
    fi
done

total_seqs=$(grep -c "^>" "$RAW_FILE")
echo "-> Collected $total_seqs viral sequences in total."

# ==========================================
# 4. DEREPLICATION (Remove duplicates / 95% ANI)
# ==========================================
echo "------------------------------------------------"
echo "-> Starting Global Clustering and Dereplication (95% ANI)..."

db_blast="$OUT_DIR/clustering/blast_db"
out_blast="$OUT_DIR/clustering/all_blast.tsv"
out_ani="$OUT_DIR/clustering/all_ani.tsv"
out_clusters="$OUT_DIR/clustering/all_clusters.tsv"
lista_reps="$OUT_DIR/clustering/unique_representatives.txt"

UNIQUE_FASTA="$OUT_DIR/viral_catalog.fasta"
CHECKV_SCRIPTS_DIR="${CHECKV_SCRIPTS_DIR:-$HOME/checkv/scripts}"

# All-vs-all BLAST
makeblastdb -in "$RAW_FILE" -dbtype nucl -out "$db_blast" > /dev/null
blastn -query "$RAW_FILE" -db "$db_blast" -outfmt '6 std qlen slen' -max_target_seqs 10000 -out "$out_blast" -num_threads 8

# Clustering with CheckV scripts
(
    cd "$CHECKV_SCRIPTS_DIR"
    python3 anicalc.py -i "$out_blast" -o "$out_ani"
    python3 aniclust.py --fna "$RAW_FILE" --ani "$out_ani" --out "$out_clusters" --min_ani 95 --min_tcov 85 --min_qcov 0
)

# Extract a single representative sequence per cluster
awk -F'\t' '{print $1}' "$out_clusters" | sort -u > "$lista_reps"
python3 "$OUT_DIR/extractor.py" "$lista_reps" "$RAW_FILE" "$UNIQUE_FASTA"

unique_total=$(grep -c "^>" "$UNIQUE_FASTA")
echo "-> Dereplication complete. Reduced from $total_seqs to $unique_total unique viral sequences."

# ==========================================
# 5. QUALITY CONTROL WITH CHECKV
# ==========================================
echo "------------------------------------------------"
echo "-> Evaluating viral database quality using CheckV..."

CHECKV_DB="${CHECKV_DB:-$HOME/checkv/checkv-db-v1.5}"

checkv end_to_end \
    "$UNIQUE_FASTA" \
    "$OUT_DIR/checkv_final" \
    -d "$CHECKV_DB" \
    -t 8 

echo "========================================================"
echo "VIRAL DATABASE READY"
echo "========================================================"
