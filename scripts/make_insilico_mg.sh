#!/usr/bin/env bash
# make_insilico_mg.sh — build a synthetic metagenome from known genomes, so the
# Efesto-vs-FeGenie metagenome benchmark has GROUND TRUTH (you know which genes
# went in). Reads -> assembly -> contigs. Feed the contigs to benchmark_fegenie.py
# with --meta; score against benchmark/truth.tsv (same genome labels).
#
# Two read simulators supported (pick whichever is installed):
#   InSilicoSeq  (iss)      — conda install -c bioconda insilicoseq
#   ART          (art_illumina)
# Assembler: megahit (fast) — conda install -c bioconda megahit
#
# Usage:
#   bash scripts/make_insilico_mg.sh <genomes_dir> <out_dir> [total_reads] [complexity]
#     genomes_dir : dir of *.fna (labelled as in truth.tsv, e.g. from prep_genomes.sh)
#     out_dir     : output dir
#     total_reads : total read pairs (default 3000000 ~ small MG)
#     complexity  : 'even' or 'lognormal' abundance (default lognormal)
#
# Output: <out_dir>/assembly/<out_dir_name>.contigs.fna  -> use as --fna_dir input
set -euo pipefail

GENOMES="${1:?genomes dir}"
OUT="${2:?out dir}"
NREADS="${3:-3000000}"
ABUND="${4:-lognormal}"
THREADS="${THREADS:-8}"
mkdir -p "$OUT/reads" "$OUT/assembly"

# 1. one multi-fasta of all source genomes (records tagged by genome label)
COMBINED="$OUT/source_genomes.fna"
: > "$COMBINED"
for f in "$GENOMES"/*.fna; do
    label="$(basename "$f" .fna)"
    awk -v L="$label" '/^>/{print ">"L"__"substr($0,2); next}{print}' "$f" >> "$COMBINED"
done
echo "[info] combined $(grep -c '^>' "$COMBINED") contigs from $(ls "$GENOMES"/*.fna | wc -l) genomes"

# 2. simulate reads
R1="$OUT/reads/mg_R1.fastq"; R2="$OUT/reads/mg_R2.fastq"
if command -v iss >/dev/null; then
    echo "[info] InSilicoSeq ($ABUND, $NREADS pairs)"
    iss generate --genomes "$COMBINED" --n_reads "$NREADS" \
        --model novaseq --abundance "$ABUND" \
        --cpus "$THREADS" --output "$OUT/reads/mg" >/dev/null
    R1="$OUT/reads/mg_R1.fastq"; R2="$OUT/reads/mg_R2.fastq"
elif command -v art_illumina >/dev/null; then
    echo "[info] ART (even coverage fallback)"
    art_illumina -ss HS25 -i "$COMBINED" -p -l 150 -f 10 -m 350 -s 50 \
        -o "$OUT/reads/mg_R" >/dev/null
    R1="$OUT/reads/mg_R1.fq"; R2="$OUT/reads/mg_R2.fq"
else
    echo "[ERROR] need 'iss' (insilicoseq) or 'art_illumina'"; exit 1
fi

# 3. assemble
command -v megahit >/dev/null || { echo "[ERROR] megahit not found"; exit 1; }
rm -rf "$OUT/assembly/megahit"
megahit -1 "$R1" -2 "$R2" -t "$THREADS" -o "$OUT/assembly/megahit" >/dev/null
cp "$OUT/assembly/megahit/final.contigs.fa" "$OUT/assembly/metagenome.fna"

echo "[done] contigs: $OUT/assembly/metagenome.fna"
echo "Next (metagenome benchmark, both tools, native + shared-ORF):"
echo "  mkdir -p $OUT/mg_in && cp $OUT/assembly/metagenome.fna $OUT/mg_in/"
echo "  python scripts/benchmark_fegenie.py --mode both --meta \\"
echo "         --fna_dir $OUT/mg_in --hmm_dir hmm_library \\"
echo "         --out $OUT/bench --threads $THREADS"
echo
echo "NOTE: ground-truth recall on a co-assembly is approximate — fragmentation and"
echo "      strain collapse lose genes. Report recall trend vs coverage, not absolute."
