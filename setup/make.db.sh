#!/usr/bin/env bash

# contam_db_setup.sh
# Downloads curated reference genomes for contamination / host-background screening
# and builds a BLAST database

set -euo pipefail
set -x

### 1. OUTPUT DIRECTORY
DB_DIR="./tick_contam_database"
mkdir -p "$DB_DIR"
cd "$DB_DIR"

### 2. GENOME LIST accession <tab> label
cat <<EOF > accessions.txt
#CF_000001405.39	Human
#F_000001635.27	Mouse
#CF_000146045.2	Yeast
#CF_000002765.6	Plasmodium
#CF_000006985.1	Mycoplasma
#CF_000012865.1	Mycoplasma
#CF_009858895.2	Ecoli
#CF_000007565.1	Pseudomonas
#CF_000008865.1	Bacillus
#CF_000002655.1	Candida
#CF_000149845.1	Aspergillus
#CF_000002335.4	Trypanosoma
#CF_000002875.2	Leishmania
#CF_000006805.1	Halobacterium
#CF_000010525.1	Methanobrevibacter
#CF_000819615.1	Arabidopsis
#CF_001433935.1	Oryza
#CF_000819025.1	PhiX

GCF_016920785.2	Ixodes_scapularis_blacklegged_tick
GCA_037355805.1	Ixodes_ricinus_castor_bean_tick
GCA_030143305.2	Amblyomma_americanum_lone_star_tick
GCF_023375885.2	Dermacentor_andersoni_Rocky_Mountain_wood_tick
GCA_002176555.1	Rhipicephalus_microplus_southern_cattle_tick
GCF_048455015.1	Haemaphysalis_longicornis_longhorned_tick

GCF_000008685.2	Borreliella_burgdorferi_B31_Lyme
GCF_001945665.1	Borreliella_mayonii_MN14_1420_Lyme
GCF_004353425.2	Borrelia_miyamotoi_relapsing_fever
GCF_000018225.1	Rickettsia_rickettsii_Sheila_Smith_RMSF
GCF_000013145.1	Ehrlichia_chaffeensis_Arkansas
GCF_000007765.2	Coxiella_burnetii_RSA493
GCF_000008985.1	Francisella_tularensis_SCHU_S4
GCA_001650135.1	Babesia_microti_Naushon
EOF

### 3. DOWNLOAD ASSEMBLY SUMMARIES

REFSEQ_SUMMARY="assembly_summary_refseq.txt"
GENBANK_SUMMARY="assembly_summary_genbank.txt"

if [ ! -f "$REFSEQ_SUMMARY" ]; then
  echo "📄 Downloading RefSeq assembly summary..."
  wget -q https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt \
    -O "$REFSEQ_SUMMARY"
fi

if [ ! -f "$GENBANK_SUMMARY" ]; then
  echo "📄 Downloading GenBank assembly summary..."
  wget -q https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt \
    -O "$GENBANK_SUMMARY"
fi

### 4. DOWNLOAD FASTAS

while IFS=$'\t' read -r acc label; do

  # Skip blank lines or comments
  [[ -z "${acc:-}" ]] && continue
  [[ "$acc" =~ ^# ]] && continue

  echo "📥 Downloading $acc ($label)"

  # Choose the right summary file based on accession prefix
  if [[ "$acc" == GCF_* ]]; then
    summary_file="$REFSEQ_SUMMARY"
    source_db="RefSeq"
  elif [[ "$acc" == GCA_* ]]; then
    summary_file="$GENBANK_SUMMARY"
    source_db="GenBank"
  else
    echo "⚠️  Skipping $acc: not a GCF_ or GCA_ assembly accession"
    continue
  fi

  # Match exact accession in column 1
  match=$(awk -F '\t' -v acc="$acc" '$1 == acc { print; exit }' "$summary_file" || true)

  if [ -z "$match" ]; then
    echo "❌ Accession $acc not found in $source_db assembly summary"
    continue
  fi

  # Column 20 is ftp_path in NCBI assembly_summary files
  ftp_path=$(echo "$match" | cut -f20)

  if [ -z "$ftp_path" ] || [ "$ftp_path" = "na" ]; then
    echo "❌ No FTP path found for $acc"
    continue
  fi

  file_name=$(basename "$ftp_path")
  fasta_file="${file_name}_genomic.fna.gz"

  echo "   Source: $source_db"
  echo "   FTP:    $ftp_path"
  echo "   File:   $fasta_file"

  wget -q -c "$ftp_path/$fasta_file" -O "$fasta_file" || {
    echo "❌ Failed to download $fasta_file"
    rm -f "$fasta_file"
    continue
  }

done < accessions.txt

### 5. CONCATENATE FASTA FILES

if ! ls *.fna.gz >/dev/null 2>&1; then
  echo "❌ No FASTA files downloaded. Exiting."
  exit 1
fi

cat *.fna.gz > all_contaminants.fna.gz
gunzip -c all_contaminants.fna.gz > all_contaminants.fna

### 6. MAKE BLAST DATABASE

makeblastdb -in all_contaminants.fna \
  -dbtype nucl \
  -parse_seqids \
  -title "CuratedContamDB" \
  -out contaminants_db

### 7. CREATE ACCESSION MAP FILE

awk -F '\t' '
  NF >= 2 && $1 !~ /^#/ && $1 != "" {
    print $1 "\t" $2
  }
' accessions.txt > accession_label_map.txt

### DONE

cat <<EOF

✅ Curated contamination BLAST database created at: $DB_DIR

Database prefix:
  contaminants_db

Combined FASTA:
  all_contaminants.fna

Accession-to-label map:
  accession_label_map.txt

Use with:
  blastn -query sample.fasta -db $DB_DIR/contaminants_db -outfmt 6

EOF
