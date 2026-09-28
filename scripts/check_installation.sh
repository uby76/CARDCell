#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STATE="$ROOT/.cardcell"
RGI_BIN="$STATE/envs/rgi/bin"
METAXA_BIN="$STATE/envs/metaxa2/bin"
METAXA_DIR="$STATE/software/Metaxa2_2.2.3"
CARD_DB="$STATE/databases/card/localDB"
failed=0
ALLOW_MISSING_CARD=0
if [[ "${1:-}" == "--allow-missing-card" ]]; then
  ALLOW_MISSING_CARD=1
elif [[ $# -gt 0 ]]; then
  echo "Usage: scripts/check_installation.sh [--allow-missing-card]" >&2
  exit 2
fi

check_exe() {
  local label="$1" path="$2"
  if [[ -x "$path" ]]; then printf 'OK\t%s\t%s\n' "$label" "$path"
  else printf 'MISSING\t%s\t%s\n' "$label" "$path"; failed=1; fi
}
check_file() {
  local label="$1" path="$2"
  if [[ -s "$path" ]]; then printf 'OK\t%s\t%s\n' "$label" "$path"
  else printf 'MISSING\t%s\t%s\n' "$label" "$path"; failed=1; fi
}

check_exe "RGI 6.0.8" "$RGI_BIN/rgi"
check_exe "KMA" "$RGI_BIN/kma"
check_exe "samtools" "$RGI_BIN/samtools"
check_exe "bamtools" "$RGI_BIN/bamtools"
check_exe "Perl" "$METAXA_BIN/perl"
check_exe "HMMER hmmscan" "$METAXA_BIN/hmmscan"
check_exe "BLAST blastn" "$METAXA_BIN/blastn"
check_exe "MAFFT" "$METAXA_BIN/mafft"
check_exe "VSEARCH" "$METAXA_BIN/vsearch"
check_exe "Metaxa2 2.2.3" "$METAXA_DIR/metaxa2"
check_file "Metaxa2 bacterial SSU HMM" "$METAXA_DIR/metaxa2_db/SSU/HMMs/B.hmm"
check_file "Metaxa2 SSU BLAST database" "$METAXA_DIR/metaxa2_db/SSU/blast.nhr"
if [[ "$ALLOW_MISSING_CARD" -eq 0 ]]; then
  check_file "CARD card.json" "$CARD_DB/card.json"
  check_file "CARD reference FASTA" "$CARD_DB/card_reference.fasta"
  check_file "CARD database metadata" "$CARD_DB/loaded_databases.json"
else
  echo "SKIPPED CARD database check (--allow-missing-card)"
fi

if [[ "$failed" -ne 0 ]]; then
  echo "Installation check failed. Rerun: bash scripts/install.sh" >&2
  exit 1
fi

if [[ "$ALLOW_MISSING_CARD" -eq 0 ]]; then
  card_version="$("$RGI_BIN/python" -c 'import json,sys; print(json.load(open(sys.argv[1])).get("card_json",{}).get("data_version","unknown"))' "$CARD_DB/loaded_databases.json")"
  echo "CARD data version: $card_version"
  echo "All required software and databases are available."
else
  echo "All required software and the Metaxa2 SSU database are available."
fi
