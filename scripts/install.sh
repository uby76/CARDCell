#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STATE="$ROOT/.cardcell"
ENV_ROOT="$STATE/envs"
SOFTWARE_ROOT="$STATE/software"
DB_ROOT="$STATE/databases/card"
LOG_ROOT="$STATE/logs"
RGI_ENV="$ENV_ROOT/rgi"
METAXA_ENV="$ENV_ROOT/metaxa2"
METAXA_DIR="$SOFTWARE_ROOT/Metaxa2_2.2.3"
METAXA_URL="https://microbiology.se/sw/Metaxa2_2.2.3.tar.gz"
METAXA_SHA256="f8f01b6f1a3f9e9968dc9438c84ce8d5a4e189e535dfbf529cbeec2897364c20"

usage() {
  cat <<'EOF'
Usage: bash scripts/install.sh [--skip-card-db] [--force-card-db]

Installs two isolated Conda environments and project-local databases under
.cardcell/. Existing valid components are reused. --force-card-db replaces only
the generated .cardcell/databases/card/localDB directory.
EOF
}

SKIP_CARD=0
FORCE_CARD=0
for arg in "$@"; do
  case "$arg" in
    --skip-card-db) SKIP_CARD=1 ;;
    --force-card-db) FORCE_CARD=1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $arg" >&2; usage >&2; exit 2 ;;
  esac
done

if command -v mamba >/dev/null 2>&1; then
  SOLVER=mamba
elif command -v conda >/dev/null 2>&1; then
  SOLVER=conda
else
  echo "ERROR: Conda or Mamba is required. Install Miniforge/Miniconda, then rerun." >&2
  exit 1
fi
for utility in curl tar awk; do
  if ! command -v "$utility" >/dev/null 2>&1; then
    echo "ERROR: Required system utility not found: $utility" >&2
    exit 1
  fi
done

sha256_file() {
  if command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$1" | awk '{print $1}'
  elif command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    echo "ERROR: shasum or sha256sum is required for download verification." >&2
    return 1
  fi
}

mkdir -p "$ENV_ROOT" "$SOFTWARE_ROOT" "$DB_ROOT" "$LOG_ROOT"

echo "CARDCell installation root: $STATE"
echo "Planned network download is normally below 1 GB in total."
echo "Metaxa2 archive: 47.8 MB; CARD canonical archive and Conda packages vary by platform."
echo "The installer will stop and report before any individual download known to exceed 1 GB."

create_env() {
  local prefix="$1" yaml="$2" label="$3"
  if [[ -x "$prefix/bin/conda" || -d "$prefix/conda-meta" ]]; then
    echo "Reuse $label environment: $prefix"
  else
    echo "Create $label environment: $prefix"
    "$SOLVER" env create --prefix "$prefix" --file "$yaml" --yes
  fi
}

create_env "$RGI_ENV" "$ROOT/envs/rgi.yml" "RGI/KMA"
create_env "$METAXA_ENV" "$ROOT/envs/metaxa2.yml" "Metaxa2 dependency"

if [[ -x "$METAXA_DIR/metaxa2" && -s "$METAXA_DIR/metaxa2_db/SSU/HMMs/B.hmm" && -s "$METAXA_DIR/metaxa2_db/SSU/blast.nhr" ]]; then
  echo "Reuse Metaxa2 2.2.3 and bundled SSU database: $METAXA_DIR"
else
  archive="$SOFTWARE_ROOT/Metaxa2_2.2.3.tar.gz"
  partial="$archive.partial"
  rm -f "$partial"
  echo "Download Metaxa2 2.2.3 (47.8 MB) from the official site"
  curl --fail --location --retry 3 --output "$partial" "$METAXA_URL"
  actual="$(sha256_file "$partial")"
  if [[ "$actual" != "$METAXA_SHA256" ]]; then
    rm -f "$partial"
    echo "ERROR: Metaxa2 SHA-256 mismatch: $actual" >&2
    exit 1
  fi
  mv "$partial" "$archive"
  tar -xzf "$archive" -C "$SOFTWARE_ROOT"
  chmod +x "$METAXA_DIR"/metaxa2*
fi

if [[ "$SKIP_CARD" -eq 0 ]]; then
  if [[ "$FORCE_CARD" -eq 1 && -d "$DB_ROOT/localDB" ]]; then
    stamp="$(date +%Y%m%d-%H%M%S)"
    mv "$DB_ROOT/localDB" "$DB_ROOT/localDB.backup-$stamp"
    echo "Previous CARD localDB preserved as localDB.backup-$stamp"
  fi
  if [[ -s "$DB_ROOT/localDB/card.json" && -s "$DB_ROOT/localDB/card_reference.fasta" ]]; then
    echo "Reuse CARD local database: $DB_ROOT/localDB"
  else
    echo "Download and load the current official CARD canonical database"
    (
      cd "$DB_ROOT"
      export PATH="$RGI_ENV/bin:$PATH"
      "$RGI_ENV/bin/rgi" auto_load --local --canonical --clean
    ) 2>&1 | tee "$LOG_ROOT/card_database_install.log"
  fi
fi

"$SOLVER" list --prefix "$RGI_ENV" --explicit > "$LOG_ROOT/rgi_environment_explicit.txt"
"$SOLVER" list --prefix "$METAXA_ENV" --explicit > "$LOG_ROOT/metaxa2_environment_explicit.txt"
cat > "$LOG_ROOT/metaxa2_source.txt" <<EOF
version=2.2.3
url=$METAXA_URL
sha256=$METAXA_SHA256
includes_database=SSU
EOF

if [[ "$SKIP_CARD" -eq 1 ]]; then
  "$ROOT/scripts/check_installation.sh" --allow-missing-card
else
  "$ROOT/scripts/check_installation.sh"
fi
echo
if [[ "$SKIP_CARD" -eq 1 ]]; then
  echo "Software installation complete; CARD was intentionally skipped."
  echo "Rerun without --skip-card-db before analysis."
else
  echo "Installation complete. Run bin/cardoap_kma16s without dependency path options."
fi
