#!/usr/bin/env bash
# Run the whole Sororin (CDCA5) analysis pipeline from scratch.
# Usage:  bash run_all.sh
set -e                          # stop on the first error
cd "$(dirname "$0")"            # run from the code/ directory

echo "[1/4] Fetching UniProt sequence + features ..."
python3 01_fetch_uniprot.py

echo "[2/4] Building STRING interaction network ..."
python3 02_string_network.py

echo "[3/4] Running IUPred2 + ANCHOR2 disorder prediction ..."
python3 03_run_disorder.py

echo "[4/4] Drawing figures ..."
python3 04_make_figures.py

echo "Done. See ../data/ for tables and ../figures/ for the PDFs."
