#!/usr/bin/env bash
# Run the whole Sororin (CDCA5) analysis pipeline from scratch.
# Usage:  bash run_all.sh
set -e                          # stop on the first error
cd "$(dirname "$0")"            # run from the code/ directory

echo "[1/4] Fetching UniProt sequence + features ..."
python3 01_fetch_uniprot.py

echo "[2/4] Building STRING interaction network ..."
python3 02_string_network.py

echo "[3/6] Running IUPred2 + ANCHOR2 disorder prediction ..."
python3 03_run_disorder.py

echo "[4/6] Drawing Figures 1-2 ..."
python3 04_make_figures.py

echo "[5/8] Predicting short linear motifs (ELM) and linking them to the network ..."
python3 05_linear_motifs.py

echo "[6/8] Drawing Figure 3 (motif/network map) ..."
python3 06_make_motif_figure.py

echo "[7/8] Summarising PSIPRED + DISOPRED3 predictions (needs the data/ files) ..."
python3 07_secondary_structure.py

echo "[8/8] Drawing Figure 4 (PSIPRED secondary structure + DISOPRED3 disorder) ..."
python3 08_make_secstruct_figure.py

echo "[9/9] Summarising FuzPred binding-mode prediction + Figure 6 ..."
python3 09_fuzpred.py

echo "Done. See ../data/ for tables and ../figures/ for the PDFs."
echo "Note: scripts 07-09 summarise web-server outputs stored in data/:"
echo "  07-08  PSIPRED/DISOPRED3 (cdca5_psipred.ss2, cdca5_disopred.comb/.pbdat)"
echo "  09     FuzPred (fuzpred_scores.tsv, fuzpred_regions.tsv, fuzpred_features.tsv)"
echo "These were produced on their respective web servers, not recomputed locally."
