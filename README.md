# Sororin (CDCA5) systems-biology pipeline — code

These four short Python scripts reproduce **every number and figure** in the
report from scratch, using only public databases (UniProt, STRING) and the
IUPred2A disorder predictor. Nothing is hard-coded from the paper: the scripts
download the data live and recompute the results.

## What each script does

| Script | Input | Output | In one line |
|--------|-------|--------|-------------|
| `01_fetch_uniprot.py` | UniProt Q96FF9 (web) | `data/cdca5.fasta`, `data/uniprot_features.tsv`, `data/uniprot_phossites.tsv` | Downloads the sequence and the annotated regions / motifs / phospho-sites. |
| `02_string_network.py` | STRING API (web) | `data/string_partners.tsv`, `data/string_edges.csv`, `data/string_stats.txt` | Builds the cohesin interaction network and computes nodes/edges/density/degree. |
| `03_run_disorder.py` | downloads IUPred2A (web, first run only) | `data/cdca5_iupred.tsv`, `data/disorder_summary.txt` | Predicts per-residue disorder (IUPred2) and disordered-binding (ANCHOR2). |
| `04_make_figures.py` | the `data/` files above | `figures/Fig1_interaction_network.pdf`, `figures/Fig2_disorder_map.pdf` | Redraws the two figures as vector PDFs. |

Run them in order (`01 → 02 → 03 → 04`); scripts 02 and 03 are independent of
each other, but 04 needs the outputs of both.

## How to run

```bash
# 1. install the dependencies (once)
pip install -r requirements.txt

# 2. run everything
bash run_all.sh
#    ...or run the four scripts by hand, in order:
#    python 01_fetch_uniprot.py
#    python 02_string_network.py
#    python 03_run_disorder.py
#    python 04_make_figures.py
```

Results appear in two sibling folders that the scripts create automatically:

```
Sororin_project/
├── code/            <- these scripts
├── data/            <- generated tables (TSV/CSV/TXT)
└── figures/         <- generated Fig1 and Fig2 (PDF)
```

## Notes / reproducibility

* **Internet is required** for scripts 01–03 (they query UniProt, STRING, and
  download the IUPred2A library from GitHub the first time). Script 04 is
  fully offline.
* The IUPred2A code + energy matrices are downloaded into `code/iupred2a/` on
  the first run of script 03 and reused afterwards. IUPred2A is the reference
  implementation of Mészáros, Erdős & Dosztányi, *Nucleic Acids Research*
  2018;46:W329 — the same method as the iupred2a.elte.hu web server.
* Expected key numbers (a quick sanity check that the run worked):
  the cohesin module has **16 nodes, 107 edges, density ≈ 0.89**; the top
  Sororin partner is **PDS5A (0.999)**; **61.5 %** of the sequence is predicted
  disordered.
* STRING and UniProt are versioned databases. If the providers update their
  data, the exact scores may shift slightly, but the qualitative picture
  (dense cohesin module, majority-disordered Sororin) is stable.
