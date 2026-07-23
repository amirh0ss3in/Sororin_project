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
| `04_make_figures.py` | the `data/` files above | `figures/Fig1_interaction_network.pdf`, `figures/Fig2_disorder_map.pdf` | Redraws Figures 1-2 as vector PDFs. |
| `05_linear_motifs.py` | ELM API (web), `data/cdca5_iupred.tsv` | `data/elm_motifs_raw.tsv`, `data/elm_motifs_curated.tsv`, `data/elm_motifs_summary.txt` | Predicts short linear motifs (ELM) and links the confident hits to the network partner each one plausibly explains (APC/C, SPOP, ERK2/MAPK1, PLK1, CDK1). |
| `06_make_motif_figure.py` | `data/elm_motifs_curated.tsv` | `figures/Fig3_motif_network_map.pdf` | Draws Figure 3: one row per network partner, showing where its predicted motif sits on the sequence. |
| `07_secondary_structure.py` | `data/cdca5_psipred.ss2`, `data/cdca5_disopred.comb` (PSIPRED web outputs) | `data/secstruct_summary.txt` | Summarises the PSIPRED (secondary structure) and DISOPRED3 (disorder) predictions per region; shows the folded C-terminal domain as the one ordered element and DISOPRED3 disorder (63.9%) matching IUPred2 (61.5%). |
| `08_make_secstruct_figure.py` | `data/cdca5_psipred.ss2`, `data/cdca5_disopred.comb` | `figures/Fig4_secondary_structure.pdf` | Draws Figure 4: PSIPRED secondary-structure cartoon over the DISOPRED3 disorder profile. |
| `09_fuzpred.py` | `data/fuzpred_scores.tsv`, `data/fuzpred_regions.tsv` (FuzPred web outputs) | `data/fuzpred_summary.txt`, `figures/Fig6_fuzpred.pdf` | Summarises the FuzPred context-dependent-binding prediction (pDO / MBM / region classes) and draws Figure 6; quantifies the "fuzzy hub" claim with the Fuxreiter-lab method. |

Run them in order (`01 → 02 → 03 → 04 → 05 → 06 → 07 → 08`); scripts 02 and 03
are independent, 04 needs both, 05 needs 03's disorder profile, 06 needs 05,
and 07–08 read the PSIPRED/DISOPRED3 output files in `data/`.

**Independent network check (Cytoscape).** The same 16-protein module was also
retrieved from STRING and analysed in Cytoscape 3.10 via the stringApp, as an
independent reconstruction using the course's network tool. Cytoscape's network
analyzer returned the same topology as `02_string_network.py` (16 nodes, 107
edges, density 0.892, clustering coefficient 0.925, CDCA5 degree 15). The
stringApp functional-enrichment retrieval (whole-genome background) returned a
STRING PPI enrichment p < 1e-16 and top terms *sister chromatid cohesion*
(GO:0007062, FDR 1.4e-24), KEGG *Cell cycle*, and Reactome sister-chromatid-
cohesion pathways. The exported node table (`data/cytoscape_node_table.csv`),
enrichment table (`data/string_enrichment.csv`) and network image
(`figures/Fig5_cytoscape_network.png`) are included; this step is a GUI
analysis, not part of the scripted pipeline.

**Secondary structure / second disorder predictor (scripts 07–08).** PSIPRED
and DISOPRED3 need a multiple-sequence-alignment step against a large database,
so they are run on the PSIPRED Protein Analysis Workbench
(`http://bioinf.cs.ucl.ac.uk/psipred/`) by submitting `data/cdca5.fasta` with
the PSIPRED 4.0 and DISOPRED3 methods selected. Their per-residue output files
(`cdca5_psipred.ss2`, `cdca5_disopred.comb`, `cdca5_disopred.pbdat`) are stored
in `data/` and read by scripts 07–08 — they are inputs, not recomputed locally.

**Fuzzy-binding prediction (script 09).** FuzPred
(`https://fuzpred.bio.unipd.it`, the Fuxreiter-lab method) was run on accession
Q96FF9; its per-residue scores and region calls (`fuzpred_scores.tsv`,
`fuzpred_regions.tsv`, `fuzpred_features.tsv`) and the AlphaFold cartoon
(`figures/AF-Q96FF9-F1.png`) are stored and read by script 09. Key result:
median pDO 0.40, 68% multimodal-binding residues, 63% context-dependent — i.e.
Sororin binds predominantly through fuzzy, multi-mode regions, with only the
C-terminal domain folding on binding.

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

* **Internet is required** for scripts 01–03 and 05 (they query UniProt,
  STRING, the IUPred2A library on GitHub, and the ELM prediction API
  respectively). Scripts 04 and 06 are fully offline; 07 and 08 read the
  PSIPRED/DISOPRED3 output files already provided in `data/` (the PSIPRED run
  itself is done once on their web server, see below).
* The ELM API (script 05) rate-limits UniProt-ID queries to 1 per 3 minutes;
  the script caches its raw download to `data/elm_motifs_raw.tsv` so re-runs
  don't re-query. Delete that file to force a fresh download.
* The IUPred2A code + energy matrices are downloaded into `code/iupred2a/` on
  the first run of script 03 and reused afterwards. IUPred2A is the reference
  implementation of Mészáros, Erdős & Dosztányi, *Nucleic Acids Research*
  2018;46:W329 — the same method as the iupred2a.elte.hu web server.
* Expected key numbers (a quick sanity check that the run worked):
  the cohesin module has **16 nodes, 107 edges, density ≈ 0.89**; the top
  Sororin partner is **PDS5A (0.999)**; **61.5 %** of the sequence is predicted
  disordered; ELM predicts **14 network-linked motifs**, including the
  annotated **KEN box**, two **SPOP degrons**, a **MAPK docking motif**, and
  a **PLK1 site**; **PSIPRED** calls the protein **78% coil** with its main
  helix at the folded **C-terminal domain (225–245)**, and **DISOPRED3**
  independently gives **63.9% disorder** (vs 61.5% by IUPred2).
* STRING and UniProt are versioned databases. If the providers update their
  data, the exact scores may shift slightly, but the qualitative picture
  (dense cohesin module, majority-disordered Sororin) is stable.
