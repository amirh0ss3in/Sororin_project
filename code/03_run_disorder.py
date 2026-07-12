#!/usr/bin/env python3
"""
03_run_disorder.py
==================
Predict, for every residue of Sororin, (a) how intrinsically DISORDERED it is
(IUPred2) and (b) how likely it is to be a DISORDERED BINDING site (ANCHOR2),
then summarise the result region by region.

WHY
---
An "intrinsically disordered" region has no fixed 3D shape -- it wobbles
between many conformations. IUPred2 scores this from the sequence alone: a
value > 0.5 means "predicted disordered". ANCHOR2 additionally flags
disordered stretches that are nonetheless able to bind a partner. The central
claim of the project -- that Sororin is a "fuzzy" hub whose partners bind in
disordered regions -- is tested with exactly these two numbers.

HOW WE RUN IT
-------------
IUPred2A is not a pip package; it is distributed as a small Python library plus
data matrices. This script downloads the reference implementation (the authors'
code, mirrored on GitHub) the first time it runs, then calls it locally. No web
server is involved, so the run is fully reproducible.

Reference: Meszaros, Erdos & Dosztanyi, Nucleic Acids Research 2018;46:W329.

OUTPUTS (written into ../data/)
-------------------------------
* cdca5_iupred.tsv    - per-residue POS / RES / IUPRED2 / ANCHOR2 table
* disorder_summary.txt- mean scores and % disordered per region / motif

Run:  python 03_run_disorder.py
"""

import os
import sys
import urllib.request

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
LIB_DIR = os.path.join(os.path.dirname(__file__), "iupred2a")   # downloaded here
os.makedirs(DATA_DIR, exist_ok=True)

# Files that make up the IUPred2A reference implementation (code + matrices).
GH_BASE = "https://raw.githubusercontent.com/doszilab/iupred2a/master/iupred2a"
LIB_FILES = [
    "iupred2a_lib.py",
    "data/anchor2_energy_matrix",
    "data/anchor2_interface_comp",
    "data/iupred2_long_energy_matrix",
    "data/iupred2_short_energy_matrix",
    "data/long_histogram",
    "data/short_histogram",
]

SEQUENCE = ("MSGRRTRSGGAAQRSGPRAPSPTKPLRRSQRKSGSELPSILPEIWPKTPSAAAVRKPIVL"
            "KRIVAHAVEVPAVQSPRRSPRISFFLEKENEPPGRELTKEDLFKTHSVPATPTSTPVPNP"
            "EAESSSKEGELDARDLEMSKKVRRSYSRLETLGSASTSTPGRRSCFGFEGLLGAEDLSGV"
            "SPVVCSKLTEVPRVCAKPWAPDMTLPGISPPPEKQKRKKKKMPEILKTELDEWAAAMNAE"
            "FEAAEQFDLLVE")

# Regions and interaction determinants to summarise (from UniProt, 1-based).
REGIONS = [
    ("N-terminal IDR",      1,  48),
    ("Central IDR",         72, 142),
    ("C-terminal IDR",      199, 222),
    ("Sororin domain (folded)", 230, 252),
    ("KEN box (APC/C)",     88, 90),
    ("FGF motif (PDS5A/B)", 166, 168),
    ("Ser209 site (ERK2)",  209, 209),
]


def download_library() -> None:
    """Fetch the IUPred2A code + data matrices once, into ./iupred2a/."""
    os.makedirs(os.path.join(LIB_DIR, "data"), exist_ok=True)
    for rel in LIB_FILES:
        dest = os.path.join(LIB_DIR, rel)
        if os.path.exists(dest):
            continue
        urllib.request.urlretrieve(f"{GH_BASE}/{rel}", dest)
    print(f"IUPred2A library ready in {LIB_DIR}")


def main() -> None:
    download_library()
    sys.path.insert(0, LIB_DIR)
    import iupred2a_lib                       # noqa: E402  (import after download)

    # IUPred2 in "long" mode gives per-residue disorder; ANCHOR2 gives the
    # per-residue disordered-binding propensity.
    iupred = iupred2a_lib.iupred(SEQUENCE, "long")[0]
    anchor = iupred2a_lib.anchor2(SEQUENCE)

    # ---- per-residue table ---------------------------------------------
    tsv = os.path.join(DATA_DIR, "cdca5_iupred.tsv")
    with open(tsv, "w") as fh:
        fh.write("POS\tRES\tIUPRED2\tANCHOR2\n")
        for i, aa in enumerate(SEQUENCE):
            fh.write(f"{i+1}\t{aa}\t{iupred[i]:.4f}\t{anchor[i]:.4f}\n")
    print(f"Wrote {tsv}")

    # ---- region summary -------------------------------------------------
    def summarise(a, b):
        vals = iupred[a-1:b]
        frac = sum(v > 0.5 for v in vals) / len(vals)
        mean = sum(vals) / len(vals)
        return mean, frac

    n_dis = sum(v > 0.5 for v in iupred)
    out = os.path.join(DATA_DIR, "disorder_summary.txt")
    with open(out, "w") as fh:
        fh.write(f"Total residues: {len(SEQUENCE)}\n")
        fh.write(f"Disordered (IUPred2 > 0.5): {n_dis} "
                 f"({100*n_dis/len(SEQUENCE):.1f}%)\n")
        fh.write(f"Mean IUPred2: {sum(iupred)/len(iupred):.3f}\n\n")
        fh.write(f"{'region':28s} {'range':>10s} {'meanIUP':>8s} {'%dis':>6s}\n")
        for name, a, b in REGIONS:
            mean, frac = summarise(a, b)
            fh.write(f"{name:28s} {a:4d}-{b:<4d} {mean:8.3f} {100*frac:5.0f}%\n")
    print(f"Wrote {out}")
    print(f"Overall: {100*n_dis/len(SEQUENCE):.1f}% of Sororin is predicted disordered")


if __name__ == "__main__":
    main()
