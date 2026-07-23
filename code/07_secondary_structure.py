#!/usr/bin/env python3
"""
07_secondary_structure.py
=========================
Summarise the secondary-structure and disorder predictions produced for
Sororin by the PSIPRED Protein Analysis Workbench (PSIPRED V4.0 for secondary
structure; DISOPRED3 for intrinsic disorder and protein-binding disorder),
and map them onto the annotated regions and the network-linked motifs.

WHY THIS SCRIPT EXISTS
-----------------------
The course teaches secondary-structure prediction as a standard analysis step
(the GOR method and a "metaserver", shown on the ExPASy tools in the p53
worked example), performed alongside disorder prediction. PSIPRED is the
metaserver-grade successor of that step. Running it lets us (i) confirm, with
a dedicated secondary-structure predictor, that the only ordered element of
Sororin is the folded C-terminal domain, and (ii) cross-check the IUPred2
disorder profile of script 03 with a second, independent disorder predictor
(DISOPRED3).

HOW THE DATA WERE OBTAINED
---------------------------
The 252-residue Sororin sequence (data/cdca5.fasta) was submitted to the
PSIPRED Workbench (http://bioinf.cs.ucl.ac.uk/psipred/) with the PSIPRED 4.0
and DISOPRED3 methods selected. Because PSIPRED/DISOPRED require a multiple-
sequence-alignment step against a large database, they run on the server
rather than locally; their per-residue output files are therefore stored in
data/ as inputs to this script (analogous to the UniProt/STRING downloads):

* data/cdca5_psipred.ss2      - PSIPRED per-residue C/H/E call + probabilities
* data/cdca5_disopred.comb    - DISOPRED3 per-residue disorder (* / .) + score
* data/cdca5_disopred.pbdat   - DISOPRED3 protein-binding disorder (^ / - / .)

References: Buchan & Jones, Nucleic Acids Res 2019;47:W402 (PSIPRED
Workbench); Jones & Cozzetto, Bioinformatics 2015;31:857 (DISOPRED3).

OUTPUTS (written into ../data/)
--------------------------------
* secstruct_summary.txt   - overall and per-region secondary-structure /
                            disorder statistics

Run:  python 07_secondary_structure.py
"""

import os

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")

REGIONS = [
    ("N-terminal IDR",          1,  48),
    ("Central IDR",             72, 142),
    ("C-terminal IDR",          199, 222),
    ("Sororin domain (folded)", 230, 252),
    ("KEN box (APC/C)",         88, 90),
    ("FGF motif (PDS5A/B)",     166, 168),
    ("MAPK docking (ERK2)",     77, 84),
    ("SPOP degron 1",           122, 126),
    ("SPOP degron 2",           155, 159),
    ("PLK1 motif",              201, 207),
]


def load_psipred():
    ss = {}
    with open(os.path.join(DATA_DIR, "cdca5_psipred.ss2")) as fh:
        for line in fh:
            p = line.split()
            if len(p) == 6 and p[0].isdigit():
                ss[int(p[0])] = p[2]          # C / H / E
    return ss


def load_disopred():
    diso = {}
    with open(os.path.join(DATA_DIR, "cdca5_disopred.comb")) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) >= 4 and p[0].isdigit():
                diso[int(p[0])] = (p[2] == "*")   # True if disordered
    return diso


def segments(ss, ch):
    segs, cur = [], None
    for i in sorted(ss):
        if ss[i] == ch:
            cur = [i, i] if cur is None else [cur[0], i]
        elif cur:
            segs.append(tuple(cur)); cur = None
    if cur:
        segs.append(tuple(cur))
    return segs


def main() -> None:
    ss = load_psipred()
    diso = load_disopred()
    n = len(ss)

    def pct(rng, pred):
        rng = list(rng)
        return 100 * sum(pred(i) for i in rng) / len(rng)

    out = os.path.join(DATA_DIR, "secstruct_summary.txt")
    with open(out, "w") as fh:
        fh.write("PSIPRED (secondary structure) + DISOPRED3 (disorder) for Sororin/CDCA5\n")
        fh.write(f"Residues: {n}\n")
        fh.write(f"PSIPRED overall: helix {pct(ss,lambda i: ss[i]=='H'):.1f}%  "
                 f"strand {pct(ss,lambda i: ss[i]=='E'):.1f}%  "
                 f"coil {pct(ss,lambda i: ss[i]=='C'):.1f}%\n")
        fh.write(f"DISOPRED3 disordered: {pct(diso,lambda i: diso[i]):.1f}% "
                 f"(IUPred2 in script 03 gave 61.5%)\n")
        fh.write(f"PSIPRED helix segments (>=4 res): "
                 f"{[s for s in segments(ss,'H') if s[1]-s[0]+1>=4]}\n")
        fh.write(f"PSIPRED strand segments: {segments(ss,'E')}\n\n")
        fh.write(f"{'region':26s}{'range':>10s}{'%helix':>8s}{'%strand':>8s}"
                 f"{'%coil':>7s}{'%disord':>9s}\n")
        for name, a, b in REGIONS:
            rng = range(a, b + 1)
            fh.write(f"{name:26s}{f'{a}-{b}':>10s}"
                     f"{pct(rng,lambda i: ss[i]=='H'):7.0f}%"
                     f"{pct(rng,lambda i: ss[i]=='E'):7.0f}%"
                     f"{pct(rng,lambda i: ss[i]=='C'):6.0f}%"
                     f"{pct(rng,lambda i: diso[i]):8.0f}%\n")
    print(f"Wrote {out}")
    print(open(out).read())


if __name__ == "__main__":
    main()
