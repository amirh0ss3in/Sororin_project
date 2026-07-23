#!/usr/bin/env python3
"""
09_fuzpred.py
=============
Summarise the FuzPred prediction of Sororin's context-dependent binding
behaviour, and draw Figure 6 (binding-mode profile + region map + AlphaFold
model). FuzPred (Miskei et al. 2020; Horvath et al. 2020) predicts, per
residue, the probability of disorder-to-order binding (pDO), its complement
disorder-to-disorder / fuzzy binding (pDD), and the multiplicity of binding
modes (MBM); regions are classified as disorder-to-order, disorder-to-disorder
or context-dependent. It is the method of the professor's own laboratory and
directly quantifies the "fuzzy hub" hypothesis.

HOW THE DATA WERE OBTAINED
---------------------------
The UniProt accession Q96FF9 was submitted to the FuzPred server
(https://fuzpred.bio.unipd.it). The per-residue scores and predicted regions
were downloaded and stored in data/ (fuzpred_scores.tsv, fuzpred_regions.tsv,
fuzpred_features.tsv); the AlphaFold cartoon exported by the server is in
figures/AF-Q96FF9-F1.png. As with PSIPRED, this is a web-server run, so these
files are inputs to the script rather than recomputed locally.

OUTPUTS
-------
* data/fuzpred_summary.txt
* figures/Fig6_fuzpred.pdf

Run:  python 09_fuzpred.py
"""

import os
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
try:
    from PIL import Image
    HAVE_PIL = True
except Exception:
    HAVE_PIL = False

HERE = os.path.dirname(__file__)
DATA_DIR = os.path.join(HERE, "..", "data")
FIG_DIR = os.path.join(HERE, "..", "figures")

REGCOL = {"Disorder-to-Order region": "#2166ac",
          "Context-dependent region": "#1a9850",
          "Disorder-to-Disorder region": "#f4a582"}
REGSHORT = {"Disorder-to-Order region": "fold-on-binding",
            "Context-dependent region": "context-dependent (fuzzy)",
            "Disorder-to-Disorder region": "stays disordered (fuzzy)"}
DETS = [("KEN", 88, 90), ("FGF", 166, 168), ("S209", 209, 209),
        ("Sororin\ndomain", 230, 252)]


def load():
    S = {}
    with open(os.path.join(DATA_DIR, "fuzpred_scores.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            S[int(r["position"])] = (float(r["pDO"]), float(r["pDD"]), float(r["MBM"]))
    R = []
    with open(os.path.join(DATA_DIR, "fuzpred_regions.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            R.append((r["type"], int(r["start"]), int(r["end"])))
    return S, R


def summary(S, R):
    n = len(S)
    pdo = np.array([S[i][0] for i in sorted(S)])
    mbm = np.array([S[i][2] for i in sorted(S)])
    lines = []
    lines.append(f"FuzPred prediction for Sororin/CDCA5 (Q96FF9), {n} residues")
    lines.append(f"pDO median {np.median(pdo):.3f}  (below 0.5 => leans disordered/fuzzy binding)")
    lines.append(f"MBM mean {mbm.mean():.3f}  median {np.median(mbm):.3f}")
    lines.append(f"multimodal residues (MBM>0.65): {100*(mbm>0.65).mean():.0f}%")
    lines.append(f"fuzzy-leaning residues (pDD>0.5): {100*(pdo<0.5).mean():.0f}%")
    for t in REGCOL:
        cov = set()
        for tt, a, b in R:
            if tt == t:
                cov.update(range(a, b + 1))
        lines.append(f"{t}: {len(cov)} res ({100*len(cov)/n:.0f}%)")
    lines.append("")
    for name, a, b in [("KEN box", 88, 90), ("FGF motif", 166, 168),
                       ("Ser209", 209, 209), ("Sororin domain", 230, 252)]:
        v = [S[i] for i in range(a, b + 1) if i in S]
        lines.append(f"{name:16s} {a}-{b}: pDO={np.mean([x[0] for x in v]):.2f} "
                     f"pDD={np.mean([x[1] for x in v]):.2f} MBM={np.mean([x[2] for x in v]):.2f}")
    out = os.path.join(DATA_DIR, "fuzpred_summary.txt")
    open(out, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nWrote {out}")


def autocrop(im):
    import numpy as np
    a = np.asarray(im.convert("RGB"))
    mask = (a < 245).any(2)
    if not mask.any():
        return im
    ys, xs = np.where(mask)
    return im.crop((xs.min(), ys.min(), xs.max() + 1, ys.max() + 1))


def figure(S, R):
    pos = np.array(sorted(S)); n = len(pos)
    pdo = np.array([S[i][0] for i in pos]); mbm = np.array([S[i][2] for i in pos])

    have_af = HAVE_PIL and os.path.exists(os.path.join(FIG_DIR, "AF-Q96FF9-F1.png"))
    nrows = 3 if have_af else 2
    hr = [4, 1.7, 3] if have_af else [4, 1.7]
    fig, axes = plt.subplots(nrows, 1, figsize=(11, 6.8 if have_af else 4.6),
                             gridspec_kw={"height_ratios": hr})
    ax, axR = axes[0], axes[1]

    ax.plot(pos, pdo, color="#2166ac", lw=1.6, label="pDO (fold-on-binding)")
    ax.plot(pos, mbm, color="#1a9850", lw=1.3, alpha=.8, label="MBM (binding-mode multiplicity)")
    ax.axhline(0.5, color="#2166ac", ls="--", lw=.6, alpha=.5)
    ax.axhline(0.65, color="#1a9850", ls="--", lw=.6, alpha=.5)
    for name, a, b in DETS:
        x = (a + b) / 2
        ax.axvspan(a, b, color="#d6604d", alpha=.18)
        ax.text(x, 1.02, name, ha="center", va="bottom", fontsize=6.4,
                color="#b2182b", fontweight="bold")
    ax.set_xlim(1, n); ax.set_ylim(0, 1.06)
    ax.set_ylabel("probability", fontsize=9)
    ax.legend(loc="lower right", fontsize=7, framealpha=.9)
    ax.set_title("FuzPred binding-mode analysis of Sororin (Fuxreiter-lab method)",
                 fontsize=10.5, pad=10)
    ax.set_xticklabels([])

    # region track: one row per FuzPred region type (they overlap, so stack them)
    axR.set_xlim(1, n); axR.set_ylim(0, 3); axR.axis("off")
    order = ["Disorder-to-Order region", "Context-dependent region", "Disorder-to-Disorder region"]
    for row, t in enumerate(order):
        y = 2 - row
        for tt, a, b in R:
            if tt == t:
                axR.add_patch(Rectangle((a, y + 0.15), b - a, 0.7, facecolor=REGCOL[t], edgecolor="none"))
        axR.text(0, y + 0.5, REGSHORT[t], ha="right", va="center", fontsize=6.3,
                 color=REGCOL[t], fontweight="bold")

    if have_af:
        axA = axes[2]; axA.axis("off")
        im = Image.open(os.path.join(FIG_DIR, "AF-Q96FF9-F1.png"))
        axA.imshow(autocrop(im))
        axA.set_title("AlphaFold model (blue = fold-on-binding helices; rest disordered)",
                      fontsize=8, pad=2)

    plt.tight_layout(h_pad=0.5)
    out = os.path.join(FIG_DIR, "Fig6_fuzpred.pdf")
    plt.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    S, R = load()
    summary(S, R)
    figure(S, R)
