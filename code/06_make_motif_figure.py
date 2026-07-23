#!/usr/bin/env python3
"""
06_make_motif_figure.py
========================
Draw Figure 3: a schematic that maps every curated, network-linked ELM motif
(from script 05) directly onto the Sororin sequence, one row per network
partner, coloured to match the partner's node colour in Figure 1. This is the
figure that makes the "linear motifs are what connect the network to the
sequence" argument visually, complementing Figure 2 (which makes the
disorder argument).

Run scripts 02, 03 and 05 first. Output goes to ../figures/.

Run:  python 06_make_motif_figure.py
"""

import os
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

plt.rcParams["pdf.fonttype"] = 42

HERE = os.path.dirname(__file__)
DATA_DIR = os.path.join(HERE, "..", "data")
FIG_DIR = os.path.join(HERE, "..", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

IDRS = [(1, 48), (72, 142), (199, 222)]

# Same colour scheme as 04_make_figures.py, so a partner's colour matches its
# node colour in Figure 1.
COLOR = {"CDC20 / APC-C": "#9970ab", "SPOP": "#8c510a",
         "MAPK1 (ERK2)": "#1a9850", "PLK1": "#1a9850", "CDK1": "#1a9850",
         "PDS5A/B (FGF, from script 01)": "#f4a582"}

# Row order (top to bottom) and the label shown on the left.
ROWS = ["CDC20 / APC-C", "SPOP", "MAPK1 (ERK2)", "PLK1", "CDK1",
        "PDS5A/B (FGF, from script 01)"]


def load_curated():
    hits = {r: [] for r in ROWS}
    with open(os.path.join(DATA_DIR, "elm_motifs_curated.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            hits[row["network_partner"]].append(
                (int(row["start"]), int(row["stop"]), row["is_annotated"] == "True"))
    hits["PDS5A/B (FGF, from script 01)"].append((166, 168, False))  # UniProt motif
    return hits


def main() -> None:
    hits = load_curated()

    fig, ax = plt.subplots(figsize=(11, 3.6))
    ax.set_xlim(1, 252)
    ax.set_ylim(0, len(ROWS) + 0.6)

    for a, b in IDRS:                                   # disordered background
        ax.axvspan(a, b, color="#fde0d2", alpha=.55, zorder=0)
    ax.add_patch(Rectangle((230, 0), 22, len(ROWS) + 0.6,
                            facecolor="#4393c3", alpha=.18, zorder=0))

    for i, row in enumerate(reversed(ROWS)):
        y = i + 1
        ax.axhline(y, color="#e0e0e0", lw=.6, zorder=1)
        for a, b, annotated in hits[row]:
            w = max(b - a, 2)
            ax.add_patch(Rectangle((a, y - 0.28), w, 0.56,
                                    facecolor=COLOR[row], edgecolor="k",
                                    lw=1.1 if annotated else 0.4, zorder=2))
        ax.text(-6, y, row, ha="right", va="center", fontsize=8.5, fontweight="bold")

    ax.set_yticks([])
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.set_xlabel("Residue position", fontsize=9)
    ax.set_title("Predicted short linear motifs (ELM) link specific network "
                  "partners to the Sororin sequence", fontsize=10.5, pad=10)

    legend = [
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#bbbbbb",
               markeredgecolor="k", markeredgewidth=1.1, markersize=10,
               label="ELM-annotated instance (experimentally known)"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#bbbbbb",
               markeredgecolor="k", markeredgewidth=.4, markersize=10,
               label="ELM-predicted, not yet experimentally annotated"),
        Rectangle((0, 0), 1, 1, facecolor="#fde0d2", alpha=.55,
                  label="intrinsically disordered region"),
        Rectangle((0, 0), 1, 1, facecolor="#4393c3", alpha=.18,
                  label="folded Sororin domain"),
    ]
    ax.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, -0.28),
              ncol=2, fontsize=7.3, frameon=False)

    plt.tight_layout()
    out = os.path.join(FIG_DIR, "Fig3_motif_network_map.pdf")
    plt.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
