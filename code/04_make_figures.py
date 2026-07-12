#!/usr/bin/env python3
"""
04_make_figures.py
==================
Re-create the two figures used in the report, as vector PDFs:

  Figure 1  -  the Sororin-centred interaction network (from string_edges.csv)
  Figure 2  -  the disorder landscape with interaction sites annotated
               (from cdca5_iupred.tsv)

Run scripts 02 and 03 first so the input data exist. Output goes to ../figures/.

Run:  python 04_make_figures.py
"""

import os
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")                     # render to file, no interactive window
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import networkx as nx

plt.rcParams["pdf.fonttype"] = 42          # editable text in the PDF

HERE = os.path.dirname(__file__)
DATA_DIR = os.path.join(HERE, "..", "data")
FIG_DIR = os.path.join(HERE, "..", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# UniProt-annotated intrinsically disordered regions (1-based).
IDRS = [(1, 48), (72, 142), (199, 222)]


# --------------------------------------------------------------------------- #
#  Figure 2 : disorder landscape                                              #
# --------------------------------------------------------------------------- #
def make_disorder_figure() -> None:
    pos, iu, an = [], [], []
    with open(os.path.join(DATA_DIR, "cdca5_iupred.tsv")) as fh:
        next(fh)
        for line in fh:
            p, _res, i, a = line.split("\t")
            pos.append(int(p)); iu.append(float(i)); an.append(float(a))
    pos, iu = np.array(pos), np.array(iu)

    fig, (axT, ax) = plt.subplots(
        2, 1, figsize=(11, 6.2),
        gridspec_kw={"height_ratios": [1, 4.2]}, sharex=True)

    # ---- top track: domain architecture + phospho-site "lollipops" ----
    axT.set_ylim(0, 1); axT.set_xlim(1, 252); axT.axis("off")
    axT.add_patch(Rectangle((1, 0.42), 251, 0.16, facecolor="#d9d9d9"))
    for a, b in IDRS:                                      # disordered = orange
        axT.add_patch(Rectangle((a, 0.42), b - a, 0.16, facecolor="#f4a582"))
    axT.add_patch(Rectangle((230, 0.40), 22, 0.20,        # folded domain = blue
                            facecolor="#4393c3", edgecolor="k", lw=.5))
    axT.text(241, 0.5, "Sororin\ndomain", ha="center", va="center",
             fontsize=6.5, color="white", fontweight="bold")
    axT.add_patch(Rectangle((88, 0.40), 3, 0.20, facecolor="#1a9850", edgecolor="k", lw=.4))
    axT.add_patch(Rectangle((165, 0.40), 4, 0.20, facecolor="#762a83", edgecolor="k", lw=.4))
    axT.text(89.5, 0.75, "KEN", ha="center", fontsize=6, color="#1a9850", fontweight="bold")
    axT.text(167, 0.22, "FGF", ha="center", fontsize=6, color="#762a83", fontweight="bold")
    for p in [21, 75, 79, 111, 115, 159]:                 # CDK1 proline-directed sites
        axT.plot([p, p], [0.6, 0.9], color="#404040", lw=.8)
        axT.plot(p, 0.9, "o", ms=3.4, color="#404040")
    axT.plot([209, 209], [0.6, 0.95], color="#d6604d", lw=1.1)   # ERK2 site, red
    axT.plot(209, 0.95, "o", ms=5, color="#d6604d")
    axT.text(209, 1.02, "S209\n(ERK2)", ha="center", va="bottom", fontsize=5.8, color="#d6604d")
    axT.text(11, 0.9, "CDK1 sites (S/T-P)", fontsize=6, color="#404040")

    # ---- main panel: IUPred2 and ANCHOR2 curves ----
    ax.axhspan(0, 0.5, color="#f7f7f7")
    for a, b in IDRS:
        ax.axvspan(a, b, color="#fde0d2", alpha=.7)
    ax.plot(pos, iu, color="#b2182b", lw=1.7, label="IUPred2 (disorder)")
    ax.plot(pos, an, color="#2166ac", lw=1.3, alpha=.85, label="ANCHOR2 (disordered binding)")
    ax.axhline(0.5, color="k", ls="--", lw=.7, alpha=.6)
    ax.text(250, 0.52, "disorder threshold 0.5", ha="right", fontsize=6, alpha=.7)
    for x, t, c, y in [(89.5, "APC/C\n(KEN box)", "#1a9850", 0.90),
                       (167, "PDS5A/B\n(FGF SLiM)", "#762a83", 0.30),
                       (209, "ERK2\n(Ser209)", "#d6604d", 0.90),
                       (241, "cohesin core\n(folded domain)", "#2166ac", 0.15)]:
        ax.annotate(t, xy=(x, iu[int(x)-1]), xytext=(x, y), fontsize=6.6,
                    ha="center", color=c, fontweight="bold",
                    arrowprops=dict(arrowstyle="->", color=c, lw=1))
    ax.set_ylim(0, 1.02); ax.set_xlim(1, 252)
    ax.set_xlabel("Residue position", fontsize=9)
    ax.set_ylabel("Prediction score", fontsize=9)
    ax.legend(loc="upper right", fontsize=7.5, framealpha=.9)

    plt.tight_layout(h_pad=0.3)
    out = os.path.join(FIG_DIR, "Fig2_disorder_map.pdf")
    plt.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")


# --------------------------------------------------------------------------- #
#  Figure 1 : interaction network                                             #
# --------------------------------------------------------------------------- #
# functional class of every node (drives the colours)
CLASS = {"CDCA5": "hub",
         "SMC3": "core", "SMC1A": "core", "RAD21": "core", "STAG1": "core", "STAG2": "core",
         "PDS5A": "reg", "PDS5B": "reg", "WAPL": "reg", "NIPBL": "reg", "MAU2": "reg",
         "CDK1": "kin", "PLK1": "kin", "AURKB": "kin", "MAPK1": "kin",
         "ESPL1": "apc", "CDC20": "apc", "SPOP": "ub"}
COLOR = {"hub": "#d6604d", "core": "#4393c3", "reg": "#f4a582",
         "kin": "#1a9850", "apc": "#9970ab", "ub": "#8c510a"}
LABEL = {"hub": "Sororin (CDCA5)", "core": "Cohesin core", "reg": "Cohesin regulators",
         "kin": "Mitotic kinases", "apc": "APC/C / separase", "ub": "SPOP (ubiquitin)"}


def make_network_figure() -> None:
    G = nx.Graph()
    with open(os.path.join(DATA_DIR, "string_edges.csv")) as fh:
        for row in csv.DictReader(fh):
            G.add_edge(row["a"], row["b"], w=float(row["score"]))
    # add the two literature-curated regulatory edges (not in high-conf STRING)
    for a, b in [("CDCA5", "MAPK1"), ("CDCA5", "SPOP")]:
        G.add_edge(a, b, w=0.0, lit=True)

    pos = nx.spring_layout(G, seed=7, k=0.9, iterations=200)
    pos["CDCA5"] = np.array([0.0, 0.0])                    # pin the hub at centre
    deg = dict(G.degree())

    fig, ax = plt.subplots(figsize=(8.6, 7))
    for u, v, d in G.edges(data=True):
        xs, ys = [pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]]
        if d.get("lit"):                                  # dashed red = literature
            ax.plot(xs, ys, color="#d6604d", lw=1.6, ls=(0, (4, 3)), alpha=.9)
        else:                                             # grey, width ~ confidence
            w = d["w"]
            ax.plot(xs, ys, color="#999999", lw=0.4 + 2.2*(w-0.4), alpha=0.35 + 0.4*(w-0.4))
    for n in G.nodes():
        c = COLOR[CLASS[n]]
        size = 650 if n == "CDCA5" else 300 + deg[n]*40
        ax.scatter(*pos[n], s=size, color=c, edgecolors="k", linewidths=.8, zorder=3)
        ax.text(*pos[n], n, fontsize=7.2 if n == "CDCA5" else 6.3, ha="center", va="center",
                zorder=4, fontweight="bold" if n == "CDCA5" else "normal",
                color="white" if CLASS[n] in ("hub", "core", "apc", "ub") else "black")

    legend = [Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR[k],
                     markeredgecolor="k", markersize=9, label=LABEL[k])
              for k in ["hub", "core", "reg", "kin", "apc", "ub"]]
    legend.append(Line2D([0], [0], color="#d6604d", lw=1.6, ls="--",
                         label="literature-curated edge"))
    ax.legend(handles=legend, loc="upper left", fontsize=7, framealpha=.95)
    ax.axis("off")
    ax.set_title("Sororin (CDCA5)-centred interaction network", fontsize=11, pad=8)

    plt.tight_layout()
    out = os.path.join(FIG_DIR, "Fig1_interaction_network.pdf")
    plt.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    make_network_figure()
    make_disorder_figure()
