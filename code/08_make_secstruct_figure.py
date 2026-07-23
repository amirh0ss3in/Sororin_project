#!/usr/bin/env python3
"""
08_make_secstruct_figure.py
===========================
Draw Figure 4: the PSIPRED secondary-structure assignment along Sororin
(helix/strand/coil cartoon) on top of the DISOPRED3 disorder profile, with
the intrinsically disordered regions and the folded C-terminal domain shaded.
The figure shows that the only extended ordered element is the C-terminal
helix (the folded Sororin domain), while the rest of the chain is coil and
predicted disordered.

Run script 07 (and have the PSIPRED/DISOPRED data files in ../data/) first.
Output goes to ../figures/.

Run:  python 08_make_secstruct_figure.py
"""

import os
import numpy as np
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
SS_COLOR = {"H": "#d6604d", "E": "#f4a582", "C": "#bfbfbf"}


def load():
    ss = {}
    with open(os.path.join(DATA_DIR, "cdca5_psipred.ss2")) as fh:
        for line in fh:
            p = line.split()
            if len(p) == 6 and p[0].isdigit():
                ss[int(p[0])] = p[2]
    dpos, dval = [], []
    with open(os.path.join(DATA_DIR, "cdca5_disopred.comb")) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) >= 4 and p[0].isdigit():
                dpos.append(int(p[0])); dval.append(float(p[3]))
    return ss, np.array(dpos), np.array(dval)


def main() -> None:
    ss, dpos, dval = load()
    n = len(ss)

    fig, (axC, ax) = plt.subplots(
        2, 1, figsize=(11, 4.6),
        gridspec_kw={"height_ratios": [1, 4]}, sharex=True)

    # ---- top: PSIPRED secondary-structure cartoon ----
    axC.set_ylim(0, 1); axC.set_xlim(1, n); axC.axis("off")
    for i in range(1, n + 1):
        axC.add_patch(Rectangle((i - 0.5, 0.35), 1, 0.3,
                                facecolor=SS_COLOR[ss[i]], edgecolor="none"))
    axC.text(n + 1, 0.5, "PSIPRED", va="center", fontsize=7, fontweight="bold")

    # ---- bottom: DISOPRED disorder profile ----
    for a, b in IDRS:
        ax.axvspan(a, b, color="#fde0d2", alpha=.55, zorder=0)
    ax.add_patch(Rectangle((230, 0), 22, 1.02, facecolor="#4393c3", alpha=.15, zorder=0))
    ax.plot(dpos, dval, color="#2166ac", lw=1.6, label="DISOPRED3 disorder")
    ax.axhline(0.5, color="k", ls="--", lw=.7, alpha=.6)
    ax.text(n, 0.52, "disorder threshold 0.5", ha="right", fontsize=6.5, alpha=.7)

    ax.annotate("folded Sororin domain\n(C-terminal helix)", xy=(238, 0.08),
                xytext=(205, 0.62), fontsize=6.8, ha="center", color="#1a5276",
                fontweight="bold",
                arrowprops=dict(arrowstyle="->", color="#1a5276", lw=1))
    ax.set_xlim(1, n); ax.set_ylim(0, 1.02)
    ax.set_xlabel("Residue position", fontsize=9)
    ax.set_ylabel("DISOPRED3 disorder", fontsize=9)

    legend = [Line2D([0], [0], color="#2166ac", lw=1.6, label="DISOPRED3 disorder"),
              Rectangle((0, 0), 1, 1, facecolor=SS_COLOR["H"], label="PSIPRED helix"),
              Rectangle((0, 0), 1, 1, facecolor=SS_COLOR["E"], label="PSIPRED strand"),
              Rectangle((0, 0), 1, 1, facecolor=SS_COLOR["C"], label="PSIPRED coil")]
    ax.legend(handles=legend, loc="upper center", ncol=4, fontsize=7,
              framealpha=.9, bbox_to_anchor=(0.5, 1.24))

    axC.set_title("PSIPRED secondary structure over the DISOPRED3 disorder "
                  "profile of Sororin", fontsize=10.5, pad=6)
    plt.tight_layout(h_pad=0.6)
    out = os.path.join(FIG_DIR, "Fig4_secondary_structure.pdf")
    plt.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
