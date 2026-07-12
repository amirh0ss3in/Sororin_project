#!/usr/bin/env python3
"""
02_string_network.py
====================
Build the Sororin (CDCA5)-centred protein-protein interaction network from the
STRING database and compute its basic topology (nodes, edges, density, degree).

BACKGROUND
----------
STRING is a database of known and predicted protein-protein associations. Each
association has a "combined score" in [0, 1] that measures how confident STRING
is that the two proteins are functionally linked. We keep only reasonably
confident links (score >= 0.4) and treat the result as a graph.

We restrict the graph to the cohesin machinery and its mitotic regulators
(the biologically meaningful neighbourhood of Sororin). Two extra regulators
that are important in the literature but NOT in STRING's high-confidence set --
ERK2 (gene MAPK1) and SPOP -- are added separately by later analysis as
"literature-curated" edges; they are intentionally left out here so that the
STRING statistics stay clean.

OUTPUTS (written into ../data/)
-------------------------------
* string_partners.tsv  - every partner of CDCA5 with its combined score
* string_edges.csv     - the edge list of the cohesin module (for the figure)
* string_stats.txt     - nodes / edges / density / degree ranking

Run:  python 02_string_network.py
"""

import os
import csv
import itertools
import urllib.request
from collections import Counter

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
os.makedirs(DATA_DIR, exist_ok=True)

STRING_API = "https://string-db.org/api"
SPECIES = "9606"          # NCBI taxonomy id for Homo sapiens
SCORE_CUTOFF = 0.4        # discard associations below this combined score

# The cohesin-centred module we analyse (gene symbols understood by STRING).
MODULE = [
    "CDCA5",                                    # Sororin itself
    "SMC3", "SMC1A", "RAD21", "STAG1", "STAG2", # cohesin core ring
    "PDS5A", "PDS5B", "WAPL", "NIPBL", "MAU2",  # cohesin regulators / loader
    "CDK1", "PLK1", "AURKB",                    # mitotic kinases
    "ESPL1", "CDC20",                           # separase, APC/C activator
]


def get_tsv(url: str):
    """Download a STRING TSV endpoint and return it as a list of rows."""
    text = urllib.request.urlopen(url).read().decode()
    return [line.split("\t") for line in text.strip().split("\n")]


def main() -> None:
    # ---- 1. all interaction partners of CDCA5 (for the partners table) ---
    rows = get_tsv(f"{STRING_API}/tsv/interaction_partners"
                   f"?identifiers=CDCA5&species={SPECIES}&limit=30")
    header = rows[0]
    name_i = header.index("preferredName_B")
    score_i = header.index("score")
    partners = sorted(((r[name_i], float(r[score_i])) for r in rows[1:]),
                      key=lambda x: -x[1])
    with open(os.path.join(DATA_DIR, "string_partners.tsv"), "w") as fh:
        fh.write("partner\tcombined_score\n")
        for name, score in partners:
            fh.write(f"{name}\t{score:.3f}\n")
    print(f"CDCA5 has {len(partners)} STRING partners; top = {partners[0]}")

    # ---- 2. the network *within* the chosen module (for the figure) -----
    ids = "%0d".join(MODULE)   # STRING wants identifiers joined by %0d
    rows = get_tsv(f"{STRING_API}/tsv/network"
                   f"?identifiers={ids}&species={SPECIES}")
    header = rows[0]
    a_i, b_i = header.index("preferredName_A"), header.index("preferredName_B")
    s_i = header.index("score")

    edges = {}                 # {(a,b): score}  -- undirected, de-duplicated
    for r in rows[1:]:
        a, b, score = r[a_i], r[b_i], float(r[s_i])
        if score >= SCORE_CUTOFF:
            edges[tuple(sorted((a, b)))] = score

    with open(os.path.join(DATA_DIR, "string_edges.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["a", "b", "score"])
        for (a, b), score in edges.items():
            w.writerow([a, b, f"{score:.3f}"])

    # ---- 3. topology statistics -----------------------------------------
    nodes = set(itertools.chain.from_iterable(edges))
    n, e = len(nodes), len(edges)
    density = 2 * e / (n * (n - 1))            # fraction of possible edges present
    degree = Counter()
    for a, b in edges:
        degree[a] += 1
        degree[b] += 1

    stats_path = os.path.join(DATA_DIR, "string_stats.txt")
    with open(stats_path, "w") as fh:
        fh.write(f"nodes\t{n}\nedges\t{e}\ndensity\t{density:.3f}\n\n")
        fh.write("degree ranking:\n")
        for name, d in degree.most_common():
            fh.write(f"  {name}\t{d}\n")
    print(f"Module: {n} nodes, {e} edges, density = {density:.3f}")
    print(f"Wrote {stats_path}")


if __name__ == "__main__":
    main()
