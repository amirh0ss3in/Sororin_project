#!/usr/bin/env python3
"""
01_fetch_uniprot.py
===================
Download the reference sequence and the curated feature annotation of human
Sororin (gene CDCA5, UniProt accession Q96FF9) from the UniProt REST API.

WHAT THIS GIVES US
------------------
* the 252-residue amino-acid sequence (saved as FASTA)
* the annotated "features": disordered regions, the C-terminal Sororin domain,
  the KEN box and FGF sequence motifs, and the experimentally observed
  phosphorylation sites (the "modified residue" features).

These are the ground-truth positions that every later script and the paper
refer to (e.g. "FGF motif at residues 166-168", "Ser209 phospho-site").

OUTPUTS (written into ../data/)
-------------------------------
* cdca5.fasta            - the protein sequence
* uniprot_features.tsv   - one row per annotated feature
* uniprot_phossites.tsv  - just the phosphorylation sites, with a flag for
                           whether they are "proline-directed" (S/T-P), i.e.
                           candidate targets of the kinases CDK1 and ERK.

Run:  python 01_fetch_uniprot.py
"""

import os
import json
import urllib.request

ACCESSION = "Q96FF9"                      # human Sororin / CDCA5
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
os.makedirs(DATA_DIR, exist_ok=True)


def fetch(url: str) -> bytes:
    """Small helper that downloads a URL and returns the raw bytes."""
    with urllib.request.urlopen(url) as response:
        return response.read()


def main() -> None:
    # ---- 1. the sequence, in FASTA format -------------------------------
    fasta = fetch(f"https://rest.uniprot.org/uniprotkb/{ACCESSION}.fasta").decode()
    with open(os.path.join(DATA_DIR, "cdca5.fasta"), "w") as fh:
        fh.write(fasta)
    seq = "".join(line for line in fasta.splitlines() if not line.startswith(">"))
    print(f"Sequence length: {len(seq)} residues")

    # ---- 2. the full entry as JSON, so we can read the features ---------
    entry = json.loads(fetch(f"https://rest.uniprot.org/uniprotkb/{ACCESSION}.json"))

    # Write every feature (chain, regions, motifs, modified residues, ...)
    feat_path = os.path.join(DATA_DIR, "uniprot_features.tsv")
    with open(feat_path, "w") as fh:
        fh.write("type\tstart\tend\tdescription\n")
        for f in entry.get("features", []):
            loc = f["location"]
            start = loc["start"].get("value")
            end = loc["end"].get("value")
            fh.write(f"{f['type']}\t{start}\t{end}\t{f.get('description', '')}\n")
    print(f"Wrote {feat_path}")

    # ---- 3. the phosphorylation sites, with the proline-directed flag ---
    # A kinase like CDK1 or ERK phosphorylates a Ser/Thr that is immediately
    # followed by a Proline (the "S/T-P" or "proline-directed" consensus).
    # We flag those, because they are the sites that switch Sororin's
    # interactions on and off during mitosis and in cancer.
    phos_path = os.path.join(DATA_DIR, "uniprot_phossites.tsv")
    with open(phos_path, "w") as fh:
        fh.write("position\tresidue\tnext_residue\tproline_directed\tdescription\n")
        for f in entry.get("features", []):
            if f["type"] != "Modified residue":
                continue
            if "Phospho" not in f.get("description", ""):
                continue
            pos = f["location"]["start"]["value"]           # 1-based position
            res = seq[pos - 1]
            nxt = seq[pos] if pos < len(seq) else "-"        # the +1 residue
            pdir = "yes" if nxt == "P" else "no"
            fh.write(f"{pos}\t{res}\t{nxt}\t{pdir}\t{f.get('description', '')}\n")
    print(f"Wrote {phos_path}")


if __name__ == "__main__":
    main()
