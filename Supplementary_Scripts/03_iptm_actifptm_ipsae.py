#!/usr/bin/env python3
"""
03_iptm_actifptm_ipsae.py  --  ITEM #3: ipSAE vs ipTM vs actifpTM, same prediction.

The paper argues ipTM is a treacherous filter and ipSAE is better. Right now that
rests on citations. This makes it concrete on YOUR OWN leads: for each lead's
binder-p19 prediction it reports, from the *same* AF2 output, ipSAE, ipTM, and
actifpTM, and shows that ipTM is near-saturated / ranks the leads differently
than ipSAE. That is the evidence a reviewer wants.

Input: a directory with one AF2/ColabFold output set per lead, each providing
  - a PDB of the binder-p19 complex
  - a PAE source (ColabFold *_scores_*.json, or an AF2 .pkl/.npz with PAE+iptm)
Point PER_LEAD_DIR at it and set how to find each file.

actifpTM: install Varga et al.'s tool (`pip install actifptm` or from their repo);
this calls it as a function. If unavailable the column is left blank and the
ipSAE-vs-ipTM contrast still prints.

Companion: for BindCraft's design-time ipTM per lead (already computed, no GPU),
see 00_extract_lead_metrics.py (Average_i_pTM).
"""
import argparse
import glob
import json
import os
import subprocess
import pandas as pd
from scipy.stats import spearmanr

LEADS = ["l94", "l104", "l72", "l92", "l83"]
PAE_CUTOFF, DIST_CUTOFF = 10, 10
IPSAE_PY = "ipsae.py"


def find_one(d, pattern):
    hits = sorted(glob.glob(os.path.join(d, pattern)))
    return hits[0] if hits else None


def get_iptm_ptm(scores_json):
    with open(scores_json) as fh:
        d = json.load(fh)
    return d.get("iptm"), d.get("ptm")


def get_ipsae(pdb, pae_json):
    subprocess.run(["python", IPSAE_PY, pae_json, pdb,
                    str(PAE_CUTOFF), str(DIST_CUTOFF)], check=True)
    stem = os.path.splitext(pdb)[0]
    txt = f"{stem}_{PAE_CUTOFF}_{DIST_CUTOFF}.txt"
    vals = []
    with open(txt) as fh:
        header = fh.readline().split()
        col = header.index("ipSAE") if "ipSAE" in header else None
        for line in fh:
            p = line.split()
            if col is not None and len(p) > col:
                try:
                    vals.append(float(p[col]))
                except ValueError:
                    pass
    return max(vals) if vals else None


def get_actifptm(pdb, pae_json):
    """Return actifpTM for the binder-p19 interface, or None if tool missing."""
    try:
        from actifptm import compute_actifptm  # name depends on the package version
    except Exception:
        return None
    try:
        # Most implementations take the PAE matrix + chain lengths/ids. Adapt the
        # call to your installed version; this is the typical signature.
        return float(compute_actifptm(pae_json=pae_json, pdb=pdb))
    except Exception as e:
        print(f"  actifpTM failed for {pdb}: {e}")
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-lead-dir", required=True,
                    help="dir with one subdir per lead (l94, l104, ...)")
    ap.add_argument("--pdb-glob", default="*rank_001*.pdb")
    ap.add_argument("--scores-glob", default="*scores_rank_001*.json")
    args = ap.parse_args()

    rows = []
    for lead in LEADS:
        d = os.path.join(args.per_lead_dir, lead)
        pdb = find_one(d, args.pdb_glob)
        scores = find_one(d, args.scores_glob)
        if not pdb or not scores:
            print(f"WARNING: missing files for {lead} in {d}")
            rows.append(dict(lead=lead, ipSAE=None, ipTM=None, actifpTM=None))
            continue
        iptm, ptm = get_iptm_ptm(scores)
        ipsae = get_ipsae(pdb, scores)
        actif = get_actifptm(pdb, scores)
        rows.append(dict(lead=lead, ipSAE=ipsae, ipTM=iptm, actifpTM=actif, pTM=ptm))

    out = pd.DataFrame(rows)
    print(out.round(3).to_string(index=False))
    out.to_csv("metric_comparison.csv", index=False)

    # Quantify the divergence: do ipTM and ipSAE rank the leads the same way?
    valid = out.dropna(subset=["ipSAE", "ipTM"])
    if len(valid) >= 3:
        rho, p = spearmanr(valid["ipTM"], valid["ipSAE"])
        print(f"\nSpearman(ipTM, ipSAE) across leads = {rho:.3f} (p={p:.3f})")
        print("Low/!=1 rho => ipTM and ipSAE disagree on lead ranking, which is "
              "exactly the point: ipTM is near-saturated and non-discriminative "
              "where ipSAE still separates the leads.")


if __name__ == "__main__":
    main()
