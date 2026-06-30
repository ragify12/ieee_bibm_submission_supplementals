#!/usr/bin/env python3
"""
00_extract_lead_metrics.py  --  ITEM #4 (hotspot RMSD, all 5 leads) and the
quick part of ITEM #3 (BindCraft's own design-time ipTM per lead).

These numbers already exist inside your BindCraft per-design statistics, so no
GPU or re-prediction is needed. This pulls the five leads and prints a table you
can drop straight into Table I (add a "RMSD" column) and into the metric-
comparison discussion.

NOTE on ipTM here: 'Average_i_pTM' is BindCraft's *design-time* AF2-multimer
ipTM (from co-folding). It is NOT the same prediction as your AF2-initial-guess
ipSAE re-score, so use it only to illustrate that ipTM and ipSAE rank the leads
differently. For a strict same-prediction ipTM-vs-ipSAE comparison, use
03_iptm_actifptm_ipsae.py instead.

Usage:
    python 00_extract_lead_metrics.py \
        --p19 final_design_stats_p19.csv \
        --full final_design_stats.csv \
        --out lead_extra_metrics.csv
"""
import argparse
import pandas as pd

# The five leads exactly as named in the BindCraft 'Design' column.
LEADS = [
    "il23_l94_s482991_mpnn8",   # p19-only
    "il23_l104_s309115_mpnn6",  # p19-only
    "il23_l72_s467308_mpnn8",   # p19-only
    "il23_l92_s675453_mpnn18",  # full heterodimer
    "il23_l83_s957191_mpnn2",   # full heterodimer
]

# Columns to surface. Average_* are the per-design averages over the 5 AF2 models.
COLS = [
    "Design",
    "Average_Hotspot_RMSD",   # <-- ITEM #4
    "Average_i_pTM",          # <-- ITEM #3 (BindCraft design-time ipTM)
    "Average_pTM",
    "Average_i_pAE",
    "Average_Binder_pLDDT",
    "Average_dG",
    "Average_dSASA",
    "Average_ShapeComplementarity",
    "Average_n_InterfaceHbonds",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p19", default="final_design_stats_p19.csv")
    ap.add_argument("--full", default="final_design_stats.csv")
    ap.add_argument("--out", default="lead_extra_metrics.csv")
    args = ap.parse_args()

    df = pd.concat(
        [pd.read_csv(args.p19), pd.read_csv(args.full)],
        ignore_index=True,
    )

    have = [c for c in COLS if c in df.columns]
    missing = [c for c in COLS if c not in df.columns]
    if missing:
        print("WARNING: columns not found (check your CSV headers):", missing)

    sub = (
        df[df["Design"].isin(LEADS)][have]
        .set_index("Design")
        .reindex(LEADS)
        .round(3)
    )
    pd.set_option("display.width", 200)
    print(sub.to_string())
    sub.to_csv(args.out)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
