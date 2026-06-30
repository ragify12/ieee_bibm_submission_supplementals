#!/usr/bin/env python3
"""
01_counter_screen.py  --  ITEM #1: selectivity counter-screen.

For every lead, predict the binder against (i) p19 [positive control],
(ii) p40 alone, and (iii) intact IL-12 (p35 + p40), then score every prediction
with the SAME ipSAE pipeline you use in the paper. Expected result if your
selectivity claim holds: high ipSAE on p19, near-zero ipSAE on p40 and IL-12.
That turns "selectivity is structurally inferred" into "selectivity is measured."

Pipeline per (lead, target):
    sequences -> ColabFold-multimer prediction -> PDB + PAE json
              -> ipSAE (your scorer)            -> table row (ipSAE, ipTM, minPAE)

Why ColabFold-multimer for the off-targets instead of AF2-initial-guess:
initial-guess needs a designed pose, and there is no designed binder<->p40 pose.
The clean specificity test is "can AF2 find ANY confident interface at all": if
it cannot from a multimer prediction, that is strong evidence of non-binding.
For the p19 positive control you can instead point run_structure_prediction()
at your AF2-initial-guess wrapper so the control matches the main text exactly.

Dependencies: ColabFold (colabfold_batch on PATH), Biopython, and either the
reference ipsae.py OR your own validated ipSAE script (see run_ipsae()).
"""
import argparse
import json
import os
import subprocess
import glob
import pandas as pd
from Bio import SeqIO
from Bio.PDB import PDBParser, PPBuilder

# ----------------------------------------------------------------------------
# CONFIG  -- fill these in for your environment
# ----------------------------------------------------------------------------
LEADS = [
    "il23_l94_s482991_mpnn8",
    "il23_l104_s309115_mpnn6",
    "il23_l72_s467308_mpnn8",
    "il23_l92_s675453_mpnn18",
    "il23_l83_s957191_mpnn2",
]

# Where to read each lead's binder sequence from. The BindCraft stats CSVs carry
# a 'Sequence' column; we read it from there so nothing is retyped by hand.
P19_CSV = "final_design_stats_p19.csv"
FULL_CSV = "final_design_stats.csv"
SEQ_COLUMN = "Sequence"          # <-- confirm this is the binder-sequence header

# 3DUH gives you p40 (chain A) and p19 (chain C) directly.
PDB_3DUH = "3duh.pdb"
P40_CHAIN = "A"
P19_CHAIN = "C"

# IL-12 p35 (IL12A) is NOT in 3DUH. Provide it yourself: UniProt P29459, mature
# chain (drop the signal peptide), or pull it from an IL-12 structure
# (e.g. PDB 1F45 / 3HMX). Paste the one-letter sequence here, or set a FASTA path.
P35_IL12A_SEQ = ""               # <-- REQUIRED for the IL-12 arm. e.g. "RNLPVAT..."
# Leave P35_IL12A_SEQ empty to skip the IL-12 arm and run only the p40 arm.

# ipSAE cutoffs (match whatever you used in the paper; reference defaults shown)
PAE_CUTOFF = 10
DIST_CUTOFF = 10
IPSAE_PY = "ipsae.py"            # path to reference ipsae.py (or your script)

OUTDIR = "counter_screen"
# ----------------------------------------------------------------------------


def chain_sequence(pdb_path, chain_id):
    """One-letter sequence of a chain from a PDB (modeled residues only)."""
    structure = PDBParser(QUIET=True).get_structure("x", pdb_path)
    ppb = PPBuilder()
    chain = structure[0][chain_id]
    return "".join(str(pp.get_sequence()) for pp in ppb.build_peptides(chain))


def load_lead_sequences():
    df = pd.concat([pd.read_csv(P19_CSV), pd.read_csv(FULL_CSV)], ignore_index=True)
    if SEQ_COLUMN not in df.columns:
        raise SystemExit(
            f"'{SEQ_COLUMN}' not in CSV. Available: "
            + ", ".join(c for c in df.columns if 'seq' in c.lower())
        )
    df = df[df["Design"].isin(LEADS)].set_index("Design")
    return {name: df.loc[name, SEQ_COLUMN] for name in LEADS}


def run_structure_prediction(jobname, sequences, outdir):
    """Predict a complex with ColabFold-multimer.

    `sequences` is a list of chains; ColabFold joins them with ':' for a complex.
    Returns (pdb_path, scores_json_path) for the top-ranked model.

    >>> To use AF2-initial-guess for the p19 positive control instead, replace
        the body of this function with a call to your initial-guess wrapper that
        emits a PDB + a PAE json/npz, and return those two paths. <<<
    """
    os.makedirs(outdir, exist_ok=True)
    fasta = os.path.join(outdir, f"{jobname}.fasta")
    with open(fasta, "w") as fh:
        fh.write(f">{jobname}\n{':'.join(sequences)}\n")

    subprocess.run(
        ["colabfold_batch", "--num-models", "5", "--num-recycle", "3",
         fasta, outdir],
        check=True,
    )

    pdbs = sorted(glob.glob(os.path.join(outdir, f"{jobname}*rank_001*.pdb")))
    jsons = sorted(glob.glob(os.path.join(outdir, f"{jobname}*scores_rank_001*.json")))
    if not pdbs or not jsons:
        raise RuntimeError(f"ColabFold output not found for {jobname} in {outdir}")
    return pdbs[0], jsons[0]


def run_ipsae(pdb_path, pae_json):
    """Call ipSAE and return the headline (max over chain directions) value.

    Replace with your own validated ipSAE script if you prefer; just return the
    single max ipSAE for the binder<->target interface.
    """
    subprocess.run(
        ["python", IPSAE_PY, pae_json, pdb_path, str(PAE_CUTOFF), str(DIST_CUTOFF)],
        check=True,
    )
    # ipsae.py writes <stem>_<pae>_<dist>.txt with an 'ipSAE' column; take the max.
    stem = os.path.splitext(pdb_path)[0]
    out_txt = f"{stem}_{PAE_CUTOFF}_{DIST_CUTOFF}.txt"
    best = float("nan")
    with open(out_txt) as fh:
        header = fh.readline().split()
        col = header.index("ipSAE") if "ipSAE" in header else None
        vals = []
        for line in fh:
            parts = line.split()
            if col is not None and len(parts) > col:
                try:
                    vals.append(float(parts[col]))
                except ValueError:
                    pass
        if vals:
            best = max(vals)
    return best


def parse_iptm(scores_json):
    with open(scores_json) as fh:
        d = json.load(fh)
    # ColabFold scores json carries 'iptm' and 'ptm'; min interface PAE is the
    # smallest off-diagonal PAE block value.
    return d.get("iptm"), d.get("ptm")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=OUTDIR)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    binders = load_lead_sequences()
    p40 = chain_sequence(PDB_3DUH, P40_CHAIN)
    p19 = chain_sequence(PDB_3DUH, P19_CHAIN)

    targets = {
        "p19": [p19],          # positive control
        "p40": [p40],          # off-target
    }
    if P35_IL12A_SEQ.strip():
        targets["IL12"] = [P35_IL12A_SEQ.strip(), p40]   # intact IL-12 = p35 + p40
    else:
        print("NOTE: P35_IL12A_SEQ empty -> skipping IL-12 arm (p40 arm still runs).")

    rows = []
    for lead, bseq in binders.items():
        for tname, tchains in targets.items():
            job = f"{lead}__vs__{tname}"
            print(f"=== {job} ===")
            try:
                pdb, scores = run_structure_prediction(job, [bseq] + tchains, args.outdir)
                ipsae = run_ipsae(pdb, scores)
                iptm, ptm = parse_iptm(scores)
                rows.append(dict(lead=lead, target=tname, ipSAE=ipsae,
                                 ipTM=iptm, pTM=ptm))
            except Exception as e:
                print(f"  FAILED {job}: {e}")
                rows.append(dict(lead=lead, target=tname, ipSAE=None,
                                 ipTM=None, pTM=None))

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(args.outdir, "counter_screen_ipsae.csv"), index=False)
    print("\n", out.to_string(index=False))
    print("\nInterpretation: p19 ipSAE high, p40/IL12 ipSAE << 0.60 => selective.")


if __name__ == "__main__":
    main()
