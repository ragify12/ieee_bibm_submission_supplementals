#!/usr/bin/env python3
"""
02_novelty_search.py  --  ITEM #2: are the leads actually novel?

A de novo design paper must show the leads are not accidental recapitulations of
known proteins. This does two independent searches per lead:

  STRUCTURE  Foldseek of the binder chain vs the PDB and the AlphaFold DB,
             reporting the best TM-score. TM-score < ~0.5 to anything natural =
             a genuinely novel fold.
  SEQUENCE   MMseqs2 (or BLAST) of the binder sequence vs UniRef, reporting the
             best percent identity and E-value. No significant hit = novel
             sequence.

Prereqs (all free, CPU-only):
  foldseek            https://github.com/steineggerlab/foldseek
  mmseqs2             https://github.com/soedinglab/MMseqs2
  Foldseek DBs:       foldseek databases PDB pdb tmp ;  foldseek databases Alphafold/UniProt afdb tmp
  MMseqs DB:          mmseqs databases UniRef50 uniref50 tmp
(For BLAST instead of MMseqs: blastp -query binders.fasta -db nr -remote ... )

Usage:
    python 02_novelty_search.py --complex-dir ./lead_complex_pdbs \
        --binder-chain B --pdb-db /db/pdb --afdb /db/afdb --seq-db /db/uniref50
"""
import argparse
import glob
import os
import subprocess
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select, PPBuilder

LEADS = {
    "il23_l94_s482991_mpnn8": "l94",
    "il23_l104_s309115_mpnn6": "l104",
    "il23_l72_s467308_mpnn8": "l72",
    "il23_l92_s675453_mpnn18": "l92",
    "il23_l83_s957191_mpnn2": "l83",
}


class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id


def extract_binder(complex_pdb, chain_id, out_pdb):
    s = PDBParser(QUIET=True).get_structure("x", complex_pdb)
    io = PDBIO()
    io.set_structure(s)
    io.save(out_pdb, ChainSelect(chain_id))
    # also return the sequence for the sequence search
    chain = s[0][chain_id]
    ppb = PPBuilder()
    return "".join(str(pp.get_sequence()) for pp in ppb.build_peptides(chain))


def best_foldseek(query_pdb, db, workdir, tag):
    out = os.path.join(workdir, f"{tag}.m8")
    tmp = os.path.join(workdir, f"tmp_{tag}")
    subprocess.run(
        ["foldseek", "easy-search", query_pdb, db, out, tmp,
         "--format-output", "query,target,evalue,alntmscore,prob",
         "--alignment-type", "1"],
        check=True,
    )
    best = None
    if os.path.exists(out):
        with open(out) as fh:
            for line in fh:
                q, t, ev, tm, prob = line.split("\t")[:5]
                tm = float(tm)
                if best is None or tm > best[1]:
                    best = (t, tm, float(ev))
    return best  # (target, TMscore, evalue) or None


def best_mmseqs(fasta, db, workdir):
    out = os.path.join(workdir, "seq_hits.m8")
    tmp = os.path.join(workdir, "tmp_seq")
    subprocess.run(
        ["mmseqs", "easy-search", fasta, db, out, tmp,
         "--format-output", "query,target,fident,evalue,bits",
         "-s", "7.5"],
        check=True,
    )
    best = {}
    if os.path.exists(out):
        with open(out) as fh:
            for line in fh:
                q, t, fident, ev, bits = line.split("\t")[:5]
                fident = float(fident)
                if q not in best or fident > best[q][1]:
                    best[q] = (t, fident, float(ev))
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complex-dir", required=True,
                    help="dir with the 5 lead complex PDBs")
    ap.add_argument("--binder-chain", default="B",
                    help="chain id of the BINDER in the complex PDBs")
    ap.add_argument("--pdb-db", required=True, help="foldseek PDB database")
    ap.add_argument("--afdb", default=None, help="foldseek AlphaFold database (optional)")
    ap.add_argument("--seq-db", required=True, help="mmseqs sequence database (e.g. uniref50)")
    ap.add_argument("--workdir", default="novelty")
    args = ap.parse_args()
    os.makedirs(args.workdir, exist_ok=True)

    fasta = os.path.join(args.workdir, "binders.fasta")
    rows = []
    binder_pdbs = {}
    with open(fasta, "w") as fh:
        for design, short in LEADS.items():
            # match a file containing the design name
            hits = glob.glob(os.path.join(args.complex_dir, f"*{design}*.pdb"))
            if not hits:
                print(f"WARNING: no complex PDB for {design}; skipping")
                continue
            bpdb = os.path.join(args.workdir, f"{short}_binder.pdb")
            seq = extract_binder(hits[0], args.binder_chain, bpdb)
            binder_pdbs[short] = bpdb
            fh.write(f">{short}\n{seq}\n")

    # sequence search (all leads at once)
    seq_best = best_mmseqs(fasta, args.seq_db, args.workdir)

    # structure search (per lead, vs PDB and optionally AFDB)
    for short, bpdb in binder_pdbs.items():
        pdb_hit = best_foldseek(bpdb, args.pdb_db, args.workdir, f"{short}_pdb")
        afdb_hit = best_foldseek(bpdb, args.afdb, args.workdir, f"{short}_afdb") if args.afdb else None
        struct_candidates = [h for h in (pdb_hit, afdb_hit) if h]
        best_struct = max(struct_candidates, key=lambda h: h[1]) if struct_candidates else None
        s = seq_best.get(short)
        rows.append(dict(
            lead=short,
            best_TMscore=round(best_struct[1], 3) if best_struct else None,
            struct_hit=best_struct[0] if best_struct else None,
            best_seq_identity=round(s[1], 3) if s else None,
            seq_hit=s[0] if s else None,
            seq_evalue=s[2] if s else None,
        ))

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(args.workdir, "novelty_summary.csv"), index=False)
    print(out.to_string(index=False))
    print("\nRead as: TMscore < ~0.5 AND no significant sequence hit => novel de novo design.")


if __name__ == "__main__":
    main()
