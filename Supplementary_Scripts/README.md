# Supplementary Scripts

Analysis scripts for: Paladugu & Alam, "Computational De Novo Design of
Miniprotein Candidate Neutralizers Targeting the IL-23 Cytokine Subunit p19,"
IEEE BIBM 2026.

---

## Run order

| # | Script | Needs | Output |
|---|--------|-------|--------|
| 1 | `00_extract_lead_metrics.py` | BindCraft CSVs only | `../Supplementary_S4_lead_metrics.csv` — hotspot RMSD + design-time ipTM/pTM/iPAE/pLDDT/ΔG/dSASA/SC/HB for the 5 leads |
| 2 | `01_counter_screen.py` | ColabFold + ipSAE + 3DUH | ipSAE of each lead vs p19 / p40 (backs `../Supplementary_S3_p40_selectivity.csv`) |
| 3 | `02_novelty_search.py` | Foldseek + MMseqs2 + DBs | best structural TM-score + sequence identity per lead |
| 4 | `03_iptm_actifptm_ipsae.py` | AF2 outputs + ipSAE + actifpTM | ipSAE vs ipTM vs actifpTM on the same prediction |
| 5 | `leads_vs_p40.py` | ColabFold + ipSAE + p40 structure | per-lead p40 ipSAE/ipTM (fast single-model screen) |
| 6 | `run_p40_zaratan.sbatch` | SLURM + above scripts | submits the full 5-model/3-recycle p40 counter-screen on Zaratan |

Scripts 1 runs anywhere (no GPU). Scripts 2-6 require a cluster.

---

## Confirmed metric values (S4)

`00_extract_lead_metrics.py` reads existing BindCraft statistics — no new
predictions needed. Confirmed output matches `../Supplementary_S4_lead_metrics.csv`:

| Lead | Hotspot RMSD (A) | ipTM | pTM | iPAE | pLDDT | dG (REU) |
|------|-----------------|------|-----|------|-------|----------|
| l94  | 1.69 | 0.86 | 0.88 | 0.18 | 0.93 | -57.5 |
| l104 | 1.50 | 0.90 | 0.90 | 0.17 | 0.96 | -84.9 |
| l72  | 1.37 | 0.91 | 0.85 | 0.20 | 0.95 | -87.3 |
| l92  | 3.17 | 0.92 | 0.85 | 0.26 | 0.90 | -77.0 |
| l83  | 1.86 | 0.88 | 0.83 | 0.26 | 0.92 | -75.0 |

Notable: l94 is best by ipSAE (0.804) but worst by design-time ipTM (0.86) --
direct evidence that ipTM is a poor discriminator at the top of the quality range.

---

## Notes

- ipsae.py calls assume the reference Dunbrack script; swap in your validated
  implementation by editing run_ipsae() / get_ipsae().
- 01_* runs only the p40 arm by default; IL-12 arm needs the p35 (IL12A)
  sequence (UniProt P29459).
- 02_* needs the binder chain ID in your complex PDBs (--binder-chain).
- None of these scripts invent data; they only orchestrate your tools and parse outputs.
