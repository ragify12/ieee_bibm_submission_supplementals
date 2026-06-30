# il23-miniprotein-binder-design — Supplementary Data

Supplementary materials for:

> Anonymous, "Computational De Novo Design of Miniprotein Candidate
> Neutralizers Targeting the IL-23 Cytokine Subunit p19,"
> IEEE BIBM 2026.

---

## Contents

| File / Folder | Description |
|---|---|
| `Supplementary_S1_sequences.fasta` | All 50 BindCraft Accepted binder sequences. Five leads marked `LEAD*`. |
| `Supplementary_S2_seed_table.csv` | 42 p19-only designs mapped to their 23 distinct backbone seeds; verifies the "23 independent design seeds" claim in Results §IV-C. |
| `Supplementary_S3_p40_selectivity.csv` | Raw output of the full 5-model / 3-recycle p40 counter-screen (Table III in the paper). All five leads return ipSAE = 0.00. |
| `Supplementary_S4_lead_metrics.csv` | Per-lead BindCraft metrics for the five characterized leads: hotspot RMSD, design-time ipTM/pTM/iPAE, binder pLDDT, ΔG, dSASA, shape complementarity, interface H-bonds. Output of `Supplementary_Scripts/00_extract_lead_metrics.py`. |
| `Supplementary_Scripts/` | Analysis scripts used to produce or verify the results above. See `Supplementary_Scripts/README.md` for run order and prerequisites. |

---

## Reproducing key numbers

- **42 of 50 ≥ 0.60 ipSAE** — recompute from S1 headers.
- **23 independent seeds** — count distinct `seed_id` rows in S2.
- **ipSAE 0.00 on p40** — read S3 directly; produced by `01_counter_screen.py`.
- **Lead biophysics (Table II)** — read S4 directly; produced by `00_extract_lead_metrics.py`.

All results are computational. Raw AlphaFold prediction files and BindCraft
trajectory outputs are available from the authors on request.
