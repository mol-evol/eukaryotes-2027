# eukaryotes-2027

Code and figures for McInerney & Martin, a reanalysis of Bernabeu, Manzano-Morales,
Marcet-Houben & Gabaldón, "Gene ancestries reveal diverse microbial associations during
eukaryogenesis," *Nature* **656**:399–405 (2026).

Every figure and reported number is computed from the authors' deposited data
(Zenodo 13897946) plus standard public reference tables. The rendered figures are in
`figures/`.

## What's in the repo

- Figure scripts: `fig1_machines.py`, `fig2_enrichment.py`, `fig3_confounds.py`,
  `figS1_peak_vs_mode.py`, `figS2_bandwidth.py`, and `excess_residuals.py`.
- Builder scripts (regenerate the processed data from the deposit): `build_master.py`,
  `enrichment_lgt.py`, `phylum_rate.py`, `myxococcota_quality.py`, `fig1a.py`.
- `_common.py`: shared loaders used by the figure scripts.
- `figures/`: the rendered manuscript and SI figures (PNG + SVG).
- `data/reported_modes.tsv`: the authors' reported "wave" modes (small reference table).

The large processed objects (`data/*.pkl`, ~20 MB) and the KEGG KO table
(`data/ko_list.txt`) are **not** committed — the `.pkl` because they are derived and
bulky, the KEGG table because its license does not permit redistribution. Both are
regenerated or downloaded in one step (below).

Requires `python3` with `numpy`, `scipy`, `matplotlib` (and `openpyxl` for the
spreadsheet). Run all scripts from the repo root.

## Reproduce the figures

The figure scripts read `data/*.pkl`. Rebuild those once from the deposit, then run the
scripts:

1. Download the deposit (Zenodo 13897946) and set the input paths at the top of each
   builder (`EX`, `TF`, `COG_DEF`, `DEPOSIT_TREES`, `DEPOSIT`).
2. Build the processed data:
   ```
   python build_master.py       # -> data/master.pkl
   python enrichment_lgt.py      # -> data/lgt_{A,B,C}.pkl
   python phylum_rate.py         # -> data/phylum_rate.pkl, data/rate_enrich_merge.pkl
   ```
3. For Fig 1 only, download the KEGG KO table once:
   ```
   curl 'https://rest.kegg.jp/list/ko' -o data/ko_list.txt
   ```
4. Render:
   ```
   python fig1_machines.py A     # -> figures/Fig1_machines_A.*   (Fig. 1; also B, C for SI)
   python fig2_enrichment.py     # -> figures/Fig2_enrichment.*   (Fig. 2)
   python fig3_confounds.py      # -> figures/Fig3_confounds.*    (Fig. 3)
   python figS1_peak_vs_mode.py  # -> figures/figS1_peak_vs_mode.* (Supp. 1)
   python figS2_bandwidth.py     # -> figures/figS2_bandwidth.*    (Supp. 2)
   python excess_residuals.py    # -> prints the excess-over-confound ratios
   ```

`cog-20.def.tab` (for `build_master.py`) is NCBI's COG-2020 definitions table:
ftp.ncbi.nih.gov/pub/COG/COG2020/data/

## Figure → script

| Figure | Script |
|---|---|
| Fig. 1 | `fig1_machines.py A` (eTOLDB A; B, C → Supp. 3, 4) |
| Fig. 2 | `fig2_enrichment.py` |
| Fig. 3 | `fig3_confounds.py` |
| Supp. 1 | `figS1_peak_vs_mode.py` |
| Supp. 2 | `figS2_bandwidth.py` |
| "positive controls" numbers | `excess_residuals.py` |

## Notes

- The analysis stands on a single database version (eTOLDB A) in the main text; the three
  versions differ only in eukaryotic taxon sampling, so pooling them is pseudo-replication.
  B and C appear in the SI as robustness replicates.
- Fig 1's six mitochondrial machines reproduce the manuscript counts exactly
  (eTOLDB A: Complex I n=26, II n=4, III n=4, IV n=4, mitoribosome large/small 17/7).

## License

MIT — see `LICENSE`.
