# UCST-IDP FRAP Analysis

This repository contains the FRAP analysis workflows used to quantify macromolecular probe transport in GRGDSPYS-based UCST-IDP hydrogels.

The repository includes two analysis workflows:

1. **WT probe-size analysis**: compares 40, 250, and 500 kDa probes in pore and matrix regions of 20 wt% GRGDSPYS hydrogels. This workflow supports Figure S11 and the corresponding Supporting Information table.
2. **Sequence-comparison analysis**: compares GRGDSPYS, GRGNSPWS, and GRGNSPFS hydrogels using 40 and 250 kDa probes in pore and matrix regions. This workflow supports Figure 5 and the corresponding Supporting Information table.

## Repository structure

```text
notebooks/
  01_FRAP_WT_probe_size_colab.ipynb
  02_FRAP_sequence_comparison_colab.ipynb

scripts/
  frap_wt_probe_size.py
  frap_sequence_comparison.py

data/
  README_data.md
  raw/          # optional: raw Excel exports from ZEN
  processed/    # optional: processed metrics/tables

results/
  figures/      # final or exported figures
  tables/       # final or exported tables

docs/
  FRAP_methods_summary.md
```

## Input data format

Each Excel sheet should contain the following columns exported from ZEN or prepared from the ZEN intensity traces:

- `Adjusted time [s]`
- `Intensity Region 1`
- `Intensity Region 2`
- `Intensity Region 3`

Region definitions used in the analysis:

- Region 1 = pore ROI
- Region 2 = matrix / continuous phase ROI
- Region 3 = unbleached reference ROI

Sheet names should include the probe size and replicate number, for example:

- `YS_40kDa_1`
- `YS_250kDa_2`
- `WS_40kDa_3`
- `FS_250kDa_1`

File names are used to infer sequence identity. The scripts recognize sequence tokens such as `GRGDSPYS`, `WT`, `YS`, `GRGNSPWS`, `WS`, `GRGNSPFS`, `FS`, `GRGASPYA`, `AS`, and `YA`.

## FRAP analysis definitions

The analysis uses double normalization:

```text
Fdn_j(t) = (I_j(t)/<I_j,pre>) / (I_ref(t)/<I_ref,pre>)
```

Recovery curves are normalized as:

```text
R_j(t) = (Fdn_j(t) - F0_j)/(Finf_j - F0_j)
```

The recovery half-time is extracted by linear interpolation at:

```text
R_j(t) = 0.5
```

Apparent diffusion coefficients are estimated using the Soumpasis relation:

```text
D_app,j = gamma*w^2/(4*t1/2,j)
```

where `gamma = 0.88` and `w = 7 µm` for a 14 µm diameter bleach ROI.

Mobile fraction is calculated as:

```text
M_j = (Finf_j - F0_j)/(1 - F0_j)
```

Relative diffusivity is calculated as:

```text
K_j = D_app,j/D_PBS
```

`K` is a relative diffusivity, not a hindrance factor that increases with restriction. Lower `K` values indicate greater diffusive restriction relative to free diffusion in PBS.

## PBS diffusion coefficients

The following PBS diffusion coefficients were used consistently in the analysis:

| Probe | D_PBS (µm²/s) |
|---|---:|
| 40 kDa | 9.56 ± 0.62 |
| 250 kDa | 7.54 ± 0.61 |
| 500 kDa | 5.70 ± 1.98 |

## Statistics

Global comparisons were performed using Kruskal-Wallis tests. Pairwise comparisons were performed using Dunn's post hoc test with Holm correction. A corrected p-value < 0.05 was considered statistically significant.

## Running the Colab notebooks

Upload the notebook to Google Colab and run all cells.

- `01_FRAP_WT_probe_size_colab.ipynb` asks for one WT Excel file containing 40, 250, and 500 kDa sheets.
- `02_FRAP_sequence_comparison_colab.ipynb` asks for multiple Excel files for WT, WS, FS, and AS if available.

Each notebook saves the analysis outputs and downloads a zipped output folder.

## Running the scripts locally

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the WT probe-size workflow:

```bash
python scripts/frap_wt_probe_size.py \
  --excel_file data/raw/WT_FRAP_raw.xlsx \
  --output_dir results/wt_probe_size
```

Run the sequence-comparison workflow:

```bash
python scripts/frap_sequence_comparison.py \
  --input_dir data/raw \
  --output_dir results/sequence_comparison
```

Or specify files explicitly:

```bash
python scripts/frap_sequence_comparison.py \
  --excel_files data/raw/WT_FRAP_raw.xlsx data/raw/WS_FRAP_raw.xlsx data/raw/FS_FRAP_raw.xlsx \
  --output_dir results/sequence_comparison
```

Use `--hide_ns` for cleaner publication figures with only significant brackets.

## Outputs

The scripts export:

- per-replicate transport metrics
- all normalized time-point traces
- summary means and standard deviations
- Kruskal-Wallis statistics
- Dunn-Holm pairwise comparisons
- FRAP recovery plots
- boxplots with statistical annotations
