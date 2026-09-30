# Gallium – Peter Sadler TOSCA spectra

Reproducible processing archive for four gallium complex spectra measured on TOSCA.

## Samples

| Compound | Dataset | Variant | MW (g mol⁻¹) | Mass (mg) |
|---|---|---|---:|---:|
| Mp50-H | Br-H | brominated, protiated | 591.85 | 171 |
| Mp50-D | Br-D | brominated, deuterated | 595.95 | 97 |
| Mp08-H | H | non-brominated, protiated | 452.07 | 83 |
| Mp08-D | D | non-brominated, deuterated | 456.10 | 43 |

The raw empty aluminium cell and short empty-cryostat measurements are also preserved.

## Final processing pipeline

1. **Preserve raw data unchanged.**
2. Apply a **centred 5-point moving average** to the empty aluminium spectrum.
3. Subtract this smoothed aluminium background from each sample.
4. Divide each corrected spectrum by its **sample mass in mg**.
5. In each mass-normalised spectrum, identify the local maximum in **350–430 cm⁻¹**.
6. Normalize each spectrum so this anchor peak has intensity **1.0**.
7. Apply a **centred 3-point moving average** to the normalized sample spectrum.
8. Re-normalize at the **same fixed anchor grid point** after smoothing so the anchor returns to exactly **1.0**.
9. Form deuteration difference spectra as **H − D** separately for the brominated and non-brominated pairs.

### Anchor positions

| Dataset | Compound | Anchor |
|---|---|---:|
| Br-H | Mp50-H | 400.098 cm⁻¹ |
| Br-D | Mp50-D | 400.098 cm⁻¹ |
| H | Mp08-H | 402.098 cm⁻¹ |
| D | Mp08-D | 400.098 cm⁻¹ |

The later anchor normalization removes the absolute mass scale mathematically, but the mass-normalised files are retained because they were an explicit stage of the analysis and remain useful for absolute-intensity comparisons.

## Cryostat data

`raw/emptycryo.dat` is retained. We explored centred **9-point** and **19-point** averages of this short run, stored in `processed/exploratory_cryo/`. No cryostat contribution is subtracted from the final spectra.

## Repository structure

```text
raw/                              original .dat files, unchanged
metadata/                         masses, formulas, anchor positions
processed/
  01_aluminium_5pt/               smoothed empty-Al background
  02_al_subtracted/               sample - smoothed empty Al
  03_mass_normalised/             stage 02 divided by mass
  04_anchor_normalised/           ~400 cm^-1 anchor set to 1
  05_three_point_smoothed/        light smoothing of stage 04
  06_final_renormalised/          final spectra, anchor reset to 1
  07_difference/                  H - D difference spectra
  exploratory_cryo/               9- and 19-point cryostat smooths; not final
plots/                            final comparison and difference PNGs
dashboard/                        static dashboard for GitHub Pages
scripts/process_spectra.py        reproducible processing script
```

## Main comparison plots

![Brominated H vs D](plots/final_brominated_pair.png)

![Non-brominated H vs D](plots/final_non_brominated_pair.png)

## Reproducing the processing

```bash
python -m pip install -r requirements.txt
python scripts/process_spectra.py
```

The script recreates the numerical processed-data stages from the files in `raw/`.

## Dashboard

Open `dashboard/index.html` locally, or publish the `dashboard/` folder with GitHub Pages. A Pages workflow is included in `.github/workflows/pages.yml`.