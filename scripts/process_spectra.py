#!/usr/bin/env python3
"""Reproduce the Gallium / Peter Sadler TOSCA processing pipeline.

Final pipeline
--------------
1. Centred 5-point arithmetic mean of the empty aluminium spectrum.
2. Subtract smoothed aluminium from each sample.
3. Divide each corrected spectrum by sample mass (mg).
4. Identify the local maximum in 350–430 cm^-1.
5. Normalize that anchor peak to 1.
6. Apply a centred 3-point arithmetic mean.
7. Re-normalize at the same fixed anchor grid point after smoothing.
8. Calculate H-D differences separately for brominated/non-brominated pairs.

The short empty-cryostat run is preserved and exploratory 9/19-point smooths
are exported, but no cryostat subtraction is used in the final spectra.
"""
from pathlib import Path
import json
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw"
PROCESSED = ROOT / "processed"
META = ROOT / "metadata"
PLOTS = ROOT / "plots"
ASSETS = ROOT / "dashboard" / "assets"

SAMPLES = {
    "Br-H": {"file": "Br-H.dat", "mass_mg": 171.0, "compound": "Mp50-H", "mw": 591.85,
             "formula": "C16 H16 Br2 Ga N2 O4,N O3", "series": "brominated", "isotope": "H"},
    "Br-D": {"file": "Br-D.dat", "mass_mg": 97.0, "compound": "Mp50-D", "mw": 595.95,
             "formula": "C16 H12 D4 Br2 Ga N2 O4,N O3", "series": "brominated", "isotope": "D"},
    "H": {"file": "H.dat", "mass_mg": 83.0, "compound": "Mp08-H", "mw": 452.07,
          "formula": "C16 H18 Ga N2 O4,N O3,H2 O", "series": "non-brominated", "isotope": "H"},
    "D": {"file": "D.dat", "mass_mg": 43.0, "compound": "Mp08-D", "mw": 456.10,
          "formula": "C16 H14 D4 Ga N2 O4,N O3,H2 O", "series": "non-brominated", "isotope": "D"},
}

def ensure_dirs():
    for p in [
        PROCESSED/"01_aluminium_5pt", PROCESSED/"02_al_subtracted",
        PROCESSED/"03_mass_normalised", PROCESSED/"04_anchor_normalised",
        PROCESSED/"05_three_point_smoothed", PROCESSED/"06_final_renormalised",
        PROCESSED/"07_difference", PROCESSED/"exploratory_cryo",
        META, PLOTS, ASSETS,
    ]:
        p.mkdir(parents=True, exist_ok=True)

def load_dat(path):
    a = np.loadtxt(path, delimiter=",", comments="#")
    return a[:, 0], a[:, 1], a[:, 2]

def centred_average(y, n):
    p = n // 2
    return np.convolve(np.pad(y, (p, p), mode="edge"), np.ones(n)/n, mode="valid")

def save_pair(x, final, keys, labels, title, filename):
    m = (x >= 50) & (x <= 1500)
    plt.figure(figsize=(11, 6))
    for k, label in zip(keys, labels):
        plt.plot(x[m], final[k][m], linewidth=1.4, label=label)
    plt.axhline(1.0, linewidth=0.8, alpha=0.25, linestyle="--")
    plt.xlim(50, 1500)
    plt.xlabel(r"Energy transfer (cm$^{-1}$)")
    plt.ylabel(r"Relative intensity (anchor peak = 1)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    out = PLOTS / filename
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    shutil.copy2(out, ASSETS / filename)

def save_difference(x, y, title, filename):
    m = (x >= 50) & (x <= 1500)
    plt.figure(figsize=(11, 5.5))
    plt.plot(x[m], y[m], linewidth=1.4)
    plt.axhline(0, linewidth=0.9, linestyle="--", alpha=0.6)
    plt.xlim(50, 1500)
    plt.xlabel(r"Energy transfer (cm$^{-1}$)")
    plt.ylabel("Relative intensity difference (H − D)")
    plt.title(title)
    plt.tight_layout()
    out = PLOTS / filename
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    shutil.copy2(out, ASSETS / filename)

def main():
    ensure_dirs()

    # Record supplied metadata.
    rows = []
    for key, s in SAMPLES.items():
        rows.append({
            "compound": s["compound"], "dataset": key, "raw_file": s["file"],
            "series": s["series"], "isotope": s["isotope"],
            "molecular_weight_g_mol": s["mw"], "formula_as_supplied": s["formula"],
            "mass_mg": s["mass_mg"], "amount_mmol": s["mass_mg"]/s["mw"],
        })
    pd.DataFrame(rows).to_csv(META/"samples.csv", index=False)

    x, y_al, e_al = load_dat(RAW/"emptyAl.dat")
    y_al_sm = centred_average(y_al, 5)
    ep = np.pad(e_al, (2, 2), mode="edge")
    e_al_sm = np.sqrt(np.convolve(ep**2, np.ones(5), mode="valid")) / 5.0

    pd.DataFrame({
        "energy_cm-1": x, "raw_intensity": y_al, "raw_uncertainty": e_al,
        "smoothed_intensity_5pt": y_al_sm,
        "smoothed_uncertainty_5pt": e_al_sm,
    }).to_csv(PROCESSED/"01_aluminium_5pt"/"emptyAl_5pt_centered.csv", index=False)

    mass_norm = {}
    for key, s in SAMPLES.items():
        _, y, e = load_dat(RAW/s["file"])
        yc = y - y_al_sm
        ec = np.sqrt(e**2 + e_al_sm**2)
        ym = yc / s["mass_mg"]
        em = ec / s["mass_mg"]
        mass_norm[key] = ym

        pd.DataFrame({
            "energy_cm-1": x, "intensity_al_subtracted": yc,
            "uncertainty_al_subtracted": ec,
        }).to_csv(PROCESSED/"02_al_subtracted"/f"{key}_al_subtracted.csv", index=False)

        pd.DataFrame({
            "energy_cm-1": x, "intensity_per_mg": ym,
            "uncertainty_per_mg": em,
        }).to_csv(PROCESSED/"03_mass_normalised"/f"{key}_mass_normalised.csv", index=False)

    # Exact sampled anchor: local maximum in 350–430 cm^-1.
    window = (x >= 350) & (x <= 430)
    anchors, first_norm, smooth3, final = {}, {}, {}, {}
    for key in SAMPLES:
        inds = np.where(window)[0]
        idx = inds[np.argmax(mass_norm[key][inds])]
        peak_y = mass_norm[key][idx]
        first_norm[key] = mass_norm[key] / peak_y
        anchors[key] = {
            "index": int(idx),
            "energy_cm-1": float(x[idx]),
            "mass_normalised_peak_height": float(peak_y),
        }
        pd.DataFrame({
            "energy_cm-1": x,
            "relative_intensity_anchor_1": first_norm[key],
        }).to_csv(PROCESSED/"04_anchor_normalised"/f"{key}_anchor_normalised.csv", index=False)

        smooth3[key] = centred_average(first_norm[key], 3)
        pd.DataFrame({
            "energy_cm-1": x,
            "relative_intensity_3pt_smoothed": smooth3[key],
        }).to_csv(PROCESSED/"05_three_point_smoothed"/f"{key}_3pt_smoothed.csv", index=False)

        anchor_after = float(smooth3[key][idx])
        anchors[key]["smoothed_anchor_before_renormalisation"] = anchor_after
        final[key] = smooth3[key] / anchor_after
        pd.DataFrame({
            "energy_cm-1": x,
            "relative_intensity_final": final[key],
        }).to_csv(PROCESSED/"06_final_renormalised"/f"{key}_final.csv", index=False)

    anchor_rows = []
    for key, s in SAMPLES.items():
        a = anchors[key]
        anchor_rows.append({
            "dataset": key, "compound": s["compound"], "mass_mg": s["mass_mg"],
            "anchor_peak_cm-1": a["energy_cm-1"],
            "mass_normalised_peak_height": a["mass_normalised_peak_height"],
            "smoothed_anchor_before_renormalisation": a["smoothed_anchor_before_renormalisation"],
        })
    pd.DataFrame(anchor_rows).to_csv(META/"normalisation_anchors.csv", index=False)

    diffs = {
        "Br-H_minus_Br-D": final["Br-H"] - final["Br-D"],
        "H_minus_D": final["H"] - final["D"],
    }
    for name, y in diffs.items():
        pd.DataFrame({
            "energy_cm-1": x,
            "relative_intensity_difference_H_minus_D": y,
        }).to_csv(PROCESSED/"07_difference"/f"{name}.csv", index=False)

    # Exploratory cryostat smoothing only.
    _, ycr, _ = load_dat(RAW/"emptycryo.dat")
    for n in (9, 19):
        pd.DataFrame({
            "energy_cm-1": x,
            "raw_cryostat_intensity": ycr,
            f"cryostat_{n}pt_centered_average": centred_average(ycr, n),
        }).to_csv(PROCESSED/"exploratory_cryo"/f"emptycryo_{n}pt_centered.csv", index=False)

    save_pair(x, final, ["Br-H","Br-D"], ["Mp50-H (Br-H)", "Mp50-D (Br-D)"],
              "Brominated pair: H vs D", "final_brominated_pair.png")
    save_pair(x, final, ["H","D"], ["Mp08-H (H)", "Mp08-D (D)"],
              "Non-brominated pair: H vs D", "final_non_brominated_pair.png")
    save_difference(x, diffs["Br-H_minus_Br-D"],
                    "Brominated pair: Br-H − Br-D", "difference_brominated.png")
    save_difference(x, diffs["H_minus_D"],
                    "Non-brominated pair: H − D", "difference_non_brominated.png")

    manifest = {
        "project": "Gallium – Peter Sadler TOSCA spectra",
        "dashboard_region_cm-1": [50, 1500],
        "final_pipeline": [
            "5-point centred average of empty aluminium",
            "subtract smoothed aluminium",
            "divide by sample mass in mg",
            "local maximum in 350–430 cm^-1",
            "normalize anchor to 1",
            "3-point centred average",
            "renormalize at same anchor grid point",
            "H-D differences",
        ],
        "anchors": {k: {kk: vv for kk, vv in a.items() if kk != "index"} for k,a in anchors.items()},
        "cryostat": {"raw_preserved": True, "exploratory_smoothing_points": [9,19],
                     "used_in_final_pipeline": False},
    }
    (META/"processing_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

if __name__ == "__main__":
    main()
