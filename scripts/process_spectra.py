#!/usr/bin/env python3
"""Reproduce the Gallium / Peter Sadler TOSCA processing pipeline.

Final relative-shape pipeline
-----------------------------
1. Centred 5-point arithmetic mean of the empty aluminium spectrum.
2. Subtract the smoothed aluminium spectrum from each sample.
3. Identify the common ~400 cm^-1 anchor as the sampled local maximum in 350–430 cm^-1.
4. Apply a centred 3-point arithmetic mean to each aluminium-subtracted sample spectrum.
5. Normalize ONCE by the 3-point-smoothed intensity at that fixed anchor position.
6. Calculate H-D differences separately for brominated and non-brominated pairs.

A mass-normalised branch is exported separately for absolute-intensity comparison, but it is
not part of the final relative-shape normalization because any constant mass scale cancels
when a spectrum is divided by its own anchor intensity.

The short empty-cryostat run is preserved and exploratory 9/19-point smooths are exported,
but no cryostat subtraction is used in the final spectra.
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
    # processed/ is fully generated: rebuild it cleanly so obsolete exploratory stages
    # from earlier versions do not remain in the repository.
    if PROCESSED.exists():
        shutil.rmtree(PROCESSED)
    for p in [
        PROCESSED/"01_aluminium_5pt",
        PROCESSED/"02_al_subtracted",
        PROCESSED/"03_mass_normalised",
        PROCESSED/"04_three_point_smoothed",
        PROCESSED/"05_final_anchor_normalised",
        PROCESSED/"06_difference",
        PROCESSED/"exploratory_cryo",
        META, PLOTS, ASSETS,
    ]:
        p.mkdir(parents=True, exist_ok=True)

def load_dat(path):
    a = np.loadtxt(path, delimiter=",", comments="#")
    return a[:, 0], a[:, 1], a[:, 2]

def centred_average(y, n):
    p = n // 2
    return np.convolve(np.pad(y, (p, p), mode="edge"), np.ones(n)/n, mode="valid")

def centred_average_uncertainty(e, n):
    p = n // 2
    ep = np.pad(e, (p, p), mode="edge")
    return np.sqrt(np.convolve(ep**2, np.ones(n), mode="valid")) / n

def save_pair(x, final, keys, labels, title, filename):
    m = (x >= 50) & (x <= 1500)
    plt.figure(figsize=(11, 6))
    for k, label in zip(keys, labels):
        plt.plot(x[m], final[k][m], linewidth=1.4, label=label)
    plt.axhline(1.0, linewidth=0.8, alpha=0.25, linestyle="--")
    plt.xlim(50, 1500)
    plt.xlabel(r"Energy transfer (cm$^{-1}$)")
    plt.ylabel(r"Relative intensity (smoothed anchor = 1)")
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

    rows = []
    for key, s in SAMPLES.items():
        rows.append({
            "compound": s["compound"], "dataset": key, "raw_file": s["file"],
            "series": s["series"], "isotope": s["isotope"],
            "molecular_weight_g_mol": s["mw"], "formula_as_supplied": s["formula"],
            "mass_mg": s["mass_mg"], "amount_mmol": s["mass_mg"]/s["mw"],
        })
    pd.DataFrame(rows).to_csv(META/"samples.csv", index=False)

    # Background: smooth empty aluminium with a centred 5-point average.
    x, y_al, e_al = load_dat(RAW/"emptyAl.dat")
    y_al_sm = centred_average(y_al, 5)
    e_al_sm = centred_average_uncertainty(e_al, 5)

    pd.DataFrame({
        "energy_cm-1": x,
        "raw_intensity": y_al,
        "raw_uncertainty": e_al,
        "smoothed_intensity_5pt": y_al_sm,
        "smoothed_uncertainty_5pt": e_al_sm,
    }).to_csv(PROCESSED/"01_aluminium_5pt"/"emptyAl_5pt_centered.csv", index=False)

    al_sub, al_sub_err = {}, {}

    # Main sample correction. Mass-normalised files are also exported, but they form
    # a parallel absolute-intensity branch rather than an input to the final normalization.
    for key, s in SAMPLES.items():
        _, y, e = load_dat(RAW/s["file"])
        yc = y - y_al_sm
        ec = np.sqrt(e**2 + e_al_sm**2)
        al_sub[key], al_sub_err[key] = yc, ec

        pd.DataFrame({
            "energy_cm-1": x,
            "intensity_al_subtracted": yc,
            "uncertainty_al_subtracted": ec,
        }).to_csv(PROCESSED/"02_al_subtracted"/f"{key}_al_subtracted.csv", index=False)

        pd.DataFrame({
            "energy_cm-1": x,
            "intensity_per_mg": yc / s["mass_mg"],
            "uncertainty_per_mg": ec / s["mass_mg"],
        }).to_csv(PROCESSED/"03_mass_normalised"/f"{key}_mass_normalised.csv", index=False)

    # Fix the anchor position using the sampled local maximum in the Al-subtracted
    # spectrum between 350 and 430 cm^-1. Dividing by mass would not change this position.
    window = (x >= 350) & (x <= 430)
    anchors = {}
    for key in SAMPLES:
        inds = np.where(window)[0]
        idx = inds[np.argmax(al_sub[key][inds])]
        anchors[key] = {"index": int(idx), "energy_cm-1": float(x[idx])}

    # Light sample smoothing FIRST, followed by the ONE normalization used in the
    # final relative spectra.
    smooth3, smooth3_err, final = {}, {}, {}
    for key in SAMPLES:
        ys = centred_average(al_sub[key], 3)
        es = centred_average_uncertainty(al_sub_err[key], 3)
        smooth3[key], smooth3_err[key] = ys, es

        pd.DataFrame({
            "energy_cm-1": x,
            "intensity_3pt_smoothed": ys,
            "uncertainty_3pt_smoothed": es,
        }).to_csv(PROCESSED/"04_three_point_smoothed"/f"{key}_3pt_smoothed.csv", index=False)

        idx = anchors[key]["index"]
        anchor_height = float(ys[idx])
        anchors[key]["smoothed_anchor_height"] = anchor_height

        yf = ys / anchor_height
        ef = es / abs(anchor_height)
        final[key] = yf

        pd.DataFrame({
            "energy_cm-1": x,
            "relative_intensity_final": yf,
            "uncertainty_scaled_by_anchor": ef,
        }).to_csv(PROCESSED/"05_final_anchor_normalised"/f"{key}_final.csv", index=False)

    anchor_rows = []
    for key, s in SAMPLES.items():
        a = anchors[key]
        anchor_rows.append({
            "dataset": key,
            "compound": s["compound"],
            "mass_mg": s["mass_mg"],
            "anchor_peak_cm-1": a["energy_cm-1"],
            "smoothed_anchor_height_before_normalisation": a["smoothed_anchor_height"],
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
        }).to_csv(PROCESSED/"06_difference"/f"{name}.csv", index=False)

    # Exploratory cryostat smoothing only; not used in the final pipeline.
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
        "final_relative_pipeline": [
            "5-point centred average of empty aluminium",
            "subtract smoothed aluminium",
            "identify anchor position as local maximum in 350–430 cm^-1",
            "3-point centred average of aluminium-subtracted sample spectrum",
            "single normalization by smoothed intensity at the fixed anchor position",
            "H-D differences",
        ],
        "mass_normalised_branch": "Exported separately for absolute-intensity comparison; not used in final relative-shape normalization.",
        "anchors": {k: {kk: vv for kk, vv in a.items() if kk != "index"} for k,a in anchors.items()},
        "cryostat": {
            "raw_preserved": True,
            "exploratory_smoothing_points": [9,19],
            "used_in_final_pipeline": False
        },
    }
    (META/"processing_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

if __name__ == "__main__":
    main()
