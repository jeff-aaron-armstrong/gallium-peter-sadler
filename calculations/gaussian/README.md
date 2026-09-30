# Gaussian reference calculations

These files were supplied as **relevant structure and vibrational-frequency calculations** for the gallium complexes measured on TOSCA.

## Mapping to the experimental families

| Calculation label | Experimental family | Calculated cation | Bromination |
|---|---|---|---|
| R1 | Mp08 | C16H18GaN2O4+ | non-brominated |
| R2 | Mp50 | C16H16Br2GaN2O4+ | dibrominated |

The mapping is based directly on the atom content of the supplied structures: R2 contains two Br atoms, whereas R1 contains none.

## Calculation details

Both frequency calculations are Gaussian 16 jobs using:

- **B3LYP**
- **CEP-31G**
- **Grimme D3 dispersion** (`empiricaldispersion=gd3`)
- `Freq=noRaman`
- charge **+1**
- multiplicity **1**

Both jobs terminated normally and contain **117 vibrational modes**, with **no imaginary frequencies**.

| System | Lowest mode (cm⁻¹) | Highest mode (cm⁻¹) | Imaginary modes |
|---|---:|---:|---:|
| R1 non-brominated | 22.3229 | 3837.7700 | 0 |
| R2 brominated | 9.0358 | 3840.9860 | 0 |

## Important scope/caveats

These Gaussian models are the **cation only**. They do **not** contain the nitrate counterion present in the experimental sample formula, nor a periodic crystal environment.

The supplied calculations are also **protiated**. No deuterium isotope masses were specified in the Gaussian frequency jobs. They can therefore serve as a useful starting point for mode assignment and for generating isotope-substituted calculations, but they are not direct D calculations.

## Files committed here

- `R1_non_brominated/R1_cat_Fr.pdb` — supplied R1 structure.
- `R1_non_brominated/frequencies_and_ir.csv` — all 117 Gaussian frequencies and IR intensities parsed from the supplied R1 log.
- `R2_brominated/R2_cat_Fr.pdb` — supplied R2 structure.
- `R2_brominated/frequencies_and_ir.csv` — all 117 Gaussian frequencies and IR intensities parsed from the supplied R2 log.

The original upload also contains the full Gaussian normal-mode logs:

- `R1_cat_Freq.log` — SHA-256 `b01c2092489b9c2c60ee0532430b1eb61b40a3984fcb8fe760080493592b8cbc`
- `R2_cat_Fr.log` — SHA-256 `591e30eb76a3fdc8f2c500f8cdb0b7068e632a4559f2dcfb2a75915b65e115e2`

Those logs include the normal-mode displacement vectors as well as the frequencies and are the most useful source files for a subsequent INS/AbINS calculation.

Original archive: `Fw__Cifs_for_Ga_for_DFT.zip`  
Archive SHA-256: `b32bb206bc5be9784b327a3db1c75c5f847323f8542f24c0e73da607d7b8c80e`
