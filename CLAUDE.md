# FTIR Data Processing Project

## Overview

Python toolkit (`ftir_tools.py`) that replaces a suite of Excel templates for analyzing
FTIR (Fourier Transform Infrared) spectra of volcanic glass samples. Calculates dissolved
H2O (wt%) and CO2 (ppm) concentrations from transmission spectra.

## Project Files

```
FTIR/
  ftir_tools.py                              # Main module (all code in one file)
  co2_reference_spectra.json                 # CO3 reference spectra (from CO2_bestfit.xls)
  pyiroglass_pca_components.json             # PCA baseline/H2Om from PyIRoGlass (596 pts)
  sample_config.yaml                         # Configuration file template for users
  calculation_guide.txt                      # Standalone calculation logic document
  requirements.txt                           # Python dependencies (numpy, scipy, matplotlib, pyyaml, ...)
  README.md                                  # Quick-start guide
  CLAUDE.md                                  # This file — project documentation

  scripts/                                   # Run scripts
    run_example.py                           # Usage example (full pipeline + step-by-step)
    run_all_models.py                        # Batch comparison: fixed/pca_shift models
    run_pca_comparison.py                    # Model comparison + CSV output

  excel_templates/                           # Excel reference sources (read-only)
    Reflectance Method Template.xlsx         # Source for Module 1 formulas
    Glass density and index of refraction Template.xls  # Source for Module 2
    CO2_bestfit.xls                          # Source for Module 4 (CO3 fitting model)
    Master Table Template.xlsx               # Source for Module 5 (concentration formulas)
    FTIR fitting overview.docx               # Overview document

  example/                                   # Example spectra and output
    example spectra/                         # Transmission spectra (6 samples)
      reflectance spectra/                   # Reflectance spectra (4 samples)
    example_configs/                         # Example YAML config files (Mode B + C)
    ref good example.txt                     # Example reflectance spectrum
    trans good example.txt                   # Example transmission spectrum

  _archive/                                  # Archived experimental scripts + temp files
  verification_scripts/                      # Excel-vs-Python verification scripts (12 files)
  PyIRoGlass/                                # Reference: PyIRoGlass source + data
  Screenshots/                               # Reference screenshots
```

## Architecture

All code lives in `ftir_tools.py`, organized into 5 modules + pipeline functions:

| Module | Function(s) | Purpose | Excel Source |
|--------|------------|---------|--------------|
| 1 | `calculate_thickness()`, `calculate_thickness_multi()`, `calculate_thickness_from_spectrum()`, `olivine_refractive_index()` | Thickness from reflectance fringes | Reflectance Method Template |
| 2 | `calculate_density()`, `get_refractive_index()` | Glass density & n from oxide composition | Glass density Template |
| 3 | `fit_h2o_peak()`, `batch_fit_h2o()` | H2O peak baseline correction | — |
| 4 | `fit_carbonate()` | CO3 doublet fitting (3 models, see below) | CO2_bestfit.xls + PyIRoGlass |
| 5 | `calculate_concentration()` | Beer-Lambert concentration calc | Master Table |
| Pipeline | `process_sample()` | Full single-sample pipeline (M1-M5) | — |
| Pipeline | `process_batch()` | Batch processing with CSV/Excel output | — |
| Config | `load_config()`, `run_from_config()`, `run_from_folder()` | YAML/JSON config file support | — |

### Helper / internal functions

- `_read_spectrum_csv()` — reads tab/comma/space-delimited spectrum files
- `_split_iron()` — splits FeOT into Fe2O3 + FeO
- `_best_fringe_run()` — finds longest consecutive run of consistent fringe intervals
- `_load_reference_spectra()` — loads CO3 reference from JSON (for fixed model)
- `_load_pca_components()` — loads PyIRoGlass PCA baseline/H2Om components from JSON
- `_interpolate_to_reference_grid()` — interpolates sample spectrum onto reference grid
- `_generate_synthetic_spectrum()` — generates synthetic test data (Module 3 tests)
- `_generate_synthetic_co3_spectrum()` — generates synthetic CO3 test data (Module 4 tests)
- `_write_summary_csv()` / `_write_summary_excel()` / `_print_summary_table()` — output helpers
- `_validate()` — runs all built-in tests (Tests 1-6)

### CLI Usage

```bash
python ftir_tools.py                     # Run validation tests
python ftir_tools.py sample.yaml         # Process single sample from config
python ftir_tools.py configs/            # Process all configs in folder
python ftir_tools.py --validate          # Run validation tests
```

## Critical Design Decisions

### 1. 100 um Normalization in CO3 Fitting (IMPORTANT)

CO2_bestfit.xls normalizes the input spectrum to 100 um thickness **before** fitting:

```
normalized_spectrum = raw_spectrum * (100 / thickness_um)
```

- In the Excel: Column B = normalized spectrum, Column D = raw spectrum
- The ratio B/D = 100/thickness_um (constant for a given sample)
- The model fits **column B** (normalized), NOT column D (raw)
- Therefore, fitted coefficients `c_CO3` and `c_1630` correspond to a 100 um thick sample

This is why the Master Table concentration formulas for CO3 and 1630 use `100` as the
thickness, not the actual measured thickness:

```python
# Total H2O (from 3570 peak): uses ACTUAL thickness
total_h2o = 1e6 * A_OH * 18.02 / (62.3 * rho * d_um)

# H2Omol 1630 (from CO3 fit): uses 100 (normalized thickness)
h2omol_1630 = 1e6 * c_1630 * 18.02 / (25.0 * rho * 100)

# CO2 ppm (from CO3 fit): uses 100 (normalized thickness)
co2_ppm = 1e10 * 44.0 * c_CO3 / (375.0 * rho * 100)
```

### 2. CO3 Fitting Model

Model: `fit = offset + bkg*E + c1630*F + cCO3*G + st*H`

- 5 coefficients solved via `numpy.linalg.lstsq`
- Reference spectra (E, F, G, H) stored in `co2_reference_spectra.json`
  - E = background/baseline
  - F = H2O molecular bend at 1630 cm-1
  - G = CO3 doublet (~1430/1515 cm-1)
  - H = slope/tilt correction
- Fitting range: 1350-1800 cm-1 (116 points from 259-point reference grid, matching Excel rows 104-219)
- R-squared calculated over fitting range only: `R2 = SS(regression) / SS(total)`
- Python lstsq gives slightly better fit than Excel Solver (lower residual SS)
- **Multicollinearity**: Design matrix condition number ~11578 (very high). E_bkg and H_line are 93% correlated. This means different coefficient combinations can produce nearly identical fitted curves (<0.6% max difference). The CO3 coefficient (most important for CO2 ppm) is the most stable (~4.6% diff vs Excel).

### 3. CO3 Doublet Fitting — Two Models

Two models are available for CO3 fitting via `model` parameter in `fit_carbonate()`:

1. `model="fixed"` (default) — matches Excel CO2_bestfit.xls
2. `model="pca_shift"` — PCA baseline + Taylor-expanded Gaussians with independent shifts

**Background**: The 1515 and 1430 cm-1 CO3 peaks arise from ν3 asymmetric stretching
of CO3²⁻ split by the local site environment. Different cation environments (Na⁺, Ca²⁺,
Mg²⁺) affect each mode differently, so the shifts are NOT necessarily equal.

**Comparison of approaches**:

| | Fixed | PCA+Shift | PyIROGlass |
|--|--|--|--|
| Peak position | Fixed reference | Independent shift per peak | Independent μ₁, μ₂ |
| Baseline | E_bkg + H_line | PCA (5 comp) | PCA (5 comp) |
| Shift DOF | 0 (sweep only) | 2 (Taylor) | 2 (MCMC) |
| Model type | Linear (lstsq) | Bounded linear | Nonlinear + MCMC |
| Speed | Instant | Instant | ~20s/spectrum |
| Fit range | 1350-1800 | 1250-2200 | 1250-2400 |
| Params | 5 | 14 | ~20 |

**Reference**: Shi et al. (2024) PyIRoGlass, doi:10.30909/vol.07.02.471501

### 3a. PCA Baseline Details (used by `model="pca_shift"`)

**Implemented** following Shi et al. (2024) PyIRoGlass approach.

**PCA training data source**: We do NOT train PCA ourselves. We use **pre-computed PCA
components** from PyIRoGlass (files `BaselineAvgPC.npz` and `H2Om1635PC.npz` from their
GitHub repo, converted to `pyiroglass_pca_components.json`). PyIRoGlass trained these on
**57 devolatilized glass spectra** (basalt-andesite, SiO2=43-60%) stored in `All_BDL.csv`.

**Model components** (12 parameters, all solved simultaneously):
- **PCA baseline (5 comp)**: 1 mean + 4 principal components — captures glass absorption
  including the Si-O upturn below 1350 cm⁻¹ that the fixed model cannot handle
- **H2Om,1635 (3 comp)**: 1 empirical mean peak + 2 PCs — replaces our single F_1630
  Gaussian, better captures the molecular H2O peak shape variations
- **CO3 (2 Gaussians)**: independent peaks at 1430 and 1515 cm⁻¹ (σ=30 cm⁻¹)
- **Linear tilt + offset**: 2 additional free parameters

**Why Gaussians instead of G_CO3 reference**: The full G_CO3 reference from Excel has
97.7% of its energy OUTSIDE the CO3 peaks (large negative background). This creates
r=0.79 correlation with PCA baseline components, causing severe degeneracy in
unconstrained fits. Two narrow Gaussians (σ=30) reduce max correlation to ~0.48,
making bounded lstsq stable.

**Solver**: `scipy.optimize.lsq_linear` (bounded linear least squares).

**Bounds** (adapted from PyIRoGlass MCMC bounds, core.py lines 625-628):
- B_mean [0, 4], BPC1 [-3, 3], BPC2 [-2, 2], BPC3 [-0.6, 0.6], BPC4 [-0.3, 0.3]
- H2Om [0, 3], H2Om_PC1/2 [-2, 2]
- CO3 amplitudes [0, inf], tilt [-0.5, 0.5], offset [-1, 3]

**Fitting range**: 1250-2200 cm⁻¹ (wider than fixed model's 1350-1800, covers Si-O upturn)

**Results** (6 example spectra, JT12 composition):
- 4/6 GOOD quality (vs 1/6 for fixed model)
- R-squared higher across all samples (PCA captures Si-O upturn)
- CO2 values differ 10-20% from fixed model (different baseline treatment)
- No false positives (correctly gives 0 CO3 for samples without CO3)

### 3b. PCA + Independent Peak Shift (`model="pca_shift"`)

PCA baseline (Section 3a) with **Taylor-expanded Gaussians** for independent
peak shift estimation. Each Gaussian gets a derivative basis function:
`G'(nu; mu, sigma) = G(nu) * (nu - mu) / sigma^2`

**14 parameters**: 5 PCA baseline + 3 H2Om + 2 CO3 amplitudes + 2 CO3 derivatives + tilt + offset

After fitting, independent shifts are recovered:
- `shift_1430 = c_deriv_1430 / c_amp_1430`
- `shift_1515 = c_deriv_1515 / c_amp_1515`

**SUSPECT flag**: If |shift| > 15 cm-1 for either peak, the result is flagged as
"SUSPECT" — the feature being fitted is likely NOT dissolved CO3 (could be interference
fringes, organic contamination, or mineral inclusions). This provides automatic quality
control that PyIRoGlass's hard bounds cannot offer.

**Comparison script**: `run_pca_comparison.py` generates figures and CSV for all 4 models
in `example/example spectra/pca_comparison/`.

### 4. H2O Peak Baseline Correction

Linear baseline through two anchor regions:

- **Low anchor**: 2200-2400 cm-1 (flat, absorption-free region)
- **High anchor**: 3700-3800 cm-1 (above the O-H peak)

The peak height is the maximum baseline-corrected absorbance in the peak region
(typically near 3570 cm-1 for OH stretch, or ~3532 cm-1 depending on composition).

Previously used 3000-3800 cm-1 with both anchors at edges. The wider 2200-3800 range
with explicit low anchor at 2200-2400 is more robust.

### 5. H2O Uncertainty Estimation

Two sources, combined in quadrature:

1. **Low anchor variation**: The 2200-2400 cm-1 range is divided into overlapping
   sub-windows (~67 cm-1 wide). Each sub-window gives a slightly different baseline
   anchor point, leading to different peak heights. The std of these peak heights
   propagates to H2O wt% uncertainty.

2. **Thickness uncertainty**: Since H2O wt% = k * A / d, relative uncertainty in
   thickness directly propagates: `dH2O/H2O = dd/d`.

Previously used a narrow (3000-3800 cm-1) vs wide (2200-3800 cm-1) baseline
comparison, but the narrow range low anchor (3000-3100 cm-1) sits on the OH
peak shoulder and is unreliable. Dropped in favor of the current approach.

### 6. Thickness Measurement — Three Input Modes

`process_sample()` supports three mutually exclusive thickness input modes:

| Mode | Parameters | Description |
|------|-----------|-------------|
| A | `thickness_measurements=[(wn_hi, wn_lo, N), ...]` | Manual fringe counts (original Excel method) |
| B | `thickness_um=27.4, thickness_unc_um=0.5` | Direct thickness value |
| C | `reflectance_csv="ref.txt", olivine_fo=85` | Auto-detect from reflectance spectrum |

**Mode C algorithm** (`calculate_thickness_from_spectrum()`):

1. Compute olivine refractive index from Fo content:
   `n = 1.858 - 0.206 * Fo` (linear interpolation, Deer et al. 1992)
2. Read reflectance spectrum, extract 1800-3000 cm-1 (clean, absorption-free)
3. Savitzky-Golay smoothing (window=31, order=3) to remove noise
4. Detect maxima and minima independently via `scipy.signal.find_peaks`
5. **Best-fringe-run selection** (`_best_fringe_run()`):
   - For each consecutive pair of extrema, compute per-interval thickness
   - Find the longest consecutive sub-sequence where all per-interval values
     are within 20% of the sub-sequence median
   - Use only that sub-region for the final thickness calculation
6. Thickness from selected maxima and minima independently, then averaged
7. Uncertainty = stdev of per-interval thicknesses within selected sub-regions

**Why best-fringe-run**: Weak or noisy spectra may have false peaks from
baseline drift or absorption features. The algorithm automatically identifies
the cleanest sub-region with the most regular fringe spacing, rejecting
outlier intervals. This dramatically reduces uncertainty (e.g. from +/-14 to
+/-1.5 um on noisy spectra).

**Uncertainty propagation**: All three modes produce `stdev_um`, which
propagates into `calculate_concentration()` via:
`sigma_C/C = sqrt((sigma_A/A)^2 + (sigma_d/d)^2 + (sigma_eps/eps)^2)`

## Physical Constants & Defaults

| Parameter | Value | Notes |
|-----------|-------|-------|
| epsilon_OH (3570) | 62.3 L/(mol*cm) | Molar absorptivity |
| epsilon_H2Omol (1630) | 25.0 L/(mol*cm) | Molar absorptivity |
| epsilon_CO3 | 375.0 L/(mol*cm) | Molar absorptivity |
| Fe3+/Fe_total | 0.15 | Default ratio for iron splitting |
| MW H2O | 18.02 g/mol | |
| MW CO2 | 44.0 g/mol | |
| MW SiO2 | 60.09 | Used in density calculation |

## Running Tests

```python
python -c "from ftir_tools import _validate; _validate()"
```

Or run `ftir_tools.py` directly (calls `_validate()` in `if __name__ == '__main__'`):

```
python ftir_tools.py
```

Tests 1-6 validate each module against known Excel results.

## Input File Format

Spectrum files are tab/comma/space-delimited text with two columns:
- Column 1: wavenumber (cm-1)
- Column 2: absorbance (or reflectance for thickness measurement)

Header lines starting with non-numeric characters are automatically skipped.
Files can use `.txt`, `.csv`, or `.dat` extensions.

## Dependencies

- numpy
- scipy (interp1d, minimize)
- matplotlib (for plotting)
- openpyxl (optional, for Excel output in process_batch)

## Analysis & Verification Scripts (preserved for reference)

These scripts were created during Excel-vs-Python verification and are kept for future use:

| Script | Purpose |
|--------|---------|
All scripts are in `verification_scripts/` folder:

| Script | Purpose |
|--------|---------|
| `extract_vba.py` | Extract OLE streams from CO2_bestfit.xls (raw olefile approach) |
| `extract_vba2.py` | Extract VBA macro source code using oletools/olevba |
| `extract_solver.py` | Deep analysis of CO32-calc sheet structure + Solver named ranges |
| `check_excel.py` | Survey all 4 Excel files: sheets, macros, formula cells |
| `check_f1630_json.py` | Verify F_1630 reference (229 pts + 30 zeros) in JSON vs Excel |
| `compare_solver.py` | First comparison (discovered Excel has DIFFERENT sample than trans example) |
| `compare_same_data.py` | Feed Excel's own spectrum into Python for true head-to-head |
| `find_excel_fitrange.py` | Pinpoint Excel's exact fitting range (rows 104-219, 116 pts) |
| `compare_exact_range.py` | Compare with matching fitting range [1350, 1800] |
| `final_comparison.py` | Definitive comparison: same data, same range, multicollinearity analysis |
| `verify_all_formulas.py` | Cross-check all 4 Excel templates against Python modules |
| `verify_density.py` | JT12 density + Master Table concentration formula match |
| `count_fringes.py` | Auto-detect interference fringes from reflectance spectrum |
| `count_fringes2.py` | Fringe detection with proper spacing for ~30 um samples |

### Key findings from verification

1. **VBA macro** (`Module1` in CO2_bestfit.xls): Only a `Clean_up` data preprocessing macro (downsamples spectrum by deleting every other row). NOT the Solver optimization logic.

2. **Solver config** (stored in Excel named ranges, not VBA):
   - `solver_opt=2` (minimize), `solver_lin=1` (linear model)
   - `solver_cvg=0.0001`, `solver_pre=1e-6`, `solver_itr=100`, `solver_tol=0.05`
   - Target cell: SS(E) in O1; Adjustable cells: coefficients in J1:J5

3. **Model structure verified** to machine precision (max diff = 1.11e-16):
   `model = offset + bkg*E + c1630*F + cCO3*G + st*H`

4. **Reference spectra (JSON vs Excel)**: Exact match (diff = 0 for all columns).

5. **Normalization verified**: Col B / Col D = 100/d_um = 3.322259 (constant).

6. **Fitting range**: Excel uses rows 104-219 = [1353.81, 1797.36] cm-1 = 116 points.
   Python default updated from (1354, 1797) to (1350, 1800) to match.

7. **Coefficient differences** are due to multicollinearity (condition # = 11578).
   Python lstsq gives mathematically optimal solution (lower SS(E)).
   Despite very different coefficients, fitted curves differ by only 0.6% max.

8. **All 4 Excel templates verified**:
   - Reflectance Method: formula exact match
   - Glass density: diff < 1.2e-6 g/cc (JT12)
   - Master Table concentrations: diff < 1e-10 (Q, R, V formulas)
   - CO2_bestfit: model structure exact, coefficients differ due to multicollinearity

## Development Conventions

- All code in a single file (`ftir_tools.py`) — no separate module files
- Tests are built-in via `_validate()`, not in separate test files
- Diagnostic plots saved as PNG alongside input files
- Use tab-delimited text for spectrum data
- Windows environment: use forward slashes in paths for bash commands
- Avoid Unicode characters in print statements (cp1252 encoding on Windows)
