# FTIR Vol-Cal

Python toolkit for analyzing dissolved H2O and CO2 in volcanic glass samples from FTIR spectra. Please check out the Web app version if it is easier to use: https://ftir-vol-cal.onrender.com/

Replaces the following Excel templates for H2O and CO2 calculation, which you can find in the repository :
- Reflectance Method Template (thickness)
- Glass density and index of refraction Template (density, refractive index)
- CO2_bestfit.xls (CO3 carbonate doublet fitting)
- Master Table Template (Beer-Lambert concentration calculations)

## Local Installation

### Requirements

- Python 3.9+

### Setup

```bash
git clone https://github.com/xinkaihe26/FTIR-Vol-cal.git
cd FTIR-Vol-cal
pip install -r requirements.txt
```

### Option 1: Web Interface (recommended)

```bash
cd web
python app.py
```

Open http://localhost:5000 in your browser.

### Option 2: Command Line

Process a single sample from a YAML config file:

```bash
python ftir_tools.py sample_config.yaml
```

Process all config files in a folder:

```bash
python ftir_tools.py configs/
```

### Option 3: Python Script

```python
from ftir_tools import process_sample

result = process_sample(
    spectrum_csv="path/to/spectrum.txt",
    composition_dict={
        "SiO2": 50.28, "TiO2": 2.58, "Al2O3": 14.78,
        "FeOT": 5.32, "MgO": 5.29, "CaO": 10.97,
        "Na2O": 2.08, "K2O": 1.32, "H2O": 0.686,
    },
    thickness_um=27.4,
    thickness_unc_um=0.5,
)

print(f"Total H2O: {result['concentration']['total_h2o_wt_pct']:.2f} wt%")
print(f"CO2: {result['concentration']['co2_ppm']:.0f} ppm")
```

See `scripts/run_example.py` for more detailed examples.

## Features

- Three thickness modes: manual fringe counts, direct value, or auto-detection from reflectance spectrum
- Three CO3 fitting models: simple lstsq, lstsq with shifted peak position, PCA baseline
- Automatic fringe correction for CO3 fitting
- Composition-dependent molar absorptivities (Shi et al. 2024)
- H2O and CO2 uncertainty estimation
- Single sample and batch processing
- Diagnostic plots for quality control

## Running Tests

```bash
python ftir_tools.py --validate
```

## Project Structure

```
ftir_tools.py                    # Main module (all code)
co2_reference_spectra.json       # CO3 reference spectra
pyiroglass_pca_components.json   # PCA components for pca_shift model
sample_config.yaml               # Config file template
requirements.txt                 # Python dependencies
web/                             # Flask web application
  app.py                         # Backend
  templates/index.html           # Frontend
  static/                        # CSS and JS
scripts/                         # Example scripts
example/                         # Example spectra and output
```

## Modules

| Module | Function | Purpose |
|--------|----------|---------|
| 1 | `calculate_thickness()` | Thickness from reflectance fringes |
| 2 | `calculate_density()`, `get_refractive_index()` | Glass density & refractive index |
| 3 | `fit_h2o_peak()` | H2O peak baseline correction |
| 4 | `fit_carbonate()` | CO3 doublet fitting (3 models) |
| 5 | `calculate_concentration()` | Beer-Lambert concentrations |
| Pipeline | `process_sample()` | Full single-sample pipeline |
| Pipeline | `process_batch()` | Batch processing with CSV/Excel output |
