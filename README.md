# FTIR Volatile Analysis Tools

Python toolkit for analyzing dissolved H2O and CO2 in volcanic glass samples from FTIR spectra.

Replaces the following Excel templates:
- Reflectance Method Template (thickness)
- Glass density and index of refraction Template (density, refractive index)
- CO2_bestfit.xls (CO3 carbonate doublet fitting)
- Master Table Template (Beer-Lambert concentration calculations)

## Installation

```bash
pip install -r requirements.txt
```

## Quick Start

```bash
python run_example.py
```

This runs a complete analysis on the example spectrum (`example/trans example.txt`)
using JT12 glass composition, printing all results and saving diagnostic plots to `example/`.

## Usage

```python
from ftir_tools import process_sample

result = process_sample(
    spectrum_csv="path/to/spectrum.txt",
    composition_dict={
        "SiO2": 50.28, "TiO2": 2.58, "Al2O3": 14.78,
        "FeOT": 5.32, "MgO": 5.29, "CaO": 10.97,
        "Na2O": 2.08, "K2O": 1.32, "H2O": 0.686,
    },
    thickness_measurements=[
        (2800, 2100, 3),   # (wavenumber_high, wavenumber_low, num_fringes)
        (3200, 2100, 5),
        (3600, 2100, 7),
    ],
)

print(f"Total H2O: {result['concentration']['total_h2o_wt_pct']:.2f} wt%")
print(f"CO2: {result['concentration']['co2_ppm']:.0f} ppm")
```

See `run_example.py` for detailed usage of individual modules.

## Running Tests

```bash
python -c "from ftir_tools import _validate; _validate()"
```

## Project Structure

```
ftir_tools.py              # Main module (all code)
co2_reference_spectra.json # CO3 reference spectra
run_example.py             # Usage example script
example/                   # Example spectra and output plots
verification_scripts/      # Excel-vs-Python verification scripts
CLAUDE.md                  # Detailed project documentation
```

## Modules

| Module | Function | Purpose |
|--------|----------|---------|
| 1 | `calculate_thickness()` | Thickness from reflectance fringes |
| 2 | `calculate_density()`, `get_refractive_index()` | Glass density & n |
| 3 | `fit_h2o_peak()` | H2O peak baseline correction |
| 4 | `fit_carbonate()` | CO3 doublet fitting (least-squares) |
| 5 | `calculate_concentration()` | Beer-Lambert concentrations |
| Pipeline | `process_sample()` | Full single-sample pipeline |
| Pipeline | `process_batch()` | Batch processing with CSV/Excel output |
