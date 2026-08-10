# GIWAXS/GIXS Indexing Workbench

A Python-based workflow for converting rendered GIWAXS/GIXS images into reciprocal-space data and comparing experimental diffraction peaks with calculated crystallographic reflections.

The repository contains two graphical programs:

- **`GIWAXS_PNG_to_NPZ_Converter.py`** converts one or more GIWAXS/GIXS PNG images into `.npz` files containing reciprocal-space axes and reconstructed intensity data.
- **`GIWAXS_GIXS_Workbench.py`** provides automatic indexing and a manual experimental/calculated peak-comparison interface using CIF, NPZ, and PNG inputs.

## Project Workflow

A typical workflow is:

1. Convert GIWAXS/GIXS PNG images to NPZ using the PNG-to-NPZ converter when original numerical reciprocal-space data is not available.
2. Load a **CIF** file containing the crystal structure.
3. Load **NPZ** data for numerical reciprocal-space analysis.
4. Load a **PNG** image as the visual experimental overlay.
5. Run the indexing calculation and review calculated versus experimental reflections.
6. Use the manual comparison tools when additional inspection of peak assignments is needed.

## Data Types

### CIF

The CIF file provides the crystal structure, lattice information, space group, and calculated crystallographic reflections used by the indexing workbench.

### NPZ

NPZ files provide numerical reciprocal-space information used by the automatic calculations.

The PNG-to-NPZ converter creates files containing at least:

- `intensity` — reconstructed 2D intensity
- `qr` — horizontal reciprocal-space axis
- `qz` — vertical reciprocal-space axis
- `mask` — validity mask

Additional conversion metadata may also be stored in the NPZ file.

### PNG

PNG files are used as the visual experimental GIWAXS/GIXS overlay in the workbench.

> **Important:** When NPZ intensity is reconstructed from a rendered PNG, the intensity values are approximate because they are recovered from image colors. Original numerical reciprocal-space data should be used when available.

## GIWAXS/GIXS Workbench

The workbench combines automatic indexing with manual peak comparison.

Main features include:

- CIF-based crystallographic reflection generation
- Automatic GIWAXS/GIXS indexing
- Experimental-to-calculated peak comparison
- Manual experimental peak selection
- Reciprocal-space coordinate comparison
- Calculated and experimental intensity comparison
- d-spacing display
- HKL reflection assignments
- Simple and Advanced table views

### Simple Table View

The default **Simple view** presents the main results in more readable terms, including:

- Peak number
- Crystal plane (HKL)
- Measured and predicted reciprocal-space coordinates
- Measured and predicted d-spacing
- Match error (Δq)
- Measured and predicted intensity
- Intensity match
- Match quality

### Advanced Table View

The **Advanced view** exposes the more detailed scientific and diagnostic fields used by the program.

Changing between Simple and Advanced view changes only how results are displayed; it does not change the indexing calculations.

## PNG-to-NPZ Converter

The converter provides a graphical interface for batch conversion of rendered GIWAXS/GIXS PNG images.

Features include:

- Selecting multiple PNG files at once
- Selecting an entire folder of PNG files
- Entering reciprocal-space limits once for a batch
- Automatic detection of the colored q-space panel
- Automatic identification of common Matplotlib colormaps
- Optional crop preview
- Saved settings between runs
- Compressed NPZ output

Converted files are saved next to the source PNG with names ending in:

```text
_qspace_axes.npz
```

## Requirements

Python 3.10 or newer is recommended.

Required Python packages:

```text
numpy>=1.24
pandas>=2.0
scipy>=1.10
matplotlib>=3.7
Pillow>=10.0
gemmi>=0.6.7
cloudpickle>=3.0
PyQt6>=6.6
```

Tkinter is also required for the PNG-to-NPZ GUI. It is normally included with standard Python installations on Windows and macOS.

## Installation

Clone or download the repository, open a terminal in the project folder, and install the dependencies:

```bash
python -m pip install -r requirements.txt
```

The two Python programs also check for missing pip-installable dependencies when they start and attempt to install them into the active Python environment.

## Running in PyCharm

1. Open the repository folder as a PyCharm project.
2. Select a Python interpreter for the project.
3. Install the packages from `requirements.txt` if they are not already installed.
4. Run either program directly:

```text
GIWAXS_PNG_to_NPZ_Converter.py
GIWAXS_GIXS_Workbench.py
```

The graphical interface will open in a separate desktop window.

## Running from a Terminal

PNG-to-NPZ converter:

```bash
python GIWAXS_PNG_to_NPZ_Converter.py
```

GIWAXS/GIXS workbench:

```bash
python GIWAXS_GIXS_Workbench.py
```

## Repository Files

```text
GIWAXS_GIXS_Workbench.py
GIWAXS_PNG_to_NPZ_Converter.py
requirements.txt
README.md
.gitignore
```

These five files are the recommended core GitHub repository.

## Notes

- The workbench is designed for GIWAXS/GIXS reciprocal-space analysis and should be used with appropriate experimental calibration and crystallographic inputs.
- A good positional match between an experimental peak and a calculated reflection supports an assignment but does not by itself establish crystallographic ground truth.
- PNG-derived numerical intensity should be treated as reconstructed/approximate data.
