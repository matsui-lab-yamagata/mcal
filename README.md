# mcal: Program for the calculation of mobility tensor for organic semiconductor crystals
[![Python](https://img.shields.io/badge/python-3.7%20or%20newer-blue)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

# Overview
`mcal.py` is a tool for calculating mobility tensors of organic semiconductors. It calculates transfer integrals and reorganization energy from crystal structures, and determines mobility tensors considering anisotropy and path continuity.

# Requirements
* Python 3.7 or newer
* NumPy
* Pandas
* Gaussian 09 or 16

# Important notice
* The path of the Gaussian must be set.

# mcal Usage Manual

## Basic Usage

```bash
python mcal.py <cif_filename> <osc_type> [options]
```

### Required Arguments

- `cif_filename`: Path to the CIF file
- `osc_type`: Organic semiconductor type
  - `p`: p-type semiconductor (uses HOMO level)
  - `n`: n-type semiconductor (uses LUMO level)

### Basic Examples

```bash
# Calculate as p-type semiconductor
python mcal.py xxx.cif p

# Calculate as n-type semiconductor
python mcal.py xxx.cif n
```

## Options

### Calculation Settings

#### `-M, --method <method>`
Specify the calculation method used in Gaussian calculations.
- **Default**: `B3LYP/6-31G(d,p)`
- **Example**: `python mcal.py xxx.cif p -M "B3LYP/6-31G(d)"`

#### `-c, --cpu <number>`
Specify the number of CPUs to use.
- **Default**: `4`
- **Example**: `python mcal.py xxx.cif p -c 8`

#### `-m, --mem <memory>`
Specify the amount of memory in GB.
- **Default**: `10`
- **Example**: `python mcal.py xxx.cif p -m 16`

#### `-g, --g09`
Use Gaussian 09 (default is Gaussian 16).
- **Example**: `python mcal.py xxx.cif p -g`

### Calculation Control

#### `-r, --read`
Read results from existing log files without executing Gaussian.
- **Example**: `python mcal.py xxx.cif p -r`

#### `--resume`
Resume calculation using existing results if log files terminated normally.
- **Example**: `python mcal.py xxx.cif p --resume`

#### `--fullcal`
Disable speedup processing using moment of inertia and distance between centers of weight, and calculate transfer integrals for all pairs.
- **Example**: `python mcal.py xxx.cif p --fullcal`

#### `--cellsize <number>`
Specify the number of unit cells to expand in each direction around the central unit cell for transfer integral calculations.
- **Default**: `2` (creates 5×5×5 supercell)
- **Examples**: 
  - `python mcal.py xxx.cif p --cellsize 1` (creates 3×3×3 supercell)
  - `python mcal.py xxx.cif p --cellsize 3` (creates 7×7×7 supercell)

### Output Settings

#### `-p, --pickle`
Save calculation results to a pickle file.
- **Example**: `python mcal.py xxx.cif p -p`

### Diffusion Coefficient Calculation Methods

#### `--mc`
Calculate diffusion coefficient tensor using Monte Carlo method.
- **Example**: `python mcal.py xxx.cif p --mc`

#### `--pde`
Calculate diffusion coefficient tensor using Partial Differential Equation method.
- **Example**: `python mcal.py xxx.cif p --pde`

## Practical Usage Examples

### Basic Calculations
```bash
# Calculate mobility of p-type xxx
python mcal.py xxx.cif p

# Use 8 CPUs and 16GB memory
python mcal.py xxx.cif p -c 8 -m 16
```

### High-Precision Calculations
```bash
# Calculate transfer integrals for all pairs (high precision, time-consuming)
python mcal.py xxx.cif p --fullcal

# Use smaller supercell for faster calculation
python mcal.py xxx.cif p --cellsize 1

# Use larger supercell to widen transfer integral calculation range
python mcal.py xxx.cif p --cellsize 3

# Use different basis set
python mcal.py xxx.cif p -M "B3LYP/6-311G(d,p)"
```

### Reusing Results
```bash
# Read from existing calculation results
python mcal.py xxx.cif p -r

# Resume interrupted calculation
python mcal.py xxx.cif p --resume

# Save results to pickle file
python mcal.py xxx.cif p -p
```

### Comparing Diffusion Coefficients
```bash
# Compare with normal calculation + Monte Carlo + PDE methods
python mcal.py xxx.cif p --mc --pde
```

## Output

### Standard Output
- Reorganization energy
- Transfer integrals for each pair
- Diffusion coefficient tensor
- Mobility tensor
- Eigenvalues and eigenvectors of mobility

## Notes

1. **Calculation Time**: Calculation time varies significantly depending on the number of molecules and cell size
2. **Memory Usage**: Ensure sufficient memory for large systems
3. **Gaussian Installation**: Gaussian 09 or Gaussian 16 is required
4. **Dependencies**: Make sure all required Python libraries are installed

## Troubleshooting

### If calculation stops midway
```bash
# Resume with --resume option
python mcal.py xxx.cif p --resume
```

### Memory shortage error
```bash
# Increase memory amount
python mcal.py xxx.cif p -m 32
```

### To reduce calculation time
```bash
# Enable speedup processing (default)
python mcal.py xxx.cif p

# Increase number of CPUs
python mcal.py xxx.cif p -c 16
``` 
