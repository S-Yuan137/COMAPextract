Spectral Process
================

Install & required Packages
---------------------------
Requires:
- astropy
- healpy
- h5py
- matplotlib
- numpy
- scipy
- tqdm

Running
-------
1. Prepare a parameter file (following example given)
1. Run `python main.py input_file.py`

Outputs
-------
- **Source:_{SourceName}.hd5** HDF5 file containing all source fits and errors for each fit. Able to hold one set of fits for each method.
- **{SourceName}_PhotometryReport_{Method}.pdf** PDF containing all summary plots for each fit alongside a summary plot of all fluxes at all frequencies fitted.
