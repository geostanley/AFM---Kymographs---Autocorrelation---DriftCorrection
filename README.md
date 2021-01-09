# Auto-correlation analysis of AFM kymographs (including drift-correction algorithm)

MATLAB scripts for the drift-correction and auto-correlation analysis of AFM kymographs.

The project is run from the **Auto_correlation_analysis_kymographs.m** script. In brief, it runs as follows:

1. Loads .spm AFM files from Bruker AFM systems, and concatonates the data chronologically.
2. Applies a drift correction algorithm (found in the **ImageFuncs** class) to the concatonated kymograph to ensure each pixel represents the same point in space with time.
3. Applies an auto-correlation functions to each pixel (defined in script).
4. Plots the results as a heatmap.

The results from kymographs recorded at the same scan speed can be loaded into **Average_kymographs_heatmaps.m** to render averaged results.

These scripts were used to for the analysis of stochastically fluctuating intrinsically disordered proteins inside native and synthetic nuclear pore complexes, for a study on nucleocytoplasmic transport. The results from this analysis can be found here: https://pubs.acs.org/doi/full/10.1021/acsnano.9b02424

## Authors

George J Stanley

## Acknowledgements

Acknowledgements to Jonathan Lanset for linspecer, and to John Iversen for freezeColors.

## DOI

[DOI Badge](https://zenodo.org/badge/latestdoi/177807176)
