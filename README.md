# Example analysis workflow for RBC dataset
This repository contains scripts and supporting data for analyzing functional and structural neuroimaging data in relation to development and psychopathology, as part of the ReproBrainChart initiative and in support of the preprint:

Shafiei G, Esper NB, et al., (2022). Reproducible Brain Charts: An open data resource for mapping brain development and its associations with mental health. bioRxiv.

Please cite appropriately and consult [https://reprobrainchart.github.io/docs/get_data](https://reprobrainchart.github.io/docs/get_data) to download and prepare the necessary datasets before using this repository.

## `code`
The [code](code/) folder contains all the code used to run the analyses and generate the figures.

## `data`
The [data](data/) folder contains the data used to run the analyses. 



## `code`
All analysis and data preparation code is located in the [`code/`](code/) folder. Code is written in Python and R.

### Functional Data Workflow (Python):
1. **Obtaining Data**: Follow instructions in [GetData](https://reprobrainchart.github.io/docs/get_data) to clone C-PAC files. Two versions are used:
   - `warning-fail` (labeled as `noqc` within scripts and filenames): Includes all scans, even failed ones.
   - `complete-artifact` (labeled as `artifact` within scripts and filenames): Excludes scans marked as "Fail".

2. **Within- and Between-Network FC Calculation**:
   - `scpt_rbc_withinbetween_rsnFC.py`: Computes within- and between-network functional connectivity matrices for 7 Yeo-Kieren resting-state functional networks .

3. **Data Preparation for R Analyses**:
   - `scpt_prepare_for_r_fc_rsn.py`: Prepares and saves the data for later use in R-based Generalized Additive Model (GAM) analysis.

4. **Visualization**:
   - `scpt_rbc_rsn7x7_visualize.py`: Generates a subset of figures based on functional data. The rest of the figures are generated within R scripts (see below).

### Structural Data Workflow (Python):
1. **Obtaining Data**: Again, follow [GetData](https://reprobrainchart.github.io/docs/get_data) to clone FreeSurfer files. Two versions used:
   - `warning-fail` (labeled as `noqc` within scripts and filenames): Includes all scans, even failed ones.
   - `complete-artifact` (labeled as `artifact` within scripts and filenames): Excludes scans marked as "Fail".

2. **Data Preparation**:
   - `scpt_prepare_for_r_fstabulate.py`: Prepares FreeSurfer data files for downstream R-based analysis.

### R Scripts:
All R scripts for statistical modeling and harmonization are located in the `code/` folder.

- `scpt_rbc_combined_fc_rsn7.R`: GAM analysis on functional connectivity data.
- `scpt_rbc_combined_fstabulate.R`: GAM analysis on structural data.
- `func_GAM_rbc.R`: Contains all R functions used in GAM analyses.
- `covbat.R`: Used for neuroimaging data harmonization per modality and analysis type (i.e., developmental and psychopathology analyses).

## `data`
The [`data/dataR`](data/dataR/) folder contains all datasets prepared for R-based analyses. These include:
- **Combined files**: Primary data files that aggregate all studies.
- **Study-specific files**: For study-specific analyses.
