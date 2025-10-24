# Example analysis workflow for RBC dataset
This repository contains scripts and supporting data for analyzing functional and structural neuroimaging data in relation to development and psychopathology, as part of the Reproducible Brain Charts (RBC) initiative and in support of the preprint:

Shafiei G, Esper NB, et al., (2025). Reproducible Brain Charts: An open data resource for mapping brain development and its associations with mental health. Neuron. DOI: [10.1016/j.neuron.2025.08.026](https://doi.org/10.1016/j.neuron.2025.08.026)

Please cite appropriately and consult [https://reprobrainchart.github.io/docs/get_data](https://reprobrainchart.github.io/docs/get_data) to download and prepare the necessary datasets before using this repository.


## `code`
All analysis and data preparation code is located in the [`code/`](code/) folder. Code is written in Python and R.

### Functional Data Workflow (Python + R):
1. **Obtaining Data**: Follow instructions in [GetData](https://reprobrainchart.github.io/docs/get_data) to clone C-PAC files. Two versions are used:
   - `warning-fail` (labeled as `noqc` within scripts and filenames): Includes all scans, even failed ones.
   - `complete-artifact` (labeled as `artifact` within scripts and filenames): Excludes scans marked as "Fail".

2. **Within- and Between-Network FC Calculation**:
   - `scpt_rbc_withinbetween_rsnFC.py`: Computes within- and between-network functional connectivity matrices for 7 Yeo-Kieren resting-state functional networks .

3. **Data Preparation for R Analyses**:
   - `scpt_prepare_for_r_fc_rsn.py`: Prepares and saves the data for later use in R-based Generalized Additive Model (GAM) analysis.

4. **Modeling and Visualization (R)**:
   - `scpt_rbc_combined_fc_rsn7.R`: Performs GAM analyses on functional connectivity data (combined and study-specific).
   - `func_GAM_rbc.R`: Provides GAM fitting functions.
   - `covbat.R`: Harmonizes functional data across sites before modeling.

5. **Additional Visualization (Python)**:
   - `scpt_rbc_rsn7x7_visualize.py`: Generates a subset of figures based on functional data. The rest of the figures are generated within R scripts (i.e., `scpt_rbc_combined_fc_rsn7.R`).

### Structural Data Workflow (Python + R):
1. **Obtaining Data**: Again, follow [GetData](https://reprobrainchart.github.io/docs/get_data) to clone FreeSurfer files. Two versions used:
   - `warning-fail` (labeled as `noqc` within scripts and filenames): Includes all scans, even failed ones.
   - `complete-artifact` (labeled as `artifact` within scripts and filenames): Excludes scans marked as "Fail".

2. **Data Preparation for R Analyses**:
   - `scpt_prepare_for_r_fstabulate.py`: Prepares FreeSurfer data files for downstream R-based analysis.

3. **Modeling and Visualization (R)**:
   - `scpt_rbc_combined_fstabulate.R`: Performs GAM analyses on structural metrics across parcels.
   - `func_GAM_rbc.R`: Provides GAM fitting functions.
   - `covbat.R`: Harmonizes structural data across sites prior to modeling.

### R Scripts:
All R scripts for statistical modeling and harmonization are located in the `code/` folder.

- `scpt_rbc_combined_fc_rsn7.R`: Performs GAM analyses on resting-state functional connectivity network pairs using either age or p-factor as predictors. Includes support for combined and study-specific datasets, and generates both summary statistics and visualizations (e.g., ribbons, smooth curves).
  
- `scpt_rbc_combined_fstabulate.R`: Conducts GAM analyses on structural imaging metrics (e.g., cortical thickness, surface area, gyrification index) across Schaefer-400 parcels. Outputs include region-wise statistics, visualizations (brain maps, histograms), and study-specific comparison plots.

- `func_GAM_rbc.R`: Contains helper functions for GAM fitting and prediction, used by the above scripts. Supports smooth and linear modeling for both functional and structural data.

- `covbat.R`: Harmonizes structural and functional neuroimaging data across study sites using the `covfam()` function from the `ComBatFamily` package. Supports inclusion of covariates such as age, sex, motion (FD), Euler number, and p-factor, and outputs harmonized matrices for use in downstream GAM analyses.

### Python Scripts:
All Python scripts used for data preparation and preliminary processing are located in the `code/` folder.

- `scpt_rbc_withinbetween_rsnFC.py`: Calculates within- and between-network functional connectivity matrices based on C-PAC outputs for the 7 Yeo-Kieren resting-state networks. Used as the initial step in the functional analysis pipeline.

- `scpt_prepare_for_r_fc_rsn.py`: Formats and filters the functional connectivity data (output of the previous script) for downstream analysis in R. Prepares datasets aligned with developmental and psychopathology models.

- `scpt_prepare_for_r_fstabulate.py`: Prepares structural imaging data (e.g., FreeSurfer outputs) for analysis in R by tabulating relevant metrics per subject and parcel, formatted for GAM input.

- `scpt_rbc_rsn7x7_visualize.py`: Generates visual representations of functional connectivity matrices across the 7-network parcellation. Includes heatmaps and QC-friendly visual summaries of network relationships.

## `data`
The [`data/dataR`](data/dataR/) folder contains all datasets prepared for R-based analyses. These include:
- **Combined files**: Primary data files that aggregate all studies.
- **Study-specific files**: For study-specific analyses.
