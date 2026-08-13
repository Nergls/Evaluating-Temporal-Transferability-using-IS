# Evaluating Temporal Transferability of Imaging Spectroscopy-Based Functional Trait Models in Alpine Grasslands

This repository accompanies the manuscript:

> *Evaluating Temporal Transferability of Imaging Spectroscopy-Based Functional Trait Models in Alpine Grasslands*

It provides the analysis workflow for evaluating the temporal transferability of imaging spectroscopy models to estimate multiple plant functional traits across independent years.

The workflow integrates spectral reflectance data, field measurements, and Partial Least Squares Regression (PLSR) to assess model performance within calibration years and during transfer to independent (prediction) years.

The evaluated traits include:

- Fresh above-ground biomass (Fresh AGB; g m⁻²)
- Dry above-ground biomass (Dry AGB; g m⁻²)
- Carbon (C%)
- Nitrogen (N%)
- Acid detergent lignin (ADL%)
- Acid detergent fiber (ADF%)
- Neutral detergent fiber (NDF%)

The repository provides the PLSR modeling and statistical analysis workflows used in the manuscript.

---


## 📁 Repository Contents

##### **Core Files**

| File | Description |
|------|-------------|
| `012_PPR2_PLSR_PRO-up2d8_debugged_published.R` | PLSR model calibration and within-year validation workflow |
| `012b_PPR2_PLSR__PRO_Revision_Updates_published.R` | Revised PLSR model-performance calculations |
| `012c_PPR2_PLSR__PRO_Revision_Updates_cross_year_published.R` | Cross-year prediction and revised performance evaluation |
| `018_CALL_cumulative_models_FUNCTION_published.R` | Workflow for running pooled multi-year calibration models |
| `027_PPR2_stats_organizer_revised_published.ipynb` | Statistical analyses, model summaries, and output organization used in the manuscript |
| `Randomly_generated_sample_data_structure.csv` | Synthetic dataset reproducing the input-data structure required to demonstrate the modeling workflow |
| `geopandas_env.yml` | Conda environment for reproducing the Python environment |



## 📦 Dependencies

##### **🐍 Python Environment **

The Python statistical-analysis workflow in 027_PPR2_stats_organizer_revised_published.ipynb requires the Python environment specified in geopandas\_env.yml.

**Create the environment:**

conda env create -f geopandas\_env.yml

conda activate geopandas\_env



## ✉️ Support

For questions or suggestions, feel free to open an issue or contact the repository maintainer.

