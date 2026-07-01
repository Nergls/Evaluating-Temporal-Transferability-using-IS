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

The notebook reproduces the statistical analyses, temporal transferability evaluation, model summaries, and figures presented in the manuscript.

---


## 📁 Repository Contents

###### **Core Files**

| File | Description |
|------|-------------|
| `027_PPR2_stats_published.ipynb` | Complete analysis workflow used in the manuscript |
| `geopandas_env.yml` | Conda environment for reproducing the Python environment |



## 📦 Dependencies

###### **🐍 Python Environment (for spatial preprocessing notebooks)**

The Python Jupyter Notebooks in this repository (files 020\_\*.ipynb to 024\_\*.ipynb) require a dedicated geospatial environment.

A reproducible Conda environment file is provided:

* geopandas\_env.yml

**Create the environment:**

conda env create -f geopandas\_env.yml

conda activate geopandas\_env



## ✉️ Support

For questions or suggestions, feel free to open an issue or contact the repository maintainer.

