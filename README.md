# Soil Desiccation and Land Subsidence Analysis in Europe — INSA Rouen Normandie x Cerema

A satellite data processing and statistical analysis study in **R** investigating the coupling between soil desiccation / total water storage anomalies (GRACE) and vertical ground motion / land subsidence (EGMS InSAR).

This project was carried out as part of the GM4 (Mathematical Engineering) curriculum at **INSA Rouen Normandie** (Academic Year 2025-2026), in collaboration with **Cerema**, supervised by **Mr. Lionel Riou-Durand** (INSA Rouen) and **Dr. Teodolina Lopez** (Cerema), and created by **Cheikh Seck**, **Mouad Sheradj Drissi**, and **Amr Zaki Salih**.

---

## Overview

Climate change and groundwater overdraft accelerate ground deformation phenomena across Europe. This project explores the physical and statistical relationship between:
1. **Total Water Storage (TWS) anomalies** measured by the gravimetric satellite mission **GRACE / GRACE-FO (G3P v1.12)**.
2. **Surface displacements** recorded at millimeter accuracy by the **European Ground Motion Service (EGMS)** using Sentinel-1 Synthetic Aperture Radar (InSAR) interferometry.

The case study focuses on the urban and agricultural area of **Caen, Normandy, France** ($49.18^\circ\text{N}, -0.37^\circ\text{E}$).

---

## Datasets

- **GRACE (G3P v1.12):** Global NetCDF gridded data with a spatial resolution of $0.5^\circ \times 0.5^\circ$ (~50 km) providing monthly equivalent water thickness (TWS) anomalies.
- **EGMS (European Ground Motion Service):** High-resolution InSAR displacement time series (CSV, 100 km grid cells, 2019–2023) recording line-of-sight and vertical ground displacements.

---

## Methodologies

Two spatial aggregation and filtering strategies were developed in **R**:

1. **GRACE Pixel Selection (`projet_S7_MSRO_MethodeCarreGrace.r`):**
   - Harmonizes coordinate reference systems (CRS: WGS84 / EPSG:4326 vs UTM/local projections).
   - Identifies the specific GRACE pixel bounding box encompassing Caen ($[49.0^\circ, 49.5^\circ]\text{N} \times [-0.5^\circ, 0.0^\circ]\text{E}$).
   - Filters all EGMS radar scatterers within this bounding box and computes temporal averages for correlation analysis.

2. **Euclidean Radial Buffer (`projet_S7_MSRO_Methode_distance_euclidienne.r`):**
   - Applies an incremental Euclidean distance radius around the coordinates of Caen.
   - Filters EGMS measurements progressively to observe spatial variance, damping of noise, and localized ground motion dynamics.
   - Generates detailed interactive map visualizations (`leaflet`) and comparative time-series plots.

---

## Project Structure

```
Projet_S7_MSRO/
├── projet_S7_MSRO_MethodeCarreGrace.r        # Spatial filtering by GRACE grid pixel
├── projet_S7_MSRO_Methode_distance_euclidienne.r # Progressive Euclidean radius filtering & mapping
├── Rapport_projetMSRO_SHERADJDRISSIMouad_SECKCheikh_ZAKISALIHAmr.pdf # Full academic project report
├── Rapport_CRS_PROJET_MSRO.pdf              # Technical note on Coordinate Reference Systems (CRS)
├── .gitignore                               # Git ignore rules for R environment
└── README.md
```

---

## Prerequisites & Installation

The project requires **R (version 4.0 or higher)** and the following scientific and geospatial libraries:

```r
install.packages(c("ncdf4", "CFtime", "sf", "leaflet", "lubridate"))
```

---

## Running the Analyses

1. Clone or download the repository.
2. Update the file paths pointing to your local GRACE NetCDF file and EGMS CSV datasets in the scripts.
3. Run the desired analysis script:
   - For the GRACE pixel boundary approach:
     ```bash
     Rscript projet_S7_MSRO_MethodeCarreGrace.r
     ```
   - For the Euclidean buffer and visualization approach:
     ```bash
     Rscript projet_S7_MSRO_Methode_distance_euclidienne.r
     ```

---

## Reports & Documentation

The repository includes the full project deliverables:
- **`Rapport_projetMSRO_SHERADJDRISSIMouad_SECKCheikh_ZAKISALIHAmr.pdf`**: Complete report covering theoretical background, CRS reprojections, temporal harmonization, data analysis, and correlation findings.
- **`Rapport_CRS_PROJET_MSRO.pdf`**: Methodological note dedicated to geographic coordinate transformations and spatial alignment.

---

## Authors & Acknowledgments

- **Cheikh Seck**
- **Mouad Sheradj Drissi**
- **Amr Zaki Salih**

### Project Advisors:
- **Dr. Teodolina Lopez** — Cerema (Centre d'études et d'expertise sur les risques, l'environnement, la mobilité et l'aménagement) / INSA Rouen Normandie
- **Mr. Lionel Riou-Durand** — INSA Rouen Normandie

---

## Institutions & Partners

- **INSA Rouen Normandie** — Department of Mathematical Engineering
- **Cerema** (End-user partnership on geotechnical risks and environmental monitoring)
