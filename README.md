# ST-Cokriging: An ArcGIS Toolbox for Spatio-Temporal Crime Prediction

An ArcGIS toolbox and Python implementation of a spatio-temporal Cokriging (ST-Cokriging) method for predicting fine-resolution crime patterns by fusing historical crime records with auxiliary remote sensing data.

## Citation

If you use this toolbox in your research, please cite:

> Yang, B., Liu, L., Lan, M., Wang, Z., Zhou, H., & Yu, H. (2020). A spatio-temporal method for crime prediction using historical crime data and transitional zones identified from nightlight imagery. *International Journal of Geographical Information Science*, 1–25. [DOI: 10.1080/13658816.2020.1737701](https://doi.org/10.1080/13658816.2020.1737701)

## Background

Accurate crime prediction helps law enforcement agencies allocate resources more effectively, identify emerging hotspots, and design targeted prevention strategies. Two broad approaches dominate the literature: methods that extrapolate from **historical crime records**, and methods that relate crime to **environmental or socioeconomic covariates** known to correlate with criminal activity.

Most geostatistical approaches to crime modeling consider only one of these data types in the space-time domain, with relatively few attempts to blend multiple data sources. This project implements a spatio-temporal Cokriging algorithm that integrates both: time-series historical crime data serve as the primary variable, while urban transitional zones derived from VIIRS nightlight imagery serve as a secondary co-variable. The result is higher-resolution, more accurate crime prediction than either data source could produce alone.

## Goal

Use geostatistical fusion methods to predict crime occurrence and hotspot locations at finer spatial resolution than the primary (coarse) crime data alone would allow, by leveraging a finer-resolution, more frequently available auxiliary data source.

## Data

- **Crime data (primary variable):** time-series incident records including crime type, location, and time. Publicly available examples include city/police department open-data portals, such as the [San Francisco Police Department Incident Reports](https://data.sfgov.org/Public-Safety/Police-Department-Incident-Reports-2018-to-Presen/wg3w-h783).
- **Auxiliary data (secondary co-variable):** demographic, economic, or remote sensing data correlated with crime occurrence — in the published study, urban transitional zones identified from VIIRS nightlight imagery. Socioeconomic covariates can be sourced from the [U.S. Census Bureau](https://www.census.gov/data.html).

## Method: ST-Cokriging Workflow

<img align="center" width="500" src="/Images/fg1.png">

*Figure 1. Data processing flowchart of the ST-Cokriging method.*

In the ST-Cokriging formulation, the primary variable is assumed to be a set of coarse-spatial-resolution images sampled at high temporal frequency, while the secondary co-variable consists of fine-spatial-resolution images sparsely sampled over time (Figure 1). For clarity, the mathematical formulation considers a single co-variable observed at multiple time points; the extension to two or more co-variables is straightforward.

**Computational considerations.** In the Cokriging linear system, the covariance matrix *C* can become very large: the co-variable's fine spatial resolution alone can drive up its dimension, and even a coarser-resolution primary variable adds a large temporal dimension from its high sampling frequency. Solving such a high-dimensional linear system directly can be computationally infeasible. A common remedy is to force small matrix entries to zero (thresholding/tapering). This implementation also takes advantage of the regularly gridded structure of the input data to enable efficient parallel computation of the Cokriging predictor and its variance.

## Installation and Environment Setup

1. **Download the source code**
   - Download the ArcGIS toolbox (`yangtoolbox_crime.tbx`) from the [`Toolbox and Codes`](./Toolbox%20and%20Codes) folder in this repository.
   - Download the accompanying Python scripts from the same folder:
     - `FittingVariog_crime.py` — links to the **FittingVariogram** tool
     - `ImageFusion_Speedy_crime.py` — links to the **STCoKriging_crime** tool
     - `SemiVariog.py` — links to the **Variogram** tool

2. **Set up the ArcMap environment**

   Enable the required extensions in ArcMap:

   <img width="300" src="/Images/fg2.png">

   Configure the geoprocessing options:

   <img width="300" src="/Images/fg3.png">

   Open the ArcToolbox window, right-click, and add a toolbox:

   <img width="300" src="/Images/fg4.png">

   Navigate to the downloaded toolbox and select it:

   <img width="300" src="/Images/fg5.png">

   Expand the toolbox — the three scripts should now appear:

   <img width="300" src="/Images/fg6.png">

   Right-click each script and open its properties:

   <img width="300" src="/Images/fg7.png">

   On the **Source** tab, link each script to its corresponding downloaded `.py` file (`ST-Cokriging` → `ImageFusion_Speedy_crime.py`, `Variogram` → `SemiVariog.py`, `FittingVariogram` → `FittingVariog_crime.py`).

## Running a Prediction

**Step I — Compute spatial and temporal semi-variograms**

- Estimate the spatio-temporal semi-variograms using the input parameters shown below. All quad-week images (13 in the published study) should be included and arranged in chronological order.
- The input spatial raster should correspond to the quad-week period with the highest crime count.
- The spatial sample ratio controls the subset of spatial samples used; depending on the study area's resolution and total pixel count, a subset of roughly 3,000–10,000 samples is recommended.
- Select an output path for the spatial and temporal semi-variogram text files.

<img width="300" src="/Images/fg8.png">

**Step II — Fit the spatial and temporal semi-variograms**

- The outputs of Step I are text files containing the empirical spatial and temporal semi-variograms.
- Supply these files and choose a fitting function based on the shape of each semi-variogram. In the example below, the spatial semi-variogram is fit with a Gaussian function and the temporal semi-variogram with an Exponential function.
- The fitted output consists of two text files describing the spatial and temporal dependence structure via nugget, sill, and range parameters.

<img width="600" src="/Images/fg11.png">

The fitted spatial and temporal semi-variograms are converted to covariance/correlation and combined into a joint spatio-temporal covariance function:

<img width="600" src="/Images/fg12.png">

**Step III — Predict using ST-Cokriging**

- Input the secondary co-variable image and the time-series primary variable (in chronological order). At least 3 time points are required for spatio-temporal prediction.
- Input the fitted spatial and temporal semi-variogram functions from Step II, along with the remaining parameters.

<img width="300" src="/Images/fg10.png">

## Applications

The predictions produced by this workflow can help law enforcement agencies allocate resources more effectively, identify potential crime hotspots at finer spatial resolution than raw crime data alone permits, and inform crime prevention strategy.

## Troubleshooting

1. Both the time-series primary variable and the co-variable must be reprojected to the same coordinate system and datum.
2. The co-variable must be at the same or finer spatial resolution than the primary variable. In this version, both the co-variable and the time-series primary variable must share the same spatial coverage and resolution.
3. On a 32-bit OS/ArcGIS installation, the maximum raster size is 6,000 × 6,000 cells.
4. The fusion process can take roughly an hour, depending on CPU and RAM. Once started, avoid running other programs concurrently — Windows may report "ArcGIS Not Responding" under load. If this occurs, continue waiting rather than terminating the process.
5. The final output has a 1-pixel-wide edge at a near-constant value. This is expected: the current version does not model edge effects and instead leaves the edge at the trend value.
