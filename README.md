# ST-Cokriging ArcGIS extension for crime prediction


## Introduction
**1. Geostatistics, Cokriging, Crime Prediction**

In geostatistics, Cokriging is a multivariate variant of Kriging technique and makes spatial predictions for a sparsely sampled variable (the primary variable) of interest, with help of one or more well-sampled ancillary variables (the secondary co-variables). Cokriging method usually results in more accurate predictions of the target primary variable than Kriging. This is because Cokriging method exploits cross-correlations between the primary variable and the secondary co-variables in addition to the spatial autocorrelation of the primary variable.

Accurate crime prediction can help allocate police resources for crime reduction and prevention. There are two popular approaches to predict criminal activities: one is based on historical crime, and the other is based on environmental variables correlated with criminal patterns. Previous research on geo-statistical modeling mainly considered one type of data in space-time domain, and few sought to blend multi-source data. In this research, we proposed a spatio-temporal Cokriging algorithm to integrate historical crime data and urban transitional zones for more accurate crime prediction. Time-series historical crime data were used as the primary variable, while urban transitional zones identified from the VIIRS nightlight imagery were used as the secondary co-variable.

<img align="center" width="500" src="/Images/fg1.png">

Figure 1. Data processing flowchart of ST-Cokriging method 

**2.	ST-Cokriging work flow**

In ST-cokriging formulation, we assume that the primary variable of interest is coarse spatial resolution Images that are sampled at a high temporal frequency (high temporal resolution), and the secondary variable (co-variable) are ﬁne spatial resolution Images that are sparsely sampled over time (low temporal resolution), as shown in Figure 1. Without loss of generality, we only consider the case with only one co-variable observed at multiple time points in the mathematical formulation of ST-Cokriging method. The extension to two or more co-variables observed at multiple time points is straightforward.

**3.	Extension implementation**
In the Cokriging linear system, the covariance matrix C is of very large size. Since the co-variable is observed at high spatial resolution, it is very likely that covariance matrix are very large. Similarly, even with smaller dimension for the primary variable observed at coarse spatial resolution, adding temporal dimension can still be large, since the primary variable is observed at high temporal frequency. Therefore, the matrix C can be of high dimension. Solving such a high-dimensional linear system can be computationally infeasible. One popular method to alleviate this diﬃculty is to force small numbers in the matrix (vector) to be zero, known as thresholding or tapering. Meanwhile, we take advantage of the feature of regularly gridded data in the application presented in Section 3, which facilitates eﬃcient parallel computing of the Cokriging predictor and variance.

## Source and Environment 

1.	Source code download
	Download the toolbox (yangtoolbox_crime.tbx) for ArcGIS from this repository folder (Toolbox and Codes). 
	Download python scripts for ST-Cokriging from this repository folder (Toolbox and Codes)
	 * **FittingVariog_crime.py** to be linked with toolbox **FittingVariogram**
	 * **ImageFusion_Speedy_crime.py** to be linked with toolbox **STCoKriging_crime**
	 * **SemiVariog.py** to be linked with toolbox **Variogram**
	
2.	Environmental setup
	- Enable extensions in ArcMap
	
	<img width="300" src="/Images/fg2.png">
	
	- Setup the Geo-processing option
	
	<img width="300" src="/Images/fg3.png">
	
	- Open Arctoolbox window, right click and add a toolbox
	
	<img width="300" src="/Images/fg4.png">
	
	- Navigate to the toolbox just downloaded and select
	
	<img width="300" src="/Images/fg5.png">
	
	- Unfold the toolbox and the scripts should be appeared in the toolbox
	
	<img width="300" src="/Images/fg6.png">
	
	- Right click the ST-Cokriging script and select properties
	
	<img width="300" src="/Images/fg7.png">
	
	- Click source tab and link the ST-Cokriging script to the toolbox
	- Do the same procedure to the semi-variogram script and FittingVariog_crime script link the downloaded code to the arctoolbox script.
	
## Prediction via toolbox
* Step I, using the spatio-temporal data to calculate spatial and temporal semi-variogram
* Estimate the spatio-temporal semi-variograms
* Input the parameters as shown in the figure below. All 13 quad-week images should be included in the calculation and arranged in chronological order. 
* The input spatial raster should be one of the quad-week period with largest number of crime. 
* The spatial sample ratio is the percentage of the subset of spatial samples, depends on the resolution and total pixels of the study ration, the subset population should be around 3,000 – 10, 000. 
* At last, select the path to store the txt files of spatial and temporal semi-variogram.

<img width="300" src="/Images/fg8.png">

* Step II Fitting spatial and temporal semi-variogram
* the output are txt files of spatial and temporal semi-variogram
* input the spatial and temporal semi-variogram files to be fitting to function
* Choose the fitting function based on the shape of the semi-varigram. In the example below, the spatial semivariogram is fit for Gaussian function, while the temporal semi-variogram is fit for Exponential fucntion
* The fitted output are two txt files, using parameters of nuggest, sill, and range to depict the spatial and temporal dependance

<img width="600" src="/Images/fg11.png">

* The spatial and temporal semivariogram will be converted to covariance/correlation, then combined to spatio-temporal covariance fucntion

<img width="600" src="/Images/fg12.png">

Step III, predict using ST-Cokriging
* In put the secondary covariable image, and time-series primary variable, the primary variable should be input in the time-series order. 
* The number of time-series should be no less than 3 for using of spatio-temporal prediction. 
* Input the fitted fucntion of spatial and temporal semi-vatiogram from Step II, and assign other parameters. 
<img width="300" src="/Images/fg10.png">

## Troubleshooting
1. Both time-series primary variable and co-variable should be re-project to same coordinates system as well as same datum.
2. Co-variable should be at finer spatial resolution or same. In this version, the coverage of the both co-variable and time-series primary variable should be same, and at same spatial resolution.
3. If using 32-bit OS software computer, the maximum number of row and column are 6000 by 6000.
4. The fusion process may take a while (about 1 hour depends on the CPU and RAM of the computer). Once it started running, please don’t use other program at same time. Otherwise window 10 might report ArcGIS no response. When no response appeared, please choose to wait rather than kill the process.
5. The final product will have 1-pixel-width edge at constant value (or very smooth). It is normal because current version did not model the edge effects, and leave the edge for showing the trend value.

## Citations
* Yang, B., Liu, L., Lan, M., Wang, Z., Zhou, H., Yu, H., (2020). A spatio-temporal method for crime prediction using historical crime data and transitional zones identified from nightlight imagery. *International Journal of Geographical Information Science*, 1-25. [DOI: 10.1080/13658816.2020.1737701](https://doi.org/10.1080/13658816.2020.1737701)








