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

## Source and Environmental 

1.	Source code download
	- Download the toolbox in ArcGIS from here.
	- Download the script for spatio-temporal semi-variogram here.
	- Download the script for ST-Cokriging here.
2.	Environmental setup
	- Enable extensions in ArcMap
	
	<img src="/Images/fg2.png">
	
	- Setup the Geo-processing option
	
	<img src="/Images/fg3.png">
	
	- Open Arctoolbox window, right click and add a toolbox
	
	<img src="/Images/fg4.png">
	
	- Navigate to the toolbox just downloaded and select
	
	<img src="/Images/fg5.png">
	
	- Unfold the toolbox and the scripts should be appeared in the toolbox
	
	<img src="/Images/fg6.png">
	
	- Right click the ST-Cokriging script and select properties
	
	<img src="/Images/fg7.png">
	
	- Click source tab and link the ST-Cokriging script to the toolbox
	- Do the same procedure to the semi-variogram script and link the downloaded code to the arctoolbox script.
	
## Troubleshooting
1. Both time-series primary variable and co-variable should be re-project to same coordinates system as well as same datum.
2. Co-variable should be at finer spatial resolution or same. In this version, the coverage of the both co-variable and time-series primary variable should be same, and at same spatial resolution.
3. If using 32-bit OS software computer, the maximum number of row and column are 6000 by 6000.
4. The fusion process may take a while (about 1 hour depends on the CPU and RAM of the computer). Once it started running, please don’t use other program at same time. Otherwise window 10 might report ArcGIS no response. When no response appeared, please choose to wait rather than kill the process.
5. The final product will have 1-pixel-width edge at constant value (or very smooth). It is normal because current version did not model the edge effects, and leave the edge for showing the trend value.

## Citations
* Yang, B., Liu, L., Lan, M., Wang, Z., Zhou, H., Yu, H., (2020). A spatio-temporal method for crime prediction using historical crime data and transitional zones identified from nightlight imagery. *International Journal of Geographical Information Science*, 1-25. [DOI: 10.1080/13658816.2020.1737701](https://doi.org/10.1080/13658816.2020.1737701)
* Yu, H.; Liu, L.; **Yang, B.;** Lan, M. Crime Prediction with Historical Crime and Movement Data of Potential Offenders Using a Spatio-Temporal Cokriging Method. *ISPRS International Journal of Geo-Information*. 2020, 9, 732. [DOI:10.3390/ijgi9120732](https://doi.org/10.3390/ijgi9120732)
* Zhou, H., Liu, L., Lan, M., **Yang, B.**, Wang, Z. (2019) Assessing the Impact of Nightlight Gradients on Street Robbery and Burglary in Cincinnati of Ohio State, USA. *Remote Sensing*. 11, 1958. [DOI: 10.3390/rs11171958](https://doi.org/10.3390/rs11171958)






