# ST-Cokriging ArcGIS extension for crime prediction


## Introduction
**1. Geostatistics and Cokriging**

In geostatistics, Cokriging is a multivariate variant of Kriging technique and makes spatial predictions for a sparsely sampled variable (the primary variable) of interest, with help of one or more well-sampled ancillary variables (the secondary co-variables). Cokriging method usually results in more accurate predictions of the target primary variable than Kriging. This is because Cokriging method exploits cross-correlations between the primary variable and the secondary co-variables in addition to the spatial autocorrelation of the primary variable.
In conventional Cokriging method, both the primary variable and co-variables are in spatial domain, and time dimension is not taken into consideration. By extending it from sole spatial domain to the spatio-temporal domain, this algorithm formulated a ST-Cokriging method that is capable of taking advantage of both spatial and temporal correlation within and between primary and co-variables to produce temporally frequent predictions for the primary variable at a high spatial resolution as the co-variables.

![image](/Images/fg1.png){:height="220px"}

Figure 1. Data processing flowchart of ST-Cokriging method 

**2.	ST-Cokriging work flow**

In ST-cokriging formulation, we assume that the primary variable of interest is coarse spatial resolution Images that are sampled at a high temporal frequency (high temporal resolution), and the secondary variable (co-variable) are ﬁne spatial resolution Images that are sparsely sampled over time (low temporal resolution), as shown in Figure 1. Without loss of generality, we only consider the case with only one co-variable observed at multiple time points in the mathematical formulation of ST-Cokriging method. The extension to two or more co-variables observed at multiple time points is straightforward.

**3.	Extension implementation**
In the Cokriging linear system, the covariance matrix C is of very large size. Since the co-variable is observed at high spatial resolution, it is very likely that M_j’s are very large. Similarly, even with smaller N_i’s for the primary variable observed at coarse spatial resolution, ∑_(i=1)^T▒N_i  can still be large, since the primary variable is observed at high temporal frequency. Therefore, the matrix C can be of high dimension. Solving such a high-dimensional linear system can be computationally infeasible. One popular method to alleviate this diﬃculty is to force small numbers in the matrix (vector) to be zero, known as thresholding or tapering. Meanwhile, we take advantage of the feature of regularly gridded data in the application presented in Section 3, which facilitates eﬃcient parallel computing of the Cokriging predictor and variance.

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


