'''
Module Name: FittingVariog

Created on Apr 21, 2014

@author: Bo Yang

All Rights Reserved
'''


#===================================================================================
#               Set Module Parameter
#===================================================================================
import arcpy
import math
import time
from numpy.core.numeric import arange

NODATA= -3.40282346639e+038
spatialVarigramTxt= arcpy.gp.GetParameterAsText(0)
temporalVarigramTxt= arcpy.gp.GetParameterAsText(1)
spatialFittingRange= int(arcpy.gp.GetParameterAsText(2))
temporalFittingRange= int(arcpy.gp.GetParameterAsText(3))
fittingModel= str(arcpy.gp.GetParameterAsText(4))
saveFittingParameter= arcpy.gp.GetParameterAsText(5)

#=====================================================================
#                define fitting model
#======================================================================

def Gaussian(t, s, p,x):
    result=(t*1.0)+1.0*s*(1-math.exp(-(x**2)/((p*1.0)**2)))
    return result

def Exponential(t,s,p,x):
    result=(t*1.0)+((s*1.0))*(1-math.exp(-x/(p*1.0)))
    return result

#=========================================================================
#                read spatial and temporal data
#=========================================================================
TxtSpatialVariogram = open(spatialVarigramTxt)
spatialDataInput=TxtSpatialVariogram.readline()
TxtSpatialVariogram.close()
spatialDataString=spatialDataInput.split(',')
spatialDataList = [float(item) for item in spatialDataString]

txtTemporalVariogram = open(temporalVarigramTxt)
temporalDataInput=txtTemporalVariogram.readline()
txtTemporalVariogram.close()
temporalDataString=temporalDataInput.split(',')
temporalDataList = [float(item) for item in temporalDataString]

#===================================================================================
#                    Fitting ST data 
#===================================================================================

def ST_gau(data,h):
    arcpy.AddMessage(str(time.strftime("%m-%d %X",time.localtime())) + "; Start Gaussian Fitting Process")
    Min = 100000
    tm,sm,pm=0,0,0
    for t in arange(0,0.0001,0.001):
        for s in arange(0.000000000000001,0.0000000001,0.000000000000001):
            for p in arange(0.1,4,0.01):
                result=0
                for x in range(1,h):
                    result=result+(Gaussian(t,s,p,x)-data[x])**2
                if result<Min:
                    Min =result
                    tm=t
                    sm=s
                    pm=p
    fittingGaussianResult = str(tm) + ";" + str(sm) + ";"+ str(pm)+ ";" + str(result)
    return fittingGaussianResult
    
def ST_exp(data,h):
    arcpy.AddMessage(str(time.strftime("%m-%d %X",time.localtime())) + "; Start Exponential Fitting Process")
    Min = 100000
    tm,sm,pm=0,0,0
    for t in arange(0,0.0001,0.001):
        for s in arange(0.000000000000001,0.0000000001,0.000000000000001):
            for p in arange(0.1,4,0.01):
                result=0
                for x in range(1,h):
                    result=result+(Exponential(t,s,p,x)-data[x])**2
                if result<Min:
                    Min =result
                    tm=t
                    sm=s
                    pm=p
    fittingExponentialResult = str(tm) + ";" + str(sm) + ";" + str(pm) + ";" + str(result)
    return fittingExponentialResult

#=============================================================================================
#                Run the Main Program
#=============================================================================================
if fittingModel=="Gaussian":
    FittingParameterTxt = open(saveFittingParameter, 'w')
    FittingParameterTxt.write("Gaussian Spaial fitting result:\n")
    FittingParameterTxt.write(ST_gau(spatialDataList,spatialFittingRange))
    FittingParameterTxt.write("\n")
    FittingParameterTxt.write("Gaussian Temporal fitting result:\n")
    FittingParameterTxt.write(ST_gau(temporalDataList,temporalFittingRange))
    FittingParameterTxt.close()
else:
    fittingModel=="Exponential"
    FittingParameterTxt = open(saveFittingParameter, 'w')
    FittingParameterTxt.write("Exponential Spaial fitting result: Nuggest; Sill; Range;Squared Residuals\n")
    FittingParameterTxt.write(ST_exp(spatialDataList,spatialFittingRange))
    FittingParameterTxt.write("\n")
    FittingParameterTxt.write("Exponential Temporal fitting result:\n")
    FittingParameterTxt.write(ST_exp(temporalDataList,temporalFittingRange))
    FittingParameterTxt.close()

