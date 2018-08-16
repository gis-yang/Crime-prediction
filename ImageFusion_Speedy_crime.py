'''
Module Name:  ImageFusion.py

Created on Apr 22, 2014

@author: Bo Yang

All Rights Reserved

Version 1.0 released in Nov 12, 2017.

#modified 2018.02.01 added window size
#Same resolution fusion 
'''
import arcpy
import numpy as np
import math
import time
from numpy import matrix, linalg,asarray, squeeze
from arcpy.sa import *


#===================================================================================
#               Set Module Parameter
#===================================================================================

NODATA = -3.40282346639e+038
threshhold_nodata = 500


inputRaster= arcpy.gp.GetParameterAsText(0)
inputRasters= arcpy.gp.GetParameterAsText(1)
parameterTextFile= arcpy.gp.GetParameterAsText(2)
fitModelType= arcpy.gp.GetParameterAsText(3)
boolQuickMode= arcpy.gp.GetParameterAsText(4)
outputFusionRaster= arcpy.gp.GetParameterAsText(5)
boolUncertainty= arcpy.gp.GetParameterAsText(6)
outputuncertaintyRaster= arcpy.gp.GetParameterAsText(7)
secondFineInput= arcpy.gp.GetParameterAsText(8)
STFactor= float(arcpy.gp.GetParameterAsText(9))
ws= int(arcpy.gp.GetParameterAsText(10))
inputSeriesRasters= inputRasters.split(';')

arcpy.env.overwriteOutput = True

#ws = 7 #window size

#================load the parameter from text file=================================

txtfile=open(parameterTextFile,'r')
allLines=txtfile.readlines()
txtfile.close()
spatialParameterStr=str(allLines[1])
temporalParameterStr=str(allLines[3])

spatialParameter = [float(s) for s in spatialParameterStr.split(';') if s[0].isdigit()]
temporalParameter = [float(s) for s in temporalParameterStr.split(';') if s[0].isdigit()]

arcpy.AddMessage("Reading the " + fitModelType + " fitting model parameters")
arcpy.AddMessage(spatialParameter)
arcpy.AddMessage(temporalParameter)

#===================================================================================
#                   input fine Raster
#===================================================================================

fine1=Raster(inputRaster)
mFineRaster= arcpy.RasterToNumPyArray(fine1)*STFactor

#===================================================================================
#                     Two fine raster combine
#===================================================================================
if str(secondFineInput) == 'true':
    fine2=Raster(secondFineInput)
    mFineRaster2= arcpy.RasterToNumPyArray(fine2)
    mFineRaster = (mFineRaster + mFineRaster2)/2

#===================================================================================
#                   input Coarse Raster
#===================================================================================
mCoarseRasters=[]
for item in inputSeriesRasters[-3:]:
    mCoarseRasters.append(arcpy.RasterToNumPyArray(Raster(str(item))))


trendSum=0.0
trendCount=0
for j in range(mCoarseRasters[-1].shape[0]):
    for i in range(mCoarseRasters[-1].shape[1]):
        StackCoarse=[]
        if 0 < mCoarseRasters[-1][j][i] < threshhold_nodata:
            trendSum=trendSum+mCoarseRasters[-1][j][i]
            trendCount=trendCount+1
mMODIStrend = trendSum/trendCount

arcpy.AddMessage("The trend for the primary image is: " + str(mMODIStrend))



#===================================================================================
#                       Create New Raster
#===================================================================================
  
  
  
  
mNewRaster=np.asarray([[0.0 for i in range(mFineRaster.shape[1])] for j in range(mFineRaster.shape[0])])
mUncertyRaster=np.asarray([[0.0 for i in range(mFineRaster.shape[1])] for j in range(mFineRaster.shape[0])])
mNewCoarse=np.asarray([[0.0 for i in range(mFineRaster.shape[1])] for j in range(mFineRaster.shape[0])])


#================================================================================================================
#                                   Precalculate all the weights
#================================================================================================================
arcpy.AddMessage("Precalculate all the weights")
smallWindowpixels = []
if str(boolQuickMode) == 'false':
    windowBegin,windowEnd = 0,ws
else:
    windowBegin,windowEnd = 0,ws
for m in range(windowBegin,windowEnd):
    for n in range(windowBegin,windowEnd):
        smallWindowpixels.append(m*ws+n)


#================================================================================================================
#                   construct the matrix
#================================================================================================================

# Define Gaussian and Exponential function
def Gaussian(t,s,p,x):
    result = t+s*(1-math.exp(-(x**2)/(p**2)))
    return result

def Exponential(t,s,p,x):
    result = t+ s*(1-math.exp(-x/p))
    return result

# Define Covariance Function
if str(fitModelType) == 'Exponential':
    def spatialExpCovariance(x):
        return (spatialParameter[0] + spatialParameter[1]) - Exponential(spatialParameter[0], spatialParameter[1], spatialParameter[2], x)
    def temporalExpCovariance(x):
        return (temporalParameter[0] + temporalParameter[1]) - Exponential(temporalParameter[0], temporalParameter[1], temporalParameter[2], x)
elif str(fitModelType) == 'Gaussian':
    def spatialExpCovariance(x):
        return (spatialParameter[0] + spatialParameter[1]) - Gaussian(spatialParameter[0], spatialParameter[1], spatialParameter[2], x)
    def temporalExpCovariance(x):
        return (temporalParameter[0] + temporalParameter[1]) - Gaussian(temporalParameter[0], temporalParameter[1], temporalParameter[2], x)


#Spatial distance of fine resolution in UL sub-matrix
def Sdis(i,j):
    return 0.4 + math.sqrt((((j/ws)-(i/ws))**2)+(((j%ws)-(i%ws))**2))

#Spatial distance of coarse resolution in T1-T7 sub-matrix
#def Sdis2(i,j):
#    return math.sqrt((((j/3)-(i/3))**2)+(((j%3)-(i%3))**2))

#Spatial distance of cross resolution in UR sub-matrix
#def SdisCro(i,j):
#    return math.sqrt((((j/27)-((i/3)*9+4))**2)+((j%27)-((i%3)*9+4))**2)


#Final Matrix Function (water)
def finalMatrixWater(totalFinePixels=range(ws**2),totalCoaPixels=[range(ws**2),range(ws**2),range(ws**2)],CenterPx = ws**2/2):
    tp=totalFinePixels
    cp=totalCoaPixels
    Rk=bool(tp)+bool(cp[0] or cp[1] or cp[2])
    matrixSize=len(tp)+len(cp[0])+len(cp[1])+len(cp[2])+Rk
    # Create spatial(upper-left) sub-matrix(729*729)
    mUL=[[0 for i in tp] for j in tp]
    for i in tp:
        for j in tp:
            mUL[tp.index(j)][tp.index(i)]= spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(0)
    # Create Diagonal T1-T7 matrixs (9*9)
    mT11=[[0 for i in cp[0]] for j in cp[0]]
    mT22=[[0 for i in cp[1]] for j in cp[1]]
    mT33=[[0 for i in cp[2]] for j in cp[2]]
    for i in cp[0]:
        for j in cp[0]:
            mT11[cp[0].index(j)][cp[0].index(i)]= spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(0)
    for i in cp[1]:
        for j in cp[1]:
            mT22[cp[1].index(j)][cp[1].index(i)]= spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(0)
    for i in cp[2]:
        for j in cp[2]:
            mT33[cp[2].index(j)][cp[2].index(i)]= spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(0)
    # Create other T13-T31 matrixes (9*9)
    mT12=[[0 for i in cp[1]] for j in cp[0]]
    mT21=[[0 for i in cp[0]] for j in cp[1]] 
    mT13=[[0 for i in cp[2]] for j in cp[0]] 
    mT31=[[0 for i in cp[0]] for j in cp[2]] 
    mT23=[[0 for i in cp[2]] for j in cp[1]]
    mT32=[[0 for i in cp[1]] for j in cp[2]]
    for j in cp[1]:
        for i in cp[2]:
            mT23[cp[1].index(j)][cp[2].index(i)] = spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(1)
            mT32[cp[2].index(i)][cp[1].index(j)] = spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(1)
    for j in cp[0]:
        for i in cp[1]:
            mT12[cp[0].index(j)][cp[1].index(i)] = spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(1)
            mT21[cp[1].index(i)][cp[0].index(j)] = spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(1)
    for j in cp[0]:
        for i in cp[2]:      
            mT13[cp[0].index(j)][cp[2].index(i)] = spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(2)
            mT31[cp[2].index(i)][cp[0].index(j)] = spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(2)
    #Create Upper Right sub-matrix(9*729)
    mUR1=[[0 for i in cp[0]] for j in tp]
    mUR2=[[0 for i in cp[1]] for j in tp]
    mUR3=[[0 for i in cp[2]] for j in tp]
    mLL1=[[0 for i in tp] for j in cp[0]]
    mLL2=[[0 for i in tp] for j in cp[1]]
    mLL3=[[0 for i in tp] for j in cp[2]]
    for j in tp:
        for i in cp[0]:
            mUR1[tp.index(j)][cp[0].index(i)]=spatialExpCovariance(Sdis(i,j))*temporalExpCovariance(1)
            mLL1[cp[0].index(i)][tp.index(j)]=spatialExpCovariance(Sdis(i,j))*temporalExpCovariance(1)
    for j in tp:
        for i in cp[1]:
            mUR2[tp.index(j)][cp[1].index(i)]=spatialExpCovariance(Sdis(i,j))*temporalExpCovariance(2)
            mLL2[cp[1].index(i)][tp.index(j)]=spatialExpCovariance(Sdis(i,j))*temporalExpCovariance(2)
    for j in tp:
        for i in cp[2]:
            mUR3[tp.index(j)][cp[2].index(i)]=spatialExpCovariance(Sdis(i,j))*temporalExpCovariance(3)
            mLL3[cp[2].index(i)][tp.index(j)]=spatialExpCovariance(Sdis(i,j))*temporalExpCovariance(3)
    #Create Lower Left Quasi-diagonal matrix(747*4)
    mDiagLeft=[[0 for i in range(matrixSize-Rk)] for j in range(Rk)]
    RkIndex=0
    if len(tp)!=0:
        for i in range(len(tp)):
            mDiagLeft[RkIndex][i]=1
        RkIndex=RkIndex+1
    if len(cp[0])!=0:
        for i in range(len(tp), len(tp)+len(cp[0])):
            mDiagLeft[RkIndex][i]=1
    if len(cp[1])!=0:
        for i in range(len(tp)+len(cp[0]), len(tp)+len(cp[0])+len(cp[1])):
            mDiagLeft[RkIndex][i]=1
    if len(cp[2])!=0:  
        for i in range(len(tp)+len(cp[0])+len(cp[1]), len(tp)+len(cp[0])+len(cp[1])+len(cp[2])):
            mDiagLeft[RkIndex][i]=1
    # Create Upper Right Quasi-diagonal matrix(4*747)
    mDiagRight=[[0 for i in range(Rk)] for j in range(matrixSize-Rk)]
    for j in range(matrixSize-Rk):
        for i in range(Rk):
            mDiagRight[j][i]=mDiagLeft[i][j]
    # Create Lower Right All zero matrix(4*4)
    mZero=[[0 for i in range(Rk)] for j in range(Rk)] 
    #Create The Final Matrix(749*749)
    mFinalMatrix=[[0 for i in range(matrixSize)] for j in range(matrixSize)]
    for j in range(len(tp)):
        mFinalMatrix[j]=mUL[j]+mUR1[j] + mUR2[j] + mUR3[j] +mDiagRight[j]
    for j in range(len(cp[0])):
        mFinalMatrix[len(tp)+j]=mLL1[j] + mT11[j] + mT12[j] + mT13[j]+ mDiagRight[len(tp)+j]
    for j in range(len(cp[1])):
        mFinalMatrix[len(tp)+len(cp[0])+j]=mLL2[j] + mT21[j] + mT22[j] + mT23[j] + mDiagRight[len(tp)+len(cp[0])+j]
    for j in range(len(cp[2])):
        mFinalMatrix[len(tp)+len(cp[0])+len(cp[1])+j]=mLL3[j] + mT31[j] + mT32[j] +mT33[j] + mDiagRight[len(tp)+len(cp[0])+len(cp[1])+j]
    for j in range(Rk):
        mFinalMatrix[len(tp)+len(cp[0])+len(cp[1])+len(cp[2])+j]=mDiagLeft[j]+mZero[j]
    mmFinal=matrix(mFinalMatrix)
    #=================================================================================
    #               CY=B   ==>   Y=(C**-1)B   C is the FinalMatrix 
    #               B is the vector from predicted location every locations
    #               Below is generating B vector
    #=================================================================================
    matrixB = [0 for i in range(matrixSize)]
    matrixB[matrixSize-1]=0.9
    matrixB[matrixSize-2]=0.1
    for i in tp:
        matrixB[tp.index(i)]=spatialExpCovariance(Sdis(i,CenterPx))*temporalExpCovariance(3)
    for i in cp[0]:
        matrixB[len(tp)+cp[0].index(i)]=spatialExpCovariance(Sdis(i,CenterPx))*temporalExpCovariance(2)
    for i in cp[1]:
        matrixB[len(tp)+len(cp[0])+cp[1].index(i)]=spatialExpCovariance(Sdis(i,CenterPx))*temporalExpCovariance(1)
    for i in cp[2]:
        matrixB[len(tp)+len(cp[0])+len(cp[1])+cp[2].index(i)]=spatialExpCovariance(Sdis(i,CenterPx))*temporalExpCovariance(0)
    m1Result=squeeze(asarray(linalg.solve(mmFinal, matrix(matrixB).T)))
    if CenterPx == ws**2/2: # only save the center pixel matrix
        np.savetxt('mFinalMatrix.out', mFinalMatrix, delimiter=',')
        np.savetxt('matrixB.out', matrixB, delimiter=',')
        np.savetxt('m1Result.out', m1Result, delimiter=',')
    if str(boolUncertainty) == 'true':
        #start uncertainty calculation
        term1=spatialExpCovariance(0)*temporalExpCovariance(0)*2
        term2=0.0
        term3=0.0
        term4=0.0
        term5=0.0
        term6=0.0
        #------------------------------
        for j in range(len(tp)):
            for i in range(len(tp)):
                term2=m1Result[j]*m1Result[i]*mUL[j][i]+term2
        #------------------------------
        for j in range(len(cp[0])):
            for i in range(len(cp[0])):
                term3=m1Result[len(tp)+j]*m1Result[len(tp)+i]*mT11[j][i]+term3
        for j in range(len(cp[1])):
            for i in range(len(cp[1])):
                term3=m1Result[len(tp)+len(cp[0])+j]*m1Result[len(tp)+len(cp[0])+i]*mT22[j][i]+term3
        for j in range(len(cp[2])):
            for i in range(len(cp[2])):
                term3=m1Result[len(tp)+len(cp[0])+len(cp[1])+j]*m1Result[len(tp)+len(cp[0])+len(cp[1])+i]*mT33[j][i]+term3
        for j in range(len(cp[0])):
            for i in range(len(cp[1])):
                term3=m1Result[len(tp)+j]*m1Result[len(tp)+len(cp[0])+i]*mT12[j][i]+term3
        for j in range(len(cp[0])):
            for i in range(len(cp[2])):
                term3=m1Result[len(tp)+j]*m1Result[len(tp)+len(cp[0])+len(cp[1])+i]*mT13[j][i]+term3
        for j in range(len(cp[1])):
            for i in range(len(cp[2])):
                term3=m1Result[len(tp)+len(cp[0])+j]*m1Result[len(tp)+len(cp[0])+len(cp[1])+i]*mT23[j][i]+term3
        for j in range(len(cp[1])):
            for i in range(len(cp[0])):
                term3=m1Result[len(tp)+len(cp[0])+j]*m1Result[len(tp)+i]*mT21[j][i]+term3
        for j in range(len(cp[2])):
            for i in range(len(cp[0])):
                term3=m1Result[len(tp)+len(cp[0])+len(cp[1])+j]*m1Result[len(tp)+i]*mT31[j][i]+term3
        for j in range(len(cp[2])):
            for i in range(len(cp[1])):
                term3=m1Result[len(tp)+len(cp[0])+len(cp[1])+j]*m1Result[len(tp)+len(cp[0])+i]*mT32[j][i]+term3
        #------------------------------
        for j in range(len(tp)):
            for i in range(len(cp[0])):
                term4=2*m1Result[j]*m1Result[len(tp)+i]*mUR1[j][i] + term4
        for j in range(len(tp)):
            for i in range(len(cp[1])):
                term4=2*m1Result[j]*m1Result[len(tp)+len(cp[0])+i]*mUR2[j][i] + term4
        for j in range(len(tp)):
            for i in range(len(cp[2])):
                term4=2*m1Result[j]*m1Result[len(tp)+len(cp[0])+len(cp[1])+i]*mUR3[j][i] + term4
        for j in range(len(tp)):
            term5 = 2* m1Result[j] * matrixB[j] + term5
        for i in range(len(cp[0])):
            term6= 2*m1Result[len(tp)+i] * matrixB[len(tp)+i] + term6
        for i in range(len(cp[1])):
            term6= 2*m1Result[len(tp)+len(cp[0])+i] * matrixB[len(tp)+len(cp[0])+i] + term6
        for i in range(len(cp[2])):
            term6= 2*m1Result[len(tp)+len(cp[0])+len(cp[1])+i] * matrixB[len(tp)+len(cp[0])+len(cp[1])+i] + term6
        uncertainValue=term1+term2+term3+term4-term5-term6
        return [m1Result, uncertainValue]
    else: 
        return [m1Result,0]

fullPixelWeights= [finalMatrixWater(totalFinePixels=range(ws**2),totalCoaPixels=[range(ws**2),range(ws**2),range(ws**2)],CenterPx = i) for i in range(ws**2)]
        
# try:
    # fullPixelWeights= [finalMatrixWater(totalFinePixels=smallWindowpixels,PixelInCenter = i) for i in range(ws**2)]
# except:
    # arcpy.AddMessage("Error in pre-calculate matrix")
arcpy.AddMessage("Pre-calculate matrix finished")
arcpy.AddMessage(str(fullPixelWeights[ws**2/2]))
#================================================================================================================
#                                      New matrix delete method
#================================================================================================================
#Construct the left C
def leftMatrixC(totalFinePixels=range(ws**2),totalCoaPixels=[range(ws**2),range(ws**2),range(ws**2)]):
    tp=totalFinePixels
    cp=totalCoaPixels
    Rk=bool(tp)+bool(cp[0] or cp[1] or cp[2])
    matrixSize=len(tp)+len(cp[0])+len(cp[1])+len(cp[2])+Rk
    #print "The size of Water matrix is", len(tp),"+",len(cp[0]),'+',len(cp[1]),'+',len(cp[2]),'+ 2=',matrixSize
    # Create spatial(upper-left) sub-matrix(729*729)
    mUL=[[0 for i in tp] for j in tp]
    for i in tp:
        for j in tp:
            mUL[tp.index(j)][tp.index(i)]= spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(2)
    # Create Diagonal T1-T7 matrixs (9*9)
    mT11=[[0 for i in cp[0]] for j in cp[0]]
    mT22=[[0 for i in cp[1]] for j in cp[1]]
    mT33=[[0 for i in cp[2]] for j in cp[2]]
    for i in cp[0]:
        for j in cp[0]:
            mT11[cp[0].index(j)][cp[0].index(i)]= spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(0)
    for i in cp[1]:
        for j in cp[1]:
            mT22[cp[1].index(j)][cp[1].index(i)]= spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(0)
    for i in cp[2]:
        for j in cp[2]:
            mT33[cp[2].index(j)][cp[2].index(i)]= spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(0)
    # Create other T13-T31 matrixes (9*9)
    mT12=[[0 for i in cp[1]] for j in cp[0]]
    mT21=[[0 for i in cp[0]] for j in cp[1]] 
    mT13=[[0 for i in cp[2]] for j in cp[0]] 
    mT31=[[0 for i in cp[0]] for j in cp[2]] 
    mT23=[[0 for i in cp[2]] for j in cp[1]]
    mT32=[[0 for i in cp[1]] for j in cp[2]]
    for j in cp[1]:
        for i in cp[2]:
            mT23[cp[1].index(j)][cp[2].index(i)] = spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(1)
            mT32[cp[2].index(i)][cp[1].index(j)] = spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(1)
    for j in cp[0]:
        for i in cp[1]:
            mT12[cp[0].index(j)][cp[1].index(i)] = spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(1)
            mT21[cp[1].index(i)][cp[0].index(j)] = spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(1)
    for j in cp[0]:
        for i in cp[2]:      
            mT13[cp[0].index(j)][cp[2].index(i)] = spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(2)
            mT31[cp[2].index(i)][cp[0].index(j)] = spatialExpCovariance(Sdis(j,i))*temporalExpCovariance(2)
    #Create Upper Right sub-matrix(9*729)
    mUR1=[[0 for i in cp[0]] for j in tp]
    mUR2=[[0 for i in cp[1]] for j in tp]
    mUR3=[[0 for i in cp[2]] for j in tp]
    mLL1=[[0 for i in tp] for j in cp[0]]
    mLL2=[[0 for i in tp] for j in cp[1]]
    mLL3=[[0 for i in tp] for j in cp[2]]
    for j in tp:
        for i in cp[0]:
            mUR1[tp.index(j)][cp[0].index(i)]=spatialExpCovariance(Sdis(i,j))*temporalExpCovariance(0)
            mLL1[cp[0].index(i)][tp.index(j)]=spatialExpCovariance(Sdis(i,j))*temporalExpCovariance(0)
    for j in tp:
        for i in cp[1]:
            mUR2[tp.index(j)][cp[1].index(i)]=spatialExpCovariance(Sdis(i,j))*temporalExpCovariance(1)
            mLL2[cp[1].index(i)][tp.index(j)]=spatialExpCovariance(Sdis(i,j))*temporalExpCovariance(1)
    for j in tp:
        for i in cp[2]:
            mUR3[tp.index(j)][cp[2].index(i)]=spatialExpCovariance(Sdis(i,j))*temporalExpCovariance(2)
            mLL3[cp[2].index(i)][tp.index(j)]=spatialExpCovariance(Sdis(i,j))*temporalExpCovariance(2)
    #Create Lower Left Quasi-diagonal matrix(747*4)
    mDiagLeft=[[0 for i in range(matrixSize-Rk)] for j in range(Rk)] 
    RkIndex=0
    if len(tp)!=0:
        for i in range(len(tp)):
            mDiagLeft[RkIndex][i]=1
        RkIndex=RkIndex+1
    if len(cp[0])!=0:
        for i in range(len(tp), len(tp)+len(cp[0])):
            mDiagLeft[RkIndex][i]=1
    if len(cp[1])!=0:
        for i in range(len(tp)+len(cp[0]), len(tp)+len(cp[0])+len(cp[1])):
            mDiagLeft[RkIndex][i]=1
#    RkIndex=RkIndex+1
    if len(cp[2])!=0:
        for i in range(len(tp)+len(cp[0])+len(cp[1]), len(tp)+len(cp[0])+len(cp[1])+len(cp[2])):
            mDiagLeft[RkIndex][i]=1
    # Create Upper Right Quasi-diagonal matrix(4*747)
    mDiagRight=[[0 for i in range(Rk)] for j in range(matrixSize-Rk)] 
    for j in range(matrixSize-Rk):
        for i in range(Rk):
            mDiagRight[j][i]=mDiagLeft[i][j]
    # Create Lower Right All zero matrix(4*4)
    mZero=[[0 for i in range(Rk)] for j in range(Rk)] 
    #Create The Final Matrix(749*749)
    mFinalMatrix=[[0 for i in range(matrixSize)] for j in range(matrixSize)] 
    for j in range(len(tp)):
        mFinalMatrix[j]=mUL[j]+mUR1[j] + mUR2[j] + mUR3[j] +mDiagRight[j]
    for j in range(len(cp[0])):
        mFinalMatrix[len(tp)+j]=mLL1[j] + mT11[j] + mT12[j] + mT13[j]+ mDiagRight[len(tp)+j]
    for j in range(len(cp[1])):
        mFinalMatrix[len(tp)+len(cp[0])+j]=mLL2[j] + mT21[j] + mT22[j] + mT23[j] + mDiagRight[len(tp)+len(cp[0])+j]
    for j in range(len(cp[2])):
        mFinalMatrix[len(tp)+len(cp[0])+len(cp[1])+j]=mLL3[j] + mT31[j] + mT32[j] +mT33[j] + mDiagRight[len(tp)+len(cp[0])+len(cp[1])+j]
    for j in range(Rk):
        mFinalMatrix[len(tp)+len(cp[0])+len(cp[1])+len(cp[2])+j]=mDiagLeft[j]+mZero[j]
    return matrix(mFinalMatrix)

mleftMatrixC = leftMatrixC(totalFinePixels=smallWindowpixels)

def rightMatrixB(totalFinePixels=range(ws**2),totalCoaPixels=[range(ws**2),range(ws**2),range(ws**2)],PixelInCenter=ws**2/2):
    tp=totalFinePixels
    cp=totalCoaPixels
    Rk=bool(tp)+bool(cp[0] or cp[1] or cp[2])
    matrixSize=len(tp)+len(cp[0])+len(cp[1])+len(cp[2])+Rk
    matrixB = [0 for i in range(matrixSize)]
    matrixB[matrixSize-1]=1
    for i in tp:
        matrixB[tp.index(i)]=spatialExpCovariance(Sdis(i,PixelInCenter))*temporalExpCovariance(3)
    for i in cp[0]:
        matrixB[len(tp)+cp[0].index(i)]=spatialExpCovariance(Sdis(i,PixelInCenter))*temporalExpCovariance(2)
    for i in cp[1]:
        matrixB[len(tp)+len(cp[0])+cp[1].index(i)]=spatialExpCovariance(Sdis(i,PixelInCenter))*temporalExpCovariance(1)
    for i in cp[2]:
        matrixB[len(tp)+len(cp[0])+len(cp[1])+cp[2].index(i)]=spatialExpCovariance(Sdis(i,PixelInCenter))*temporalExpCovariance(0)
    return matrix(matrixB)

mrightMatrixB = [rightMatrixB(totalFinePixels=smallWindowpixels,totalCoaPixels=[range(ws**2),range(ws**2),range(ws**2)],PixelInCenter = i) for i in range(ws**2)]

def solveDelMatrix(fineDelete =np.array([]), coarseDelete = np.array([[],[],[]]),PixelInCenter=ws**2/2):
    deleteNum = np.concatenate((fineDelete,
                coarseDelete[0]+len(smallWindowpixels),
                coarseDelete[1] + (ws**2)+len(smallWindowpixels),
                coarseDelete[2] + (ws**2) + (ws**2)+len(smallWindowpixels)))
    if len(fineDelete)== len(smallWindowpixels):
        deleteNum = np.concatenate((fineDelete,
                                    coarseDelete[0]+len(smallWindowpixels),
                                    coarseDelete[1] + (ws**2)+len(smallWindowpixels),
                                    coarseDelete[2] + (ws**2) + (ws**2)+len(smallWindowpixels),
                                    [(ws**2) + (ws**2) + (ws**2)+len(smallWindowpixels)]))
    mleftMatrixDel = np.delete(np.delete(mleftMatrixC,deleteNum,0),deleteNum,1)
    mrightMatrixDel = np.delete(mrightMatrixB[PixelInCenter],deleteNum,1)
    return squeeze(asarray(linalg.solve(matrix(mleftMatrixDel), matrix(mrightMatrixDel).T)))


#==========================================================================================
#                             fusion Speed up
#==========================================================================================
def fuseRasterSpeedy(Raster=mFineRaster,coaRasters=mCoarseRasters):
    badPixelCount = 0
    for j in range(Raster.shape[0]):
#        if (j % (Raster.shape[0]/100) == 0):
#            arcpy.AddMessage(str(int(j*1.0/Raster.shape[0]*100)) + "%, Row Number: "+ str(j) +", "+ str(round(badPixelCount*100.0/(j*Raster.shape[1]),2)) + "% missing pixels now")
        for i in range(Raster.shape[1]):
            # if j-ws/2 <0: j = ws/2
            # if j+ws/2 > Raster.shape[0]-1: j = Raster.shape[0]-ws/2
            # if i-ws/2 <0: i = ws/2
            # if i+ws/2 > Raster.shape[0]-1: i = Raster.shape[0]-ws/2
            validPixelsValue = Raster[j-ws/2+windowBegin:j-ws/2+windowEnd,i-ws/2+windowBegin:i-ws/2+windowEnd].flatten()
            validCoaPixelsValue = np.array([coaRasters[k][j-ws/2:j-ws/2+ws,i-ws/2:i-ws/2+ws].flatten() for k in range(3)])
            if validPixelsValue[validPixelsValue<-threshhold_nodata].any() or validCoaPixelsValue[validCoaPixelsValue<-threshhold_nodata].any():
                validPixels = np.where(validPixelsValue>-threshhold_nodata)[0].tolist()
                valid_Pixels = np.where(validPixelsValue<-threshhold_nodata)[0]
                validCoaPixels0 = np.where(validCoaPixelsValue[0]>-threshhold_nodata)[0].tolist()
                validCoaPixels1 = np.where(validCoaPixelsValue[1]>-threshhold_nodata)[0].tolist()
                validCoaPixels2 = np.where(validCoaPixelsValue[2]>-threshhold_nodata)[0].tolist()
                if str(boolUncertainty) == 'true':
                    mUncertyRaster[j][i] = finalMatrixWater([item for item in smallWindowpixels if item not in valid_Pixels.tolist()],[validCoaPixels0,validCoaPixels1,validCoaPixels2],ws**2/2)[1]
                if len(validCoaPixels0)+len(validCoaPixels1)+len(validCoaPixels2)<3:
                    continue  
                valid_CoaPixels0 = np.where(validCoaPixelsValue[0]<-threshhold_nodata)[0]
                valid_CoaPixels1 = np.where(validCoaPixelsValue[1]<-threshhold_nodata)[0]
                valid_CoaPixels2 = np.where(validCoaPixelsValue[2]<-threshhold_nodata)[0]
                returnValue=solveDelMatrix(valid_Pixels,[valid_CoaPixels0,valid_CoaPixels1,valid_CoaPixels2],ws**2/2)#change here if want to apply changing-support to uncertainty
                pixelValue32=np.float32(np.concatenate((validPixelsValue[validPixelsValue>-threshhold_nodata], validCoaPixelsValue[0][validCoaPixelsValue[0]>-threshhold_nodata], validCoaPixelsValue[1][validCoaPixelsValue[1]>-threshhold_nodata], validCoaPixelsValue[2][validCoaPixelsValue[2]>-threshhold_nodata])))
                coefficient32=np.float32(returnValue[:len(validPixels)+len(validCoaPixels0)+len(validCoaPixels1)+len(validCoaPixels2)])
                badPixelCount+=1
            else:
                returnValue = fullPixelWeights[ws**2/2]
                pixelValue32=np.float32(np.concatenate((validPixelsValue, validCoaPixelsValue[0], validCoaPixelsValue[1], validCoaPixelsValue[2])))
                coefficient32=np.float32(returnValue[0][:len(validPixelsValue)+len(validCoaPixelsValue[0])+len(validCoaPixelsValue[1])+len(validCoaPixelsValue[2])])
                if str(boolUncertainty) == 'true':
                    mUncertyRaster[j][i] = fullPixelWeights[ws**2/2][1]
            mNewRaster[j][i]=VectorMultiple(pixelValue32,coefficient32).sum()

#===================================================================================
#                                 Cuda model
#=================================================================================== 
def VectorMultiple(pixelValue,coefficient):
    return pixelValue * coefficient

#==========================================================================================================
#                                   save raster
#==========================================================================================================


def saveRaster(Raster=mNewRaster, Name=outputFusionRaster, refer=Raster(inputRaster), cellSizeDiv=1):
#save raster
    fineBottom = float(arcpy.GetRasterProperties_management(refer, "BOTTOM").getOutput(0))
    fineLeft = float(arcpy.GetRasterProperties_management(refer, "LEFT").getOutput(0))
    fineXsize = float(arcpy.GetRasterProperties_management(refer, "CELLSIZEX").getOutput(0))/cellSizeDiv
    fineYsize = float(arcpy.GetRasterProperties_management(refer, "CELLSIZEY").getOutput(0))/cellSizeDiv
    arcpy.NumPyArrayToRaster(Raster,arcpy.Point(fineLeft,fineBottom),fineXsize, fineYsize).save(Name)

#==========================================================================================
#                Images detrend(trend is the first day raster)
#==========================================================================================  

def mFineFixedDetrend(Raster=mFineRaster, trend=mMODIStrend):
    for j in range(Raster.shape[0]):
        for i in range(Raster.shape[1]):
            if abs(Raster[j][i])<threshhold_nodata:
                Raster[j][i]=Raster[j][i] - trend
  
  
def mCoarseFixedDetrend(Raster=mCoarseRasters,trend=mMODIStrend):
    for j in range(mCoarseRasters[0].shape[0]):
        for i in range(mCoarseRasters[0].shape[1]):
            for k in range(len(mCoarseRasters)):
                if abs(Raster[k][j][i])<threshhold_nodata:
                    Raster[k][j][i]=Raster[k][j][i]-trend


#======================================================================================
#               run the fusion programe over the fine raster
#======================================================================================

def finalRun_detrend_Rsmp():
    arcpy.AddMessage(time.strftime("%m-%d %X",time.localtime())+" Detrending the coarse images")
    mCoarseFixedDetrend(mCoarseRasters)  
    #
    arcpy.AddMessage(time.strftime("%m-%d %X",time.localtime())+" Detrending the fine image")
    #mFineFixedDetrend(mFineRaster)
    #
    arcpy.AddMessage(time.strftime("%m-%d %X",time.localtime())+" Fusing the images")
    fuseRasterSpeedy()
    #
    arcpy.AddMessage(time.strftime("%m-%d %X",time.localtime())+" Adding the trend to the image")
    mFineFixedDetrend(mNewRaster, -(mMODIStrend))
    #
    arcpy.AddMessage(time.strftime("%m-%d %X",time.localtime())+" Saving the result")
    
    saveRaster()
    #
    arcpy.AddMessage(time.strftime("%m-%d %X",time.localtime())+outputFusionRaster+ " has finished")
#try:
finalRun_detrend_Rsmp()
#except:
    #arcpy.AddMessage("Unexpected error")

if str(boolUncertainty) == 'true':
    saveRaster(mUncertyRaster, outputuncertaintyRaster)



