'''
Module Name: SemiVariog

Created on Apr 18, 2014

@author: Bo Yang

All Rights Reserved
'''
import arcpy
import numpy
import math
import random
from arcpy.sa import *
import time


#===================================================================================
#               Set Module Parameter
#===================================================================================


NODATA= -3.40282346639e+038
inputRaster= arcpy.gp.GetParameterAsText(0)
sampleRatio= float(arcpy.gp.GetParameterAsText(1))
boolDetrend= arcpy.gp.GetParameterAsText(2)
inputRasters= arcpy.gp.GetParameterAsText(3)
spatialDistance= int(arcpy.gp.GetParameterAsText(4))
outputTxtSpatialVariogram= arcpy.gp.GetParameterAsText(5)
outputTxtTemporalVariogram= arcpy.gp.GetParameterAsText(6)
temporalOnly= arcpy.gp.GetParameterAsText(7)

inputSeriesRasters= inputRasters.split(';')

mFineRaster=arcpy.RasterToNumPyArray(Raster(inputRaster))

mCoarseRasters=[]
for item in inputSeriesRasters:
    mCoarseRasters.append(arcpy.RasterToNumPyArray(Raster(str(item))))

#===================================================================================
#               Detrended module
#===================================================================================

def imagesDetrend():
    arcpy.AddMessage(str(time.strftime("%m-%d %X",time.localtime())) + "; Start detrend")
#---------------------------Calculate the trend-----------------------------------------------------
    mMODIStrend = numpy.asarray([[0.0 for i in range(mCoarseRasters[0].shape[1])] for j in range(mCoarseRasters[0].shape[0])])   # this is the trend raster
    for j in range(mMODIStrend.shape[0]):
        for i in range(mMODIStrend.shape[1]):
            StackCoarse=[]
            for k in range(len(mCoarseRasters)):
                if abs(mCoarseRasters[k][j][i])<10000:
                    StackCoarse.append(mCoarseRasters[k][j][i])
            if StackCoarse!=[]:
                mMODIStrend[j][i]=sum(StackCoarse)/len(StackCoarse)
#---------------------------Coase Images detrend----------------------------------------------------
    for j in range(mCoarseRasters[0].shape[0]):
        for i in range(mCoarseRasters[0].shape[1]):
            for k in range(len(mCoarseRasters)):
                if abs(mCoarseRasters[k][j][i])<10000:
                    mCoarseRasters[k][j][i]=mCoarseRasters[k][j][i]-mMODIStrend[j][i]
#------------------------Fine resolution Data detrend---------------------------------------------
    for j in range(mFineRaster.shape[0]):
        for i in range(mFineRaster.shape[1]):
            for n in range(mMODIStrend.shape[0]):
                for m in range(mMODIStrend.shape[1]):
                    if (j/(mFineRaster.shape[0]/mMODIStrend.shape[0])==n and i/(mFineRaster.shape[0]/mMODIStrend.shape[0])==m and abs(mFineRaster[j][i])<500000):
                        mFineRaster[j][i]=mFineRaster[j][i] - mMODIStrend[n][m]
    arcpy.AddMessage(str(time.strftime("%m-%d %X",time.localtime())) + "; Finish detrend")


#================================================================================
#                    Spatial Semi-variogram
#================================================================================
def SpatialSemiV(Raster=mFineRaster, DISTANCE=spatialDistance, SpatialSampleRatio=sampleRatio):
    arcpy.AddMessage(str(time.strftime("%m-%d %X",time.localtime())) + "; Start SpatialSemiV Calculation")
    spatialSV=[[0,0] for i in range (DISTANCE)]
    spatialVariogram = [0.0 for i in range (DISTANCE)]
    allPixelType=[]
    for j in range(Raster.shape[0]):
        for i in range(Raster.shape[1]):
            if (abs(Raster[j][i])<10000):
                allPixelType.append([j,i,Raster[j][i]])
    distanceBetweenPixel=0.0
    vriogramValue=0.0
    arcpy.AddMessage("Selected " + str(int(SpatialSampleRatio*len(allPixelType))) + " from total " + str(len(allPixelType)) + " points for Spatial calculation")
    allPixelType=random.sample(allPixelType,int(SpatialSampleRatio*len(allPixelType)))
    for j in range(len(allPixelType)):
        for i in range(j):         
            distanceBetweenPixel=math.sqrt(((allPixelType[j][0])-(allPixelType[i][0]))**2 + ((allPixelType[j][1])-(allPixelType[i][1]))**2)
            vriogramValue=(allPixelType[j][2]-allPixelType[i][2])**2
            if distanceBetweenPixel+0.5<DISTANCE:
                spatialSV[int(distanceBetweenPixel+0.5)][0]=spatialSV[int(distanceBetweenPixel+0.5)][0]+1
                spatialSV[int(distanceBetweenPixel+0.5)][1]=spatialSV[int(distanceBetweenPixel+0.5)][1]+vriogramValue
    for k in range(DISTANCE):
        if spatialSV[k][0]!=0:
            spatialVariogram[k]=0.5*(spatialSV[k][1]/spatialSV[k][0])
    arcpy.AddMessage(str(time.strftime("%m-%d %X",time.localtime())) + "SpatialSemiV Calculation Complete")
    TxtSpatialVariogram = open(outputTxtSpatialVariogram, 'w')
    TxtSpatialVariogram.write(str(spatialVariogram).strip('[]'))
    TxtSpatialVariogram.close()

#================================================================================
#                   Coarse resolution Temporal Semi-variogram
#================================================================================
def timeSemiV(h=len(mCoarseRasters)):
    arcpy.AddMessage(str(time.strftime("%m-%d %X",time.localtime())) + "; Start Temporal SemiV Calculation")
    timeSV=[0 for i in range(h)]
    EntiretotalCount = 0
    for distance in range(h):
        totalCount = 0
        myResult = 0
        for j in range(mCoarseRasters[0].shape[0]):
            for i in range(mCoarseRasters[0].shape[1]):
                for k in range(len(mCoarseRasters) - distance):
                    if (abs(mCoarseRasters[k][j][i]) < 10000 and abs(mCoarseRasters[k+distance][j][i]) < 10000):
                        myResult = myResult + (mCoarseRasters[k][j][i] - mCoarseRasters[k+distance][j][i])**2
                        totalCount = totalCount + 1
        EntiretotalCount =+ totalCount
        if totalCount!=0:
            timeSV[distance] = 0.5 * myResult/totalCount
    arcpy.AddMessage("Selected " + str(EntiretotalCount) + " points for Temporal calculation")
    arcpy.AddMessage(str(time.strftime("%m-%d %X",time.localtime())) + "Temporal SemiV Calculation Complete")
    TxtTemporalVariogram = open(outputTxtTemporalVariogram, 'w')
    TxtTemporalVariogram.write(str(timeSV).strip('[]'))
    TxtTemporalVariogram.close()


#================================================================================
#                   Coarse resolution Temporal Semi-variogram
#================================================================================
def SpatialCoaSemi(Raster = mCoarseRasters, CoaDISTANCE=spatialDistance,h=len(mCoarseRasters)):
    arcpy.AddMessage(str(time.strftime("%m-%d %X",time.localtime())) + "; Start SpatialSemiV Calculation")
    spatialCoaSV=[[0,0] for i in range(CoaDISTANCE)]
    spatialCoaVariogram = [0.0 for i in range(CoaDISTANCE)]
    allPixelType=[]
    for k in range(h):
        for j in range(Raster[k].shape[0]):
            for i in range(Raster[k].shape[1]):
                if (abs(Raster[k][j][i])<10000):
                    allPixelType.append([j,i,Raster[k][j][i]])
    distanceBetweenPixel=0.0
    vriogramValue=0.0
    for j in range(len(allPixelType)):
        for i in range(j):         
            distanceBetweenPixel=math.sqrt(((allPixelType[j][0])-(allPixelType[i][0]))**2 + ((allPixelType[j][1])-(allPixelType[i][1]))**2)
            vriogramValue=(allPixelType[j][2]-allPixelType[i][2])**2
            if distanceBetweenPixel+0.5<CoaDISTANCE:
                spatialCoaSV[int(distanceBetweenPixel+0.5)][0]=spatialCoaSV[int(distanceBetweenPixel+0.5)][0]+1
                spatialCoaSV[int(distanceBetweenPixel+0.5)][1]=spatialCoaSV[int(distanceBetweenPixel+0.5)][1]+vriogramValue
    for k in range(CoaDISTANCE):
        if spatialCoaSV[k][0]!=0:
            spatialCoaVariogram[k]=0.5*(spatialCoaSV[k][1]/spatialCoaSV[k][0])
    arcpy.AddMessage(str(time.strftime("%m-%d %X",time.localtime())) + "Temporal SemiV Calculation Complete")
    TxtTemporalVariogram = open(outputTxtTemporalVariogram, 'w')
    TxtTemporalVariogram.write(str(spatialCoaSV).strip('[]'))
    TxtTemporalVariogram.close()
#=============================================================================================
#                Run the Main Program
#=============================================================================================
if boolDetrend == "true":
    imagesDetrend()
if temporalOnly =="false":
    SpatialSemiV()
timeSemiV()
#SpatialCoaSemi()
