# -*- coding: utf-8 -*-
"""
Created on Wed May 29 10:52:06 2019

@author: sst
"""

import numpy   as np
import netCDF4 as nc  
import numpy.ma as ma
import struct
import time, os, sys                     # call current time for timestamp
from writenetcdf import writenetcdf # from ufz
from readnetcdf  import readnetcdf  # from ufz
from fread       import fread       # from ufz
from autostring  import astr        # from ufz
from pyproj      import Proj
import scipy.misc
import scipy.ndimage
from netCDF4 import Dataset
import glob
import re
import gdal
from MOD11_HDF_READ import getSceneMetadata
from MOD11_HDF_READ import Save_array_tiff
from MOD11_HDF_READ import getInputStructure
from skimage import io



#LAIjan = np.full([3600,3600], 0)
#LAIfeb = np.full([3600,3600], 0)
#LAImar = np.full([3600,3600], 0)
#LAIapr = np.full([3600,3600], 0)
#LAImay = np.full([3600,3600], 0)
#LAIjun = np.full([3600,3600], 0)
#LAIjul = np.full([3600,3600], 0)
#LAIaug = np.full([3600,3600], 0)
#LAIsep = np.full([3600,3600], 0)
#LAIoct = np.full([3600,3600], 0)
#LAInov = np.full([3600,3600], 0)
#LAIdec = np.full([3600,3600], 0)


#LAIjan = makeMapMonthCorrect(1,2003,2012)
#LAIfeb = makeMapMonthCorrect(2,2003,2012)
#LAImar = makeMapMonthCorrect(3,2003,2012)
#LAIapr = makeMapMonthCorrect(4,2003,2012)
#LAImay = makeMapMonthCorrect(5,2003,2012)
#LAIjun = makeMapMonthCorrect(6,2003,2012)
#LAIjul = makeMapMonthCorrect(7,2003,2012)
#LAIaug = makeMapMonthCorrect(8,2003,2012)
#LAIsep = makeMapMonthCorrect(9,2003,2012)
#LAIoct = makeMapMonthCorrect(10,2003,2012)
#LAInov = makeMapMonthCorrect(11,2003,2012)
#LAIdec = makeMapMonthCorrect(12,2003,2012)
#
#
#allLAIone = np.nansum(np.dstack((LAIjan, LAIfeb)),2)
#allLAItwo = np.nansum(np.dstack((LAImar, LAIapr)),2)
#allLAIthree = np.nansum(np.dstack((LAImay, LAIjun)),2)
#allLAIfour = np.nansum(np.dstack((LAIjul, LAIaug)),2)
#allLAIfive = np.nansum(np.dstack((LAIsep, LAIoct)),2)
#allLAIsix = np.nansum(np.dstack((LAInov, LAIdec)),2)
#
#result = np.divide(allLAI, 12)
##allLAI = np.nansum(np.dstack((LAIjan, LAIfeb, LAImar, LAIapr, LAImay, LAIjun, LAIjul, LAIaug, LAIsep, LAIoct, LAInov, LAIdec)),2)
#result = np.divide(allLAI, 12)
#
#   try:                
#       referenceMetadata=getSceneMetadata('Y:\\Gorka\\EUROPA\\MCD15A2.005\\MCD15A2.0052015129.mosaic.tif')
#       inputFieldNames=getInputStructure()
#       proj=referenceMetadata.get('proj')
#       geotransform=referenceMetadata.get('gt') 
#       #outputfile= 'F:\\Temporal\\Fg_Image.mosaiccacaca.tif' 
#       outputfile= 'E:\\TSEB_Europe\\Output_Data\\averageMaps\\LAI_TOTAL.mosaic.tif'
#       
#       Save_array_tiff(result, outputfile, proj,geotransform)            #Save the FG for the stations
#       print("saved")
#   except:
#       print("something went wrong")
#   print("success")
#   
#   
#realMonth = ""
    
def runForAll():
    makeMapMonthCorrect(1,2003,2015)
    makeMapMonthCorrect(2,2003,2015)
    makeMapMonthCorrect(3,2003,2015)
    makeMapMonthCorrect(4,2003,2015)
    makeMapMonthCorrect(5,2003,2015)
    makeMapMonthCorrect(6,2003,2015)
    makeMapMonthCorrect(7,2003,2015)
    makeMapMonthCorrect(8,2003,2015)
    makeMapMonthCorrect(9,2003,2015)
    makeMapMonthCorrect(10,2003,2015)
    makeMapMonthCorrect(11,2003,2015)
    makeMapMonthCorrect(12,2003,2015)
    
def makeMapMonthCorrect(month, yearFrom, yearTo):
    
    #month = month number in format xx
   #yearFrom/yearTo = year in format xxxx doesn't include yearTo.
   import datetime 
   tm = datetime.datetime(2014, 1, 1) + datetime.timedelta(165 - 1)    
   correctDate = tm.strftime('%Y_%m_%d')
    
    
   #yearFrom = 2003
   #yearTo = 2014
   #month = 7
   count = 0
   
   
   #zLength = yearTo - yearFrom + 1
       
   
   totalArrayH = np.full([3600,3600, 31], np.nan)
   totalArrayLE = np.full([3600,3600, 31], np.nan)
   
   #totalHYears = np.full([3600,3600, zLength], np.nan)
   #totalLEYears = np.full([3600,3600, zLength], np.nan)
   
   totalCountH = np.full([3600,3600], 0)
   totalCountLE = np.full([3600,3600],0)
   totalCountLST = np.full([3600,3600],0)
   totalCountLAI = np.full([3600,3600],0)
   totalLE = np.full([3600,3600], 0)
   totalLEstd = np.full([3600,3600], 0)
   totalH = np.full([3600,3600], 0)
   totalHstd = np.full([3600,3600], 0)
   totalLST = np.full([3600,3600], 0)
   totalLAI = np.full([3600,3600], 0)
   
   
   
   
   for year in range(yearFrom, yearTo):
       #print month
       
       totalArrayH = np.full([3600,3600, 31], np.nan)
       totalArrayLE = np.full([3600,3600, 31], np.nan)
       count = 0
       for day in range(1, 32):
           #print day
           
           
           #month = 5
           #day = 5
           #year = 2004
           if(month<10):
               month = '0'+str(month)
           month = str(month) 
           
           realMonth = month
           if(day<10):
               day = '0'+str(day)
           day = str(day)    
           
           print ('X:\\TSEB_Europe\\Output_Data\\netcdf_results\\InputsTSEB_'+str(year)+'_'+month+'_'+day, 'H_C')           
           try:
               
               H_C = readnetcdf('X:\\TSEB_Europe\\Output_Data\\netcdf_results\\InputsTSEB_'+str(year)+'_'+month+'_'+day+'.nc', 'H_C')
               H_S = readnetcdf('X:\\TSEB_Europe\\Output_Data\\netcdf_results\\InputsTSEB_'+str(year)+'_'+month+'_'+day+'.nc', 'H_S')
               LE_C = readnetcdf('X:\\TSEB_Europe\\Output_Data\\netcdf_results\\InputsTSEB_'+str(year)+'_'+month+'_'+day+'.nc', 'LE_C')  
               LE_S = readnetcdf('X:\\TSEB_Europe\\Output_Data\\netcdf_results\\InputsTSEB_'+str(year)+'_'+month+'_'+day+'.nc', 'LE_S')
               #LAI  = readnetcdf('X:\\TSEB_Europe\\Output_Data\\netcdf_results\\InputsTSEB_'+str(year)+'_'+month+'_'+day+'.nc', 'LAI')
               #LST  = readnetcdf('X:\\TSEB_Europe\\Output_Data\\netcdf_results\\InputsTSEB_'+str(year)+'_'+month+'_'+day+'.nc', 'LST')
               
               
               LEavg = io.imread('E:\TSEB_Europe\\Output_Data\\averageMaps\\AverageLE_Month_' +month+'.mosaic.tif')                
               
               Havg = io.imread('E:\TSEB_Europe\\Output_Data\\averageMaps\\AverageH_Month_' +month+'.mosaic.tif')
               print("found files....") 
               
               
               
               #----------------------------------------------------------------------
               H = np.add(H_S, H_C)
               H[H>800]=np.nan
               H[H<0]=np.nan
               
               totalCountH[np.where(H>=0)]+= 1
               
               
               totalH = np.nansum(np.dstack((totalH, H)),2)  
               
               
               LE = np.add(LE_S, LE_C)
               LE[LE>800]=np.nan
               LE[LE<0]=np.nan
               
               
               HminusMean = np.subtract(H, Havg)
               Hsquared = np.square(HminusMean)
               
               LEminusMean = np.subtract(LE, LEavg)
               LEsquared = np.square(LEminusMean)
               
               totalCountLE[np.where(LE>=0)]+= 1
               
               totalLE = np.nansum(np.dstack((totalLE, LE)),2)
               
               totalLEstd = np.nansum(np.dstack((totalLEstd, LEsquared)),2) #adding up all the values minus the mean and squared
               totalHstd = np.nansum(np.dstack((totalHstd, Hsquared)),2)
               
               #-----------------------------------------------------------------------------                
               
               #LST[LST>400]=np.nan
               #LST[np.isinf(LST)]=np.nan
               #LST[np.isneginf(LST)]=np.nan
               #LST[LST<100]=np.nan
               #LAI[LAI>8]=np.nan
               #LAI[LAI<0]=np.nan
               
               #totalCountLAI[np.where(LAI>=0)]+= 1
               #totalCountLST[np.where(LST>=0)]+= 1
               
               #totalLST = np.nansum(np.dstack((totalLST, LST)),2)
               #totalLAI = np.nansum(np.dstack((totalLAI, LAI)),2)
               
               #--------------------------------------------------------------------------------
#               t = np.zeros([3,3])
#               y = np.zeros([3,3])
#               
#               t[2,2] = 5
#               t[1,1] = 4
#               y[2,2] = np.nan
#               y[1,1] = 2
#               
#               g = np.subtract(t, y)
#               l = np.square(g)
#               u = np.sqrt(l)
#               
#               month = "07"
#               im = io.imread('E:\TSEB_Europe\\Output_Data\\averageMaps\\AverageLE_Month_' +month+'.mosaic.tif')
                #fgArray = im[:, :, 20]
           except:
               print('exceeded days probably')
           count= count + 1    
       
       print("averaging for the month")
       
   print("averaging across the years")    
   averageH = np.divide(totalH, totalCountH)
   averageLE = np.divide(totalLE, totalCountLE)
   
   averageLEstd = np.divide(totalLEstd, totalCountLE) #find averag of the std
   averageLEstdSquare = np.sqrt(averageLEstd) #square the stds
   
   averageHstd = np.divide(totalHstd, totalCountH) #find averag of the std
   averageHstdSquare = np.sqrt(averageHstd) #square the stds
   
   #averageLST = np.divide(totalLST, totalCountLST)
   #averageLAI = np.divide(totalLAI, totalCountLAI)
   
   
   
     
#   print realMonth
#   print ("saving total monthly")
#   if(realMonth=="01"):
#       LAIjan = averageLAI
#       return LAIjan
#   if(realMonth=="02"):
#       LAIfeb = averageLAI
#       return LAIfeb
#   if(realMonth=="03"):
#       LAImar = averageLAI
#       return LAImar
#   if(realMonth=="04"):
#       LAIapr = averageLAI 
#       return LAIapr
#   if(realMonth=="05"):
#       LAImay = averageLAI
#       return LAImay
#   if(realMonth=="06"):
#       LAIjun = averageLAI
#       return LAIjun
#   if(realMonth=="07"):
#       LAIjul = averageLAI
#       return LAIjul
#   if(realMonth=="08"):
#       LAIaug = averageLAI
#       return LAIaug
#   if(realMonth=="9"):
#       LAIsep = averageLAI
#       return LAIsep
#   if(realMonth=="10"):
#       LAIoct = averageLAI
#       return LAIoct
#   if(realMonth=="1"):
#       LAInov = averageLAI
#       return LAInov
#   if(realMonth=="1"):
#       LAIdec = averageLAI
#       return LAIdec
#   print ("saved total monthly")   
   #averageET = np.add(averageH, averageLE)
   

#   try:                
#       referenceMetadata=getSceneMetadata('Y:\\Gorka\\EUROPA\\MCD15A2.005\\MCD15A2.0052015129.mosaic.tif')
#       inputFieldNames=getInputStructure()
#       proj=referenceMetadata.get('proj')
#       geotransform=referenceMetadata.get('gt') 
#       #outputfile= 'F:\\Temporal\\Fg_Image.mosaiccacaca.tif' 
#       outputfile= 'E:\\TSEB_Europe\\Output_Data\\averageMaps\\LAI_Month_' + month +'.mosaic.tif'
#       
#       Save_array_tiff(averageLAI, outputfile, proj,geotransform)            #Save the FG for the stations
#       print("saved")
#   except:
#       print("something went wrong")
#   print("success")
#   
#   try:                
#       referenceMetadata=getSceneMetadata('Y:\\Gorka\\EUROPA\\MCD15A2.005\\MCD15A2.0052015129.mosaic.tif')
#       inputFieldNames=getInputStructure()
#       proj=referenceMetadata.get('proj')
#       geotransform=referenceMetadata.get('gt') 
#       #outputfile= 'F:\\Temporal\\Fg_Image.mosaiccacaca.tif' 
#       outputfile= 'E:\\TSEB_Europe\\Output_Data\\averageMaps\\LST_Month_' + month +'.mosaic.tif'
#       
#       Save_array_tiff(averageLST, outputfile, proj,geotransform)            #Save the FG for the stations
#       print("saved")
#   except:
#       print("something went wrong")
#   print("success")
#   
   try:                
       referenceMetadata=getSceneMetadata('Y:\\Gorka\\EUROPA\\MCD15A2.005\\MCD15A2.0052015129.mosaic.tif')
       inputFieldNames=getInputStructure()
       proj=referenceMetadata.get('proj')
       geotransform=referenceMetadata.get('gt') 
       #outputfile= 'F:\\Temporal\\Fg_Image.mosaiccacaca.tif' 
       outputfile= 'E:\\TSEB_Europe\\Output_Data\\averageMaps\\LE_std_' + month +'.mosaic.tif'
       
       Save_array_tiff(averageLEstdSquare, outputfile, proj,geotransform)            #Save the FG for the stations
       print("saved")
   except:
       print("something went wrong")
   print("success")
   
   try:                
       referenceMetadata=getSceneMetadata('Y:\\Gorka\\EUROPA\\MCD15A2.005\\MCD15A2.0052015129.mosaic.tif')
       inputFieldNames=getInputStructure()
       proj=referenceMetadata.get('proj')
       geotransform=referenceMetadata.get('gt') 
       #outputfile= 'F:\\Temporal\\Fg_Image.mosaiccacaca.tif' 
       outputfile= 'E:\\TSEB_Europe\\Output_Data\\averageMaps\\H_std_' + month +'.mosaic.tif'
       
       Save_array_tiff(averageHstdSquare, outputfile, proj,geotransform)            #Save the FG for the stations
       print("saved")
   except:
       print("something went wrong")
   print("success")
               
#   try:                
#       referenceMetadata=getSceneMetadata('Y:\\Gorka\\EUROPA\\MCD15A2.005\\MCD15A2.0052015129.mosaic.tif')
#       inputFieldNames=getInputStructure()
#       proj=referenceMetadata.get('proj')
#       geotransform=referenceMetadata.get('gt') 
#       #outputfile= 'F:\\Temporal\\Fg_Image.mosaiccacaca.tif' 
#       outputfile= 'E:\\TSEB_Europe\\Output_Data\\averageMaps\\CountLE_Month_' + month +'.mosaic.tif'
#       
#       Save_array_tiff(totalCountLE, outputfile, proj,geotransform)            #Save the FG for the stations
#       print("saved")
#   except:
#       print("something went wrong")
#   print("success")
#   
#   try:                
#       referenceMetadata=getSceneMetadata('Y:\\Gorka\\EUROPA\\MCD15A2.005\\MCD15A2.0052015129.mosaic.tif')
#       inputFieldNames=getInputStructure()
#       proj=referenceMetadata.get('proj')
#       geotransform=referenceMetadata.get('gt') 
#       #outputfile= 'F:\\Temporal\\Fg_Image.mosaiccacaca.tif' 
#       outputfile= 'E:\\TSEB_Europe\\Output_Data\\averageMaps\\AverageH_Month_' + month +'.mosaic.tif'
#       
#       Save_array_tiff(averageH, outputfile, proj,geotransform)            #Save the FG for the stations
#       print("saved")
#   except:
#       print("something went wrong")
#   print("success")
#
#   try:                
#       referenceMetadata=getSceneMetadata('Y:\\Gorka\\EUROPA\\MCD15A2.005\\MCD15A2.0052015129.mosaic.tif')
#       inputFieldNames=getInputStructure()
#       proj=referenceMetadata.get('proj')
#       geotransform=referenceMetadata.get('gt') 
#       #outputfile= 'F:\\Temporal\\Fg_Image.mosaiccacaca.tif' 
#       outputfile= 'E:\\TSEB_Europe\\Output_Data\\averageMaps\\AverageLE_Month_' + month +'.mosaic.tif'
#       
#       Save_array_tiff(averageLE, outputfile, proj,geotransform)            #Save the FG for the stations   
#       print("saved")
#   except:
#       print("something went wrong")
#   print("success") 
#
#   try:                
#       referenceMetadata=getSceneMetadata('Y:\\Gorka\\EUROPA\\MCD15A2.005\\MCD15A2.0052015129.mosaic.tif')
#       inputFieldNames=getInputStructure()
#       proj=referenceMetadata.get('proj')
#       geotransform=referenceMetadata.get('gt') 
#       #outputfile= 'F:\\Temporal\\Fg_Image.mosaiccacaca.tif' 
#       outputfile= 'E:\\TSEB_Europe\\Output_Data\\averageMaps\\H+LE_Month_' + month +'.mosaic.tif'
#       
#       Save_array_tiff(averageET, outputfile, proj,geotransform)            #Save the FG for the stations   
#       print("saved")
#   except:
#       print("something went wrong")
#   print("success")                 
   
def makeMapMonth(month, yearFrom, yearTo):
    
   #month = month number in format xx
   #yearFrom/yearTo = year in format xxxx doesn't include yearTo.
   import datetime 
   tm = datetime.datetime(2014, 1, 1) + datetime.timedelta(165 - 1)    
   correctDate = tm.strftime('%Y_%m_%d')
    
    
   yearFrom = 2003
   yearTo = 2014
   month = 3
   count = 0
   
   
   zLength = yearTo - yearFrom + 1
       
   
   totalArrayH = np.full([3600,3600, 31], np.nan)
   totalArrayLE = np.full([3600,3600, 31], np.nan)
   
   totalHYears = np.full([3600,3600, zLength], np.nan)
   totalLEYears = np.full([3600,3600, zLength], np.nan)
   
   for year in range(yearFrom, yearTo):
       #print month
       
       totalArrayH = np.full([3600,3600, 31], np.nan)
       totalArrayLE = np.full([3600,3600, 31], np.nan)
       count = 0
       for day in range(1, 32):
           #print day
           
           
           #month = 5
           #day = 5
           #year = 2004
           if(month<10):
               month = '0'+str(month)
           month = str(month) 
           
           if(day<10):
               day = '0'+str(day)
           day = str(day)    
           
           print ('X:\\TSEB_Europe\\Output_Data\\netcdf_results\\InputsTSEB_'+str(year)+'_'+month+'_'+day, 'H_C')           
           try:
               
               H_C = readnetcdf('X:\\TSEB_Europe\\Output_Data\\netcdf_results\\InputsTSEB_'+str(year)+'_'+month+'_'+day+'.nc', 'H_C')
               H_S = readnetcdf('X:\\TSEB_Europe\\Output_Data\\netcdf_results\\InputsTSEB_'+str(year)+'_'+month+'_'+day+'.nc', 'H_S')
               LE_C = readnetcdf('X:\\TSEB_Europe\\Output_Data\\netcdf_results\\InputsTSEB_'+str(year)+'_'+month+'_'+day+'.nc', 'LE_C')  
               LE_S = readnetcdf('X:\\TSEB_Europe\\Output_Data\\netcdf_results\\InputsTSEB_'+str(year)+'_'+month+'_'+day+'.nc', 'LE_S')
               
               totalH = np.add(H_S, H_C)
               totalArrayH[:,:,count] = totalH
               
               print(count)
               totalLE = np.add(LE_S, LE_C)
               totalArrayLE[:,:,count] = totalLE
           except:
               print('exceeded days probably')
           count= count + 1    
       
       print("averaging for the month")
       avgH = np.nanmean(totalArrayH, axis=2)
       avgLE = np.nanmean(totalArrayLE, axis=2)
       
       totalHYears[:,:,yearTo-yearFrom] = avgH
       totalLEYears[:,:,yearTo-yearFrom] = avgLE
   print("averaging across the years")    
   actualAverageH = np.nanmean(totalHYears, axis=2)
   actualAverageLE = np.nanmean(totalLEYears, axis=2)

     
   totalET = np.add(actualAverageH, actualAverageLE)             
   
   try:                
       referenceMetadata=getSceneMetadata('Y:\\Gorka\\EUROPA\\MCD15A2.005\\MCD15A2.0052015129.mosaic.tif')
       inputFieldNames=getInputStructure()
       proj=referenceMetadata.get('proj')
       geotransform=referenceMetadata.get('gt') 
       #outputfile= 'F:\\Temporal\\Fg_Image.mosaiccacaca.tif' 
       outputfile= 'E:\\TSEB_Europe\\Output_Data\\averageMaps\\AverageH_Month_' + month +'.mosaic.tif'
       
       Save_array_tiff(actualAverageH, outputfile, proj,geotransform)            #Save the FG for the stations
       print("saved")
   except:
       print("something went wrong")
   print("success")

   try:                
       referenceMetadata=getSceneMetadata('Y:\\Gorka\\EUROPA\\MCD15A2.005\\MCD15A2.0052015129.mosaic.tif')
       inputFieldNames=getInputStructure()
       proj=referenceMetadata.get('proj')
       geotransform=referenceMetadata.get('gt') 
       #outputfile= 'F:\\Temporal\\Fg_Image.mosaiccacaca.tif' 
       outputfile= 'E:\\TSEB_Europe\\Output_Data\\averageMaps\\AverageLE_Month_' + month +'.mosaic.tif'
       
       Save_array_tiff(actualAverageLE, outputfile, proj,geotransform)            #Save the FG for the stations   
       print("saved")
   except:
       print("something went wrong")
   print("success") 

   try:                
       referenceMetadata=getSceneMetadata('Y:\\Gorka\\EUROPA\\MCD15A2.005\\MCD15A2.0052015129.mosaic.tif')
       inputFieldNames=getInputStructure()
       proj=referenceMetadata.get('proj')
       geotransform=referenceMetadata.get('gt') 
       #outputfile= 'F:\\Temporal\\Fg_Image.mosaiccacaca.tif' 
       outputfile= 'E:\\TSEB_Europe\\Output_Data\\averageMaps\\totalET_Month_' + month +'.mosaic.tif'
       
       Save_array_tiff(totalET, outputfile, proj,geotransform)            #Save the FG for the stations   
       print("saved")
   except:
       print("something went wrong")
   print("success")                   
           
def makeMapMonthOldVersion(dayFrom, dayTo, yearFrom, yearTo):
    
   #month = month number in format xx
   #yearFrom/yearTo = year in format xxxx doesn't include yearTo.
   import datetime 
   tm = datetime.datetime(2014, 1, 1) + datetime.timedelta(165 - 1)    
   correctDate = tm.strftime('%Y_%m_%d')
    
    
   yearFrom = 2003
   yearTo = 2008
   month = 02
   dayFrom = 1
   dayTo = 30
   count= 0
   zLength = dayTo - dayFrom + 1
       
   
   totalArrayH = np.full([3600,3600,zLength], np.nan)
   totalArrayLE = np.full([3600,3600, zLength], np.nan)
   hc = []
   for year in range(yearFrom, yearTo):
       #print month
       for day in range(dayFrom, dayTo):
           #print day
           
           print(count)
           print ('Y:\\Gorka\\Temporal\\InputsTSEB_'+str(year)+'_' + str(day)+'.nc', 'H_C')           
           try:
               
               print("here")
               month = 7
               print('X:\\TSEB_Europe\\Output_Data\\netcdf_results\\InputsTSEB_2003_0' + str(month) + '_' + str(day)+'.nc', 'H_C')
               H_C = readnetcdf('X:\\TSEB_Europe\\Output_Data\\netcdf_results\\InputsTSEB_2003_0' + str(month) + '_' + str(day)+'.nc', 'H_C')
               H_S = readnetcdf('X:\\TSEB_Europe\\Output_Data\\netcdf_results\\InputsTSEB_2003_0' + str(month) + '_' + str(day)+'.nc', 'H_S')
               LE_C = readnetcdf('X:\\TSEB_Europe\\Output_Data\\netcdf_results\\InputsTSEB_2003_0' + str(month) + '_' + str(day)+'.nc', 'LE_C')
               LE_S = readnetcdf('X:\\TSEB_Europe\\Output_Data\\netcdf_results\\InputsTSEB_2003_0' + str(month) + '_' + str(day)+'.nc', 'LE_S')
               print("there")
               #H_S[np.isnan(H_S)]=0
               #H_C[np.isnan(H_C)]=0
               #LE_S[np.isnan(LE_S)]=0
               #LE_C[np.isnan(LE_C)]=0
               
               #NEED CHECK TO SEE IF 1 NUMBER is NaN AND ONE ISNT. THEN ADDING TOGETHER GIVES US A NaN NO MATTER WHAT
               totalH = np.add(H_S, H_C)
               totalArrayH[:,:,count] = totalH
               
               
               #totalArrayH = np.add(totalArrayH, totalH)
               
               totalLE = np.add(LE_S, LE_C)
               totalArrayLE[:,:,count] = totalLE
               
               
               
               #totalArrayLE = np.add(totalArrayLE, totalLE)
#               hc = readnetcdf('Y:\\Gorka\\Temporal\\InputsTSEB_'+str(year)+'_'+str(month)+'_'+str(day), 'H_C')
#               hs = readnetcdf('Y:\\Gorka\\Temporal\\InputsTSEB_'+str(year)+'_'+str(month)+'_'+str(day), 'H_S')
#               lec = readnetcdf('Y:\\Gorka\\Temporal\\InputsTSEB_'+str(year)+'_'+str(month)+'_'+str(day), 'LE_C')  
#               les = readnetcdf('Y:\\Gorka\\Temporal\\InputsTSEB_'+str(year)+'_'+str(month)+'_'+str(day), 'LE_S')
               
               
           except:
               print('exceeded days probably')  
           count= count + 1
           
           
   #averageH = totalArrayH/count
   #averageLE = totalArrayLE/count
   
   #totalArrayLE[1000,1000,0] = np.nan
   #totalArrayLE[1000,1000,1] = 15
   #totalArrayLE[1000,1000,2] = 5         
    
   avgH = np.nanmean(totalArrayH, axis=2)
   avgLE = np.nanmean(totalArrayLE, axis=2)   

   #testLE = avgLE[1000,1000]      
   try:                
       referenceMetadata=getSceneMetadata('Y:\\Gorka\\EUROPA\\MCD15A2.005\\MCD15A2.0052015129.mosaic.tif')
       inputFieldNames=getInputStructure()
       proj=referenceMetadata.get('proj')
       geotransform=referenceMetadata.get('gt') 
       #outputfile= 'F:\\Temporal\\Fg_Image.mosaiccacaca.tif' 
       outputfile= 'E:\\TSEB_Europe\\Output_Data\\averageMaps\\averageH_Month6.mosaic.tif'
       
       Save_array_tiff(avgH, outputfile, proj,geotransform)            #Save the FG for the stations    
   except:
       print("something went wrong")
   print("success")

   try:                
       referenceMetadata=getSceneMetadata('Y:\\Gorka\\EUROPA\\MCD15A2.005\\MCD15A2.0052015129.mosaic.tif')
       inputFieldNames=getInputStructure()
       proj=referenceMetadata.get('proj')
       geotransform=referenceMetadata.get('gt') 
       #outputfile= 'F:\\Temporal\\Fg_Image.mosaiccacaca.tif' 
       outputfile= 'E:\\TSEB_Europe\\Output_Data\\averageMaps\\averageLE_Month6.mosaic.tif'
       
       Save_array_tiff(avgLE, outputfile, proj,geotransform)            #Save the FG for the stations    
   except:
       print("something went wrong")
   print("success")                     
     
#How to get around NaN issue. Maybe when we add arrays just replace NaN with 0 and then it's all fixed?
#op = np.zeros([4,4]) 
#ol = np.zeros([4,4]) 
#lk = np.zeros([4,4]) 
#count = np.zeros([4,4])
#
#op[1,1] = np.nan
#ol[1,1] = 5
#lk[1,1] = 10 
#
#count[op!=np.nan] += 1
#test = np.nansum(np.dstack((op, ol)),2)
#count[np.where(op>=0)] += 1
#
#
#m = np.zeros([4,4,3])
#m[:,:,0] = op
#m[:,:,1] = ol
#m[:,:,2] = lk
#avg = np.nanmean(m, axis=2)
#
#     
#
#total = np.add(op, ol)
#test = sum(sum(totalArray,[])) 
#dividedList = totalArray/10
#newList = [x / 7 for x in totalArray]
#sum(sum(x,[]))  
            