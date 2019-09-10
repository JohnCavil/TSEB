# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 13:46:20 2017

@author: gmg, hlan

Gets the FG data just for the stations. Modified from previous FG script that uses full 3600x3600 arrays.
"""

import pandas
import glob
import numpy as np
import gdal
from itertools import compress
from MOD11_HDF_READ import getSceneMetadata
from MOD11_HDF_READ import Save_array_tiff
from MOD11_HDF_READ import getInputStructure
#loads the excel file with thefluxnet sites
fileName='Y:\\Gorka\\EUROPA\\FluxnetSitesCoordinates\\FluxnetSitesCoordinates.xlsx'
print("starting")
def readLAIArray(array,row,col,ID):
    print('starting')
    import netCDF4
    import numpy as np

    LAI = array[row,col,:]   
   
    return LAI
#    Rn_sw_soil = fh.variables['Rn_sw_soil'][row,col]   
#    Rn_sw_veg = fh.variables['Rn_sw_veg'][row,col]  
    
    
    
ExcelFile=pandas.read_excel(fileName,sheetname=0)
Rows=ExcelFile.RowPos
Cols=ExcelFile.ColPos
IGBP=ExcelFile.IGBP
StationID=ExcelFile.SITE_ID
StationName=ExcelFile.SITE_NAME


#opens the tif file that contains the mean maps of LAI europe for each day of the year.


image='Y:\\Gorka\\Temporal_2\\MeanEuropeanLAI_mosaic.tif'

fid=gdal.Open(image,gdal.GA_ReadOnly)
array=np.zeros([fid.RasterXSize,fid.RasterYSize,fid.RasterCount])   #make an array with the size of the   
count=0

print (fid.RasterCount)

for band in range(fid.RasterCount):
    array[:,:,count]=fid.GetRasterBand(count+1).ReadAsArray()   
    count=count+1    
    print('doing this')
    
    
###########################################################################################################################################  
    
# Reads the IGBP calssification map.
IGBP_Class_Map='Y:\\Gorka\\EUROPA\\MCD12Q1\\MCD12Q12008-01-01'
landcover=np.fromfile(IGBP_Class_Map,dtype='uint8')
landcover=np.reshape(landcover, [3600,3600])*1.0
landcover[landcover>=17]=np.nan
landcover[landcover==0]=np.nan
 # IGBP classification
#                       Everg   Everg.  Decid.  Decid.                                                                          Crop
#                       Needle  Broad   Needle  Broad   Mixed   Closed  Open    Woody           Grass   Wet     Crop            Veg     Snow
#                 Water Forest  Forest  Forest  Forest  Forest  shrubs  shrubs  Savanna Savanna land    lands   lands  Urban    Mosaic  Ice     Baren   
landCoverClasses= [0,   1,      2,      3,      4,      5,      6,      7,      8,      9,      10,     11,     12,    13,      14,     15,     16 ];  # land cover class, from MOD12Q1 
globalMaxLaiArray=np.zeros([landcover.shape[0],landcover.shape[1]])  #map of max LAI for each landcover

globalMaxLaiArray[landcover==0]=0
GEvNeFor=4
GEvBrFor=4
GDeNeFor=4
GDeBrFor=4
GMixForr=4
GCloShr =4 #5 
GOpeShr=4 #5
GWooSav=4 #5
GSav=4      #4
GGrass=4 #5
GWetLand=4 #5
GCropLand  =2 #10
GUrban=4 #5
GCropVegMosaic=4 #6
GSnow=4 #0
GBaren=4 #2
globalMaxLaiArray[landcover==0]=0
globalMaxLaiArray[landcover==1]=GEvNeFor
globalMaxLaiArray[landcover==2]=GEvBrFor
globalMaxLaiArray[landcover==3]=GDeNeFor
globalMaxLaiArray[landcover==4]=GDeBrFor
globalMaxLaiArray[landcover==5]=GMixForr
globalMaxLaiArray[landcover==6]=GCloShr  
globalMaxLaiArray[landcover==7]=GOpeShr
globalMaxLaiArray[landcover==8]=GWooSav
globalMaxLaiArray[landcover==9]=GSav
globalMaxLaiArray[landcover==10]=GGrass
globalMaxLaiArray[landcover==11]=GWetLand
globalMaxLaiArray[landcover==12]=GCropLand  
globalMaxLaiArray[landcover==13]=GUrban
globalMaxLaiArray[landcover==14]=GCropVegMosaic
globalMaxLaiArray[landcover==15]=GSnow
globalMaxLaiArray[landcover==16]=GBaren
globalMaxLaiArray=globalMaxLaiArray*10. 

#############################################################################################################################
   
#creates an array with the number of stations and the doys to place the extractions
output_LAI=np.zeros([46,len(StationID)]) #creates array of size 46 x number of stations (60)
output_LAI_Stations=np.zeros([len(StationID),46])
globalMaxLaiArray_Out_put=np.zeros([len(StationID)])
LandCoverVector=np.zeros([len(StationID)])  #creates array of 60
ID_all= [''  for x in xrange(len(Rows))]
IGBP_all= [''  for x in xrange(len(Rows))]
i=0


for Station in range(0,len(Rows)):
    Row=Rows[Station]
    Col=Cols[Station]
    ID=StationID[Station]
    IGBP_Class=IGBP[Station]
    LAI=readLAIArray(array,Row,Col,ID)
    LAI[LAI==255]=np.nan
    output_LAI_Stations[i,:]=LAI
    output_LAI[:,i]=LAI
    globalMaxLaiArray_Out_put[i]=globalMaxLaiArray[Row,Col]
    #print(Row,Col,LAI,ID)
    ID_all[i]=ID
    IGBP_all[i]=IGBP_Class
    
    LandCoverVector[i]=landcover[Row,Col]
    
#    output_LAI[:,i]=LAI
    i=i+1
    
    
#Removes the stations with all nans in case there are some. LAI_Clean = array of all stations with their LAI over time (in size 365/8=45)

#LAI_Clean=output_LAI[:,~np.all(np.isnan(output_LAI), axis=0)]
LAI_Clean=output_LAI
LAI_Clean[np.isnan(LAI_Clean)]=0

#LAI_Clean2=output_LAI
#LAI_Clean2[np.isnan(LAI_Clean2)]=0

ID_Clean=[~np.all(np.isnan(output_LAI), axis=0)]
COmpressedData=list(compress(ID_all,ID_Clean[0]))


np.save('Y:\\Gorka\\EUROPA\\FluxnetSitesCoordinates\\output_LAI_stations.npy',output_LAI_Stations)
np.save('Y:\\Gorka\\EUROPA\\FluxnetSitesCoordinates\\output_LandCover_stations.npy',LandCoverVector)
np.save('Y:\\Gorka\\EUROPA\\FluxnetSitesCoordinates\\output_globalMaxLaiArray.npy',globalMaxLaiArray_Out_put)


output_LAI_Stations=np.load('Y:\\Gorka\\EUROPA\\FluxnetSitesCoordinates\\output_LAI_stations.npy')
LandCoverVector=np.load('Y:\\Gorka\\EUROPA\\FluxnetSitesCoordinates\\output_LandCover_stations.npy')
globalMaxLaiArray_Out_put=np.load('Y:\\Gorka\\EUROPA\\FluxnetSitesCoordinates\\output_globalMaxLaiArray.npy')


#clean the landcover array
LandCover_Stations_Clean=LandCoverVector[~np.all(np.isnan(output_LAI), axis=0)]

#clean the max LAI for the stations
globalMaxLaiArray_Out_put_Clean=globalMaxLaiArray_Out_put[~np.all(np.isnan(output_LAI), axis=0)]

globalMaxLaiArray_Out_put_Clean[12]=50 #This is just a one-off case because station 12 has water as its landcover by mistake

#Creates Fg Array
fgArray=np.zeros([array.shape[0],array.shape[1],array.shape[2]])
fgArray_stations=np.zeros([LAI_Clean.shape[0],LAI_Clean.shape[1]])

#fgArray_stations=np.zeros([7,10])
"""
Search for the maximum value of the LAI and gets its possition. There are two possible situations
in which the maximum is in the first part of the year, therefore it searches for the minimum in the second.
And the opposite, in which the max is in the second part and therfore the minimun needs to be found in the
first part of the year.
"""
#LAI_max_pos=np.argmax(LAI_Clean, axis=0) #finds maximum in horizonal columns (for each station)
#LAI_min_pos=np.argmin(LAI_Clean, axis=0) #finds minimum in horizonal columns (for each station)

#45 = 365/8? What if 46?
maxIntervalPosition=(0,45)
minxIntervalPosition=(0,45)        

maxArray=np.amax(array[:,:,maxIntervalPosition[0]:maxIntervalPosition[1]],axis=2)
minArray=np.amin(array[:,:,minxIntervalPosition[0]:minxIntervalPosition[1]],axis=2)
maxPos=np.argmax(array[:,:,maxIntervalPosition[0]:maxIntervalPosition[1]], axis=2)
minPos=np.argmin(array[:,:,minxIntervalPosition[0]:minxIntervalPosition[1]], axis=2)
#maxPos=maxPos+maxIntervalPosition[0]
#minPos=minPos+minxIntervalPosition[0]

#IDs=np.zeros([LAI_Clean.shape[0],LAI_Clean.shape[1]])
#IDs2=np.zeros([LAI_Clean.shape[0],LAI_Clean.shape[1]])

#IDs=np.zeros(LAI_Clean.shape[1])
#IDs2=np.zeros(LAI_Clean.shape[1])

IDs=np.zeros((3600, 3600))
IDs2=np.zeros((3600, 3600))

BreakArray=np.zeros([array.shape[0],array.shape[1]])

IDs[maxPos>minPos]=1
IDs2[maxPos<=minPos]=1

#IDs[LAI_max_pos>LAI_min_pos]=1
#IDs2[LAI_max_pos<LAI_min_pos]=1

#IDs[np.isnan(LandCover_Stations_Clean)]=0 #remove the stations that dont have landcover type.
#IDs2[np.isnan(LandCover_Stations_Clean)]=0





#                                                              ------THIS IS THE CODE FOR THE STATION VERSION-------
#-------------------------------------------------------------------------------------------------------------------------------------------

count=-1                                                            #'count' keeps track of station number
for x in range(0,array.shape[0]):
    for y in range(0,array.shape[1]): 
        FIndPos=None
        count=count+1
        
        if IDs[x,y]==1:                                                  #if station==1 then it means that the maximum LAI position occurs after the minimum
            if ~np.isnan(landcover[x,y])== True:                         #checking that a landcover exists for the station
                #print 'landcover is good'
                #initialising variables
              
                minPos_vec=None
                maxPos_vec=None
                LAI_i_plus_one=None
                LAI_minimum=None
                LAI_i=None
            
                
                
                minPos_vec=minPos[x,y]
                maxPos_vec=maxPos[x,y]
                
                
                for i in range(minPos_vec,maxPos_vec-1):        #now we iterate through the positions between the minimum and maximum LAI for that station
                    
    
                    LAI_i=array[x,y,minPos_vec-1]/10.
                    LAI_i_plus_one=array[x,y,i+1]/10.
                    
                    max_LAI_Pixel=maxArray[x,y]/10.  #does nothing?
                    test1=LAI_i_plus_one-LAI_i
                    test2=100*((LAI_i_plus_one-LAI_i)/LAI_i)
                    
                    #print(i)
                    #print(count)
                    
                    if FIndPos==None:
                        if(test2>20) & (test1>0.2):                      #tests whether we have increased 0.2 between the minimum LAI and whereever we're at. Also checks if we have increased 20%.
                            FIndPos=i                                    #finds the postion for where we need to break. We break just before i+1. So one time step before the tests are actually fulfilled.
                            BreakArray[x,y]=FIndPos                        
                            test1=None
                            test2=None
                            
                          
    #                print('the LAI minimum is: ') 
    #                print(LAI_minimum)
    #                print('------------------------------------------')
    #                print('the LAI at point is: ') 
    #                print(LAI_i_plus_one)
    #                print('------------------------------------------')
                    
        if IDs2[x,y]==1:
            if ~np.isnan(landcover[x,y])== True:
                #print 'landcover is true'
                minPos_vec=None
                LAI_i=None
                LAI_i_plus_one=None
            
                minPos_vec=minPos[x,y]
                max2Pos=np.argmax(array[x,y,minPos_vec:45], axis=0)
                
                for i in range(minPos_vec, 46):                      #iterate from the minimum to the end of the year.
                    
                    LAI_i=array[x,y,minPos_vec]/10.
                    if i <= 44:                                          #if we're not near the end of the year   
                        
                        LAI_i_plus_one=array[x,y,i+1]/10.
                    else:                                                #if we're in time step 45 or 46:
    
                        LAI_i_plus_one=array[x,y,0]/10.
                        max_LAI_Pixel=maxArray[x,y]/10.
                        test1=LAI_i_plus_one-LAI_i
                        test2=100*((LAI_i_plus_one-LAI_i)/LAI_i)
                    
                    if FIndPos==None:                       
                        if(test2>20) & (test1>0.2):
                            
                            FIndPos=i
                            BreakArray[x,y]=FIndPos
                            test1=None
                            test2=None
                            
                if FIndPos==None:                                        #if we still havent found the position then it means that the rbeakpoint (FindPos) is somewhere between time step 0 and the maximum LAI.
                    count_from_zero=0
                    for i in range(0,maxPos_vec-1):    #iterate through from 0 to the maximum LAI point.
                        
                        
                        LAI_i=array[x,y,minPos_vec]/10.
                        LAI_i_plus_one=array[x,y,i+1]/10.
                        max_LAI_Pixel=maxArray[x,y]/10.
                        test1=LAI_i_plus_one-LAI_i
                        test2=100*((LAI_i_plus_one-LAI_i)/LAI_i)
                        #print(count)
                        #print(LAI_minimum)
                        #print(LAI_plus_one)
                        #print('tests')
                        #print(test1)
                        #print(test2)
                        
                        if FIndPos==None:
                           
                           
                           if(test2>20) & (test1>0.2):
                               
                               FIndPos=i
                               BreakArray[x,y]=FIndPos
                                
                                #FIndPos=None
                               test1=None
                               test2=None
                
               
        #print(FindPos)  
        #print(x,y,test1,test2,FIndPos)       
    
              
        for step in range(0,array.shape[2]):                       #iterates through the timesteps (the x axis in LAI_Clean array)
        
            tst = landcover[x,y]                            
            
            fgArray[x,y,step]=(array[x,y,step])/(globalMaxLaiArray[x,y])
             #sets the FG value for the stations according to equation 13
            if (tst!=1 and tst!=3 and tst!=14 and tst!=16):    
                #print ("runin")
                #print BreakArray[34, 1515]
                #print FIndPos
                #print maxPos[34, 1515]
                if (~np.isnan(BreakArray[x,y])):             #test if the landcovers need to use eq 14. These 4 only need to use eq 13.
                        
                        
                   if (~np.isnan(BreakArray[x,y])):
                       
                       if (int(BreakArray[x,y])<maxPos[x,y]):         #need if statement for landcover here. Landcovers that are not dependant on seasons should should only use eq 13.
                    #            print('count is:')
                    #            print(count)
                    #            print('the position is:')
                    #            print(FindPos)
                              #sets the FG value for the stations according to equation 14
                           
                            for i in range(int(BreakArray[x,y]),maxPos[x,y]):
                                fgArray[x,y,i]=(maxArray[x,y]/globalMaxLaiArray[x,y])*(1-np.exp(-2*array[x,y,i]))
                            
                       if (int(BreakArray[x,y])>=maxPos[x,y]):                                  #if the position we found is after the max LAI, then we need to only iterate until time step 46 and not to the max LAI point like before.
                                    
                            #print('breakArray is larger than max' )        
                            for i in range(int(BreakArray[x,y]),46):
                                fgArray[x,y,i]=(maxArray[x,y]/globalMaxLaiArray[x,y])*(1-np.exp(-2*array[x,y,i]))  
                            for i in range(0,maxPos[x,y]):
                                fgArray[x,y,i]=(maxArray[x,y]/globalMaxLaiArray[x,y])*(1-np.exp(-2**array[x,y,i]))                         
                                    
    print x
    print "need 3600"  
    timeTest = np.zeros([5,5])
    timeTest[4,4] = x
    Path = 'E:\\HenrikTSEB\\Fg_stations_only_output\\'
    
    np.save(Path + 'fgArray' , timeTest)                                  
#if ~np.isnan(landcover[2,2968])== True:
#    print 'landcover is true'
#testID = IDs2[34,1515]  
#testID  = 1  
##testmaxpos = maxPos[4,3033]
##testminpos = minPos[4,3211]
#testLandcover = landcover[34, 1515]
#testarraySlice = array[34,1515, :]
#testidddd = IDs2[1012,818]    
#testid = np.argmax(IDs2)
##
##
##fgArray[fgArray==np.inf]= np.nan
#testfgArraySlice2 = fgArray[34,1515,:]
#testglobalMax = globalMaxLaiArray[4, 3033]
#arraySlice[np.isnan(fgArraySlice)] = 0
    
#countx=-1
#county=-1
#countz=-1
#
#
#q1=0
#q2=0
#q3=0
#    
#fgArray_stations_as_map=np.zeros([7,10,46]) #fgArray_stations_as_map=np.zeros([8,7,46])
#for station in fgArray_stations:                #makes the 46x59 fgArray_stations array into a 8x7x46 fgArray_stations_as_map array to mimic a map.
#    countx=countx+1
#    county=county+1
#    countz=countz+1
#    
#    fgArray_stations[fgArray_stations>1]=1      #sets the fg value to 1 if it's over 1 since that's not possible.
#    
#    
#    wolo=-1
#    for i in range(0, 10):
#        
#       
#        #countx=countx+1
#        q1=q1+1
#
#        for n in range(0, 7):
#           q2=q2+1 
#           #county=county+1
#           wolo=wolo+1
#           
#          
#           
#           for m in range(0, 46): 
#               #countz=countz+1
#               q3=q3+1
#               
#               if(wolo<61):
#                   fgArray_stations_as_map[n,i,m] = fgArray_stations[m, wolo]
#               
#               
              

               
               
    


print "now saving"
 
fgArray[fgArray==np.inf]=np.nan           
fgArray[fgArray>1]=1  


try:                
    referenceMetadata=getSceneMetadata('Y:\\Gorka\\EUROPA\\MCD15A2.005\\MCD15A2.0052015129.mosaic.tif')
    inputFieldNames=getInputStructure()
    proj=referenceMetadata.get('proj')
    geotransform=referenceMetadata.get('gt') 
    #outputfile= 'F:\\Temporal\\Fg_Image.mosaiccacaca.tif' 
    outputfile= 'E:\\HenrikTSEB\\Fg_stations_only_output\\Fg_ImageEuropeMapCorrectFfourCtwo.mosaic.tif'
    Save_array_tiff(fgArray, outputfile, proj,geotransform)            #Save the FG for the stations    
except:
    print "something went wrong"
print "success"
Path = 'E:\\HenrikTSEB\\Fg_stations_only_output\\'
    
np.save(Path + 'fgArray' , fgArray)            


#MaxCrop = np.zeros([3600,3600])
#MaxCrop = maxArray[np.where(landcover==12)]  
#
#MaxForest = np.zeros([3600,3600])
#
#
#MaxForest = maxArray[np.where(landcover==1)] 
#MaxForest2 = maxArray[np.where(landcover==2)]
#MaxForest3 = maxArray[np.where(landcover==3)] 
#MaxForest4 = maxArray[np.where(landcover==4)] 
#MaxForest5 = maxArray[np.where(landcover==5)] 
# 
#
#maxForestTotal = []
#maxForestTotal.extend(MaxForest)
#maxForestTotal.extend(MaxForest2)
#maxForestTotal.extend(MaxForest3)
#maxForestTotal.extend(MaxForest4)
#maxForestTotal.extend(MaxForest5)


         
#referenceMetadata=getSceneMetadata('Y:\\Gorka\\EUROPA\\MCD15A2.005\\MCD15A2.0052015129.mosaic.tif')
#inputFieldNames=getInputStructure()
#proj=referenceMetadata.get('proj')
#geotransform=referenceMetadata.get('gt') 
##outputfile= 'F:\\Temporal\\Fg_Image.mosaiccacaca.tif' 
#outputfile= 'E:\HenrikTSEB\Fg_stations_only_output\\Fg_ImageEuropeMap.mosaic.tif'
#Save_array_tiff(fgArray_stations, outputfile, proj,geotransform) 

#Save_array_tiff(fgArray, outputfile, proj,geotransform)           
                
                
# -----------------------------------------------------------------------------------------------------------------------------------------               
                
                
              
 #                                                      -----THIS IS THE CODE FOR THE PIXEL VERSION-------               
#
#for x in range(0,LAI_Clean.shape[0]):
#    for y in range(0,LAI_Clean.shape[1]):
#for x in range(1000,1010):
#     for y in range(1000,1010):
#for x in range(0,array.shape[0]):
#    for y in range(0,array.shape[1]):             
#            
#            
#    
#    
#        FIndPos=None
#        if IDs[x,y]==1:
#            if ~np.isnan(landcover[x,y])== True:
#                minPos_vec=None
#                LAI_i=None
#                LAI_i_plus_one=None
#                test=-100
#                
#                minPos_vec=minPos[x,y]
#                maxPos_vec=maxPos[x,y]
#                for i in range(minPos_vec,maxPos_vec-1):
#                    LAI_i=array[x,y,minPos_vec-1]/10.
#                    LAI_i_plus_one=array[x,y,i+1]/10.
#                    max_LAI_Pixel=maxArray[x,y]/10.  #does nothing?
#                    test1=LAI_i_plus_one-LAI_i
#                    test2=100*((LAI_i_plus_one-LAI_i)/LAI_i)
#                           
#                    if FIndPos==None:
#                        if (test2>20) & (test1>0.2):
#                            FIndPos=i
#                            BreakArray[x,y]=FIndPos
#                            test1=None
#                            test2=None
#        if IDs2[x,y]==1:
#            if ~np.isnan(landcover[x,y])== True:
#                minPos_vec=None
#                LAI_i=None
#                LAI_i_plus_one=None
#                
#                minPos_vec=minPos[x,y]
#                max2Pos=np.argmax(array[x,y,minPos_vec:45], axis=0)
#                for i in range(minPos_vec,46):
#                    LAI_i=array[x,y,minPos_vec]/10.
#                    if i <=44:
#                        LAI_i_plus_one=array[x,y,i+1]/10.
#                    else:
#                        LAI_i_plus_one=array[x,y,0]/10.
#                    max_LAI_Pixel=maxArray[x,y]/10.
#                    test1=LAI_i_plus_one-LAI_i
#                    test2=100*((LAI_i_plus_one-LAI_i)/LAI_i)
#                    if FIndPos==None:
#                        if (test2>20) & (test1>0.2):
#                            FIndPos=i
#                            BreakArray[x,y]=FIndPos
#                            test1=None
#                            test2=None
#                if FIndPos==None:
#                    for i in range(0,maxPos_vec-1):
#                        LAI_i=array[x,y,minPos_vec]/10.
#                        LAI_i_plus_one=array[x,y,i+1]/10.
#                        max_LAI_Pixel=maxArray[x,y]/10.
#                        test1=LAI_i_plus_one-LAI_i
#                        test2=100*((LAI_i_plus_one-LAI_i)/LAI_i)
#                        
#                        if FIndPos==None:
#                            if (test2>20) & (test1>0.2):
#                            
#                                FIndPos=i
#                                BreakArray[x,y]=FIndPos
#                            
#                            #FIndPos=None
#                                test1=None
#                                test2=None
#                            
#                                #print('found ID 2',x,y,test)
#                                #FIndPos=None
#                                #test=None
##                except:
##                    minPos_vec=0
#        #print(x,y,test1,test2,FIndPos)
#        for step in range(0,array.shape[2]):
#            fgArray[x,y,step]=(array[x,y,step])/(globalMaxLaiArray[x,y])
#        if (~np.isnan(BreakArray[x,y])):
#            if (int(BreakArray[x,y])<maxPos[x,y]): 
#                for i in range(int(BreakArray[x,y]),maxPos[x,y]):
#                    fgArray[x,y,i]=(maxArray[x,y]/globalMaxLaiArray[x,y])*(1-np.exp(-2*array[x,y,i]))    
#            if (int(BreakArray[x,y])>=maxPos[x,y]):
#                for i in range(int(BreakArray[x,y]),46):
#                    fgArray[x,y,i]=(maxArray[x,y]/globalMaxLaiArray[x,y])*(1-np.exp(-2*array[x,y,i]))  
#                for i in range(0,maxPos[x,y]):
#                    fgArray[x,y,i]=(maxArray[x,y]/globalMaxLaiArray[x,y])*(1-np.exp(-2**array[x,y,i])) 
#                    
# ---------------------------------------------------------------------------------------------------------------------------                   
#                    
#                for i in range(int(BreakArray[x,y]),46):
#                    fgArray[x,y,i]=(maxArray[x,y]/globalMaxLaiArray[x,y])*(1-np.exp(-2*array[x,y,i]))
#                for i in range(0,maxPos_vec):
#                    fgArray[x,y,i]=(maxArray[x,y]/globalMaxLaiArray[x,y])*(1-np.exp(-2*array[x,y,i]))
##                for step in range(0,array.shape[2]):
###                    fgArray[x,y,step]=(array[x,y,step])/(maxArray[x,y])
###                for i in range(int(BreakArray[x,y]),45):
###                    if landcover[x,y] !=0:
###                        fgArray[x,y,i]=(maxArray[x,y]/globalMaxLaiArray[x,y])*(1-np.exp(-2*array[x,y,i]))
###                    if fgArray[x,y,i]>1:
###                        fgArray[x,y,i]=1
###                for i in range(0,maxPos_vec):
###                    if landcover[x,y] !=0:
###                        fgArray[x,y,i]=(maxArray[x,y]/globalMaxLaiArray[x,y])*(1-np.exp(-2*array[x,y,i]))
###                    if fgArray[x,y,i]>1:
###                        fgArray[x,y,i]=1
###                else:
###                    fgArray[x,y,i]
###                BreakArray[x,y]=minPos
###                for step in range(0,array.shape[2]):
###                    fgArray[x,y,step]=(array[x,y,step])/(maxArray[x,y])
###                for i in range(0,int(maxPos[x,y])):
###                    if landcover[x,y] !=0:
###                        fgArray[x,y,i]=(maxArray[x,y]/globalMaxLaiArray[x,y])*(1-np.exp(-2*array[x,y,i]))
###                    if fgArray[x,y,i]>1:
###                        fgArray[x,y,i]=1
###                for i in range(int(BreakArray[x,y]),45):
###                    if landcover[x,y] !=0:
###                        fgArray[x,y,i]=(maxArray[x,y]/globalMaxLaiArray[x,y])*(1-np.exp(-2*array[x,y,i]))
###                    if fgArray[x,y,i]>1:
###                        fgArray[x,y,i]=1
###                else:
###                    fgArray[x,y,i]
####                else:
###                    BreakArray[x,y]=1000
##               # print(x,y,minPos,maxPos[x,y],'yujuuuuu')
##
###for Station in range(0,len(Rows)):
###    Row=Rows[Station]
###    Col=Cols[Station]
###    ID=StationID[Station]
###    IGBP_Class=IGBP[Station]
###    LAI=readLAIArray(array,Row,Col,ID)
###    LAI[LAI==255]=np.nan
###    output_LAI_Stations[0,i,:]=LAI
###    output_LAI[:,i]=LAI
###    globalMaxLaiArray_Out_put[i]=globalMaxLaiArray[Row,Col]
###    print(Row,Col,LAI,ID)
###    ID_all[i]=ID
###    IGBP_all[i]=IGBP_Class
###    #    output_LAI[i,ii]=LAI
###    i=i+1           
#
##referenceMetadata=getSceneMetadata('Y:\\Gorka\\EUROPA\\MCD15A2.005\\MCD15A2.0052015129.mosaic.tif')
##inputFieldNames=getInputStructure()
##proj=referenceMetadata.get('proj')
##geotransform=referenceMetadata.get('gt') 
###outputfile= 'F:\\Temporal\\Fg_Image.mosaiccacaca.tif' 
##outputfile= 'E:\HenrikTSEB\Fg_stations_only_output\\Fg_Image.mosaic.tif'
##Save_array_tiff(fgArray_stations, outputfile, proj,geotransform) 
#
##Save_array_tiff(fgArray, outputfile, proj,geotransform)
#
#
#
#
####FInd the thirdposition that is necesary
###BreakArray=np.zeros([array.shape[0],array.shape[1]])
###for x in range(0,array.shape[0]):
###    for y in range(0,array.shape[1]):
###        if array[x,y,0] !=0:
###            test=-100
###            for i in range(int(minPos[x,y]),int(maxPos[x,y])):
###                if test < 20:
###                    LAI=array[x,y,i+1]
###                    test=((LAI-minArray[x,y])/minArray[x,y]*100)
###                    #print(test)
###                    if test >=20:
###                        BreakArray[x,y]=i
###
###BreakArray[BreakArray==0]=np.nan
###THis part of the code creates two plots
##import matplotlib.pyplot as plt
##import numpy as np
##plt.figure(figsize=(10,10))
##plt.imshow(LAI_Clean, interpolation='nearest',cmap=plt.cm.jet)
###plt.xticks(np.arange(0,len(COmpressedData)), COmpressedData,rotation=90)
##plt.xticks(np.arange(0,len(ID_all)), ID_all,rotation=90)
##plt.colorbar()
###plt.yticks(np.arange(0,5), ['F', 'G', 'H', 'I', 'J'])
##plt.savefig('E:\\SPACE\\Figures_Papers\\LAI_Of_EC_Stations_all.png')
##plt.show()
##
##
##
##Ids=list(compress(IGBP_all,ID_Clean[0]))
##
##import matplotlib.pyplot as plt
##import numpy as np
##plt.figure(figsize=(15,10))
##plt.imshow(LAI_Clean, interpolation='nearest',cmap=plt.cm.jet)
##plt.xticks(np.arange(0,len(Ids)), Ids,rotation=90)
##plt.colorbar()
##
##plt.savefig('E:\\SPACE\\Figures_Papers\\LAI_Of_EC_StationsIGBP.png')
##plt.show()
##
##
