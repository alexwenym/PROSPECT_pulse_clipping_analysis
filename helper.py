import ROOT
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import uproot
import h5py
import os

def get_xy_from_TGraph(graph): 
    
    n = graph.GetN()
    
    x = np.asarray([])
    y = np.asarray([])

    for i in range(n): 
        x = np.append(x, graph.GetPointX(i))
        y = np.append(y, graph.GetPointY(i))
        
    return x, y


def waveform_integral(xarray,yarray, intmin, intmax): 
    
    xarraycut = np.asarray([xarray[i] for i in range(len(xarray)) if intmin<xarray[i]<intmax])
    yarraycut = np.asarray([yarray[i] for i in range(len(xarray)) if intmin<xarray[i]<intmax])
    
    ymin_extra = np.interp(intmin, xarray, yarray)
    ymax_extra = np.interp(intmax, xarray, yarray)
    
    xarraycut = np.concatenate([[intmin],xarraycut,[intmax]])
    yarraycut = np.concatenate([[ymin_extra],yarraycut,[ymax_extra]])
    
    return np.trapz(yarraycut, x=xarraycut)

def get_time(xarray, yarray):
    
    half_max = np.amax(yarray)/2
    
    xarraycut = np.asarray([xarray[i] for i in range(len(xarray)) if -20<xarray[i]<20])
    yarraycut = np.asarray([yarray[i] for i in range(len(xarray)) if -20<xarray[i]<20])
    
    return np.interp(half_max, yarraycut, xarraycut, left=float('NaN'), right=float('NaN')) 

def get_time_index(xarray, yarray): 
    
    half_max = np.amax(yarray)/2
    maxindex = np.argmax(yarray)
    
    xarraycut = xarray[0:maxindex+1]
    yarraycut = yarray[0:maxindex+1]
    
    half_max_index = np.where(yarraycut-half_max > 0, yarraycut-half_max, np.inf).argmin()
    
    return half_max_index
    

def get_total_area(xarray,yarray): 
    
    minx = np.amin(xarray)
    maxx = np.amax(xarray)
    
    return waveform_integral(xarray,yarray, intmin=minx, intmax=maxx) 

def get_total_area_discretepts(xarray,yarray): 
    
    maxindex = np.argmax(yarray)
    #maxindex = get_time_index(xarray, yarray)
    
    beginindex = maxindex - 3
    endindex = maxindex + 25
    
    xarrayint = xarray[beginindex:endindex]
    yarrayint = yarray[beginindex:endindex]
    
    return np.trapz(yarrayint, x=xarrayint)
      

def get_PSD(xarray, yarray): 
    
    time = get_time(xarray, yarray)
    
    total = waveform_integral(xarray,yarray, intmin=-12+time, intmax=200+time)
    tail  =  waveform_integral(xarray,yarray, intmin=44+time, intmax=200+time)
    
    return tail/total

def get_PSD_discretepts(xarray, yarray):
    
    timeindex = get_time_index(xarray, yarray) 
    
    totalstart = -3
    totalend = 49
    tailstart = 10
    tailend = totalend
    
    total = np.trapz(yarray[timeindex+totalstart:timeindex+totalend], x=xarray[timeindex+totalstart:timeindex+totalend])
    #tail = np.trapz(yarray[timeindex+10:timeindex+50], x=xarray[timeindex+10:timeindex+50])
    tail = np.trapz(yarray[timeindex+tailstart:timeindex+tailend], x=xarray[timeindex+tailstart:timeindex+tailend])
    
    return tail/total
    
    
    






    