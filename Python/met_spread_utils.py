import numpy as np
from netCDF4 import Dataset
import sys
import matplotlib.pyplot as plt
import os
import pandas as pd # TODO: Lots of training

def get_cmip5_nc(var="tas", rcp="45", model="ACCESS1-0", minX=0., maxX=360., minY=-90., maxY=90., spatial_mean=True):
    """
    This function will be for loading a particular CMIP5 NetCDF file, which will be spatially subset 
    by this function. These CMIP5 model output to be loaded have been regridded to the "Common grid" 
    using $ cdo remapbil. Defualts to taking the spatial mean of all global data for the var argument. 
    Change area loaded using the minX, maxX, etc. arguments. 

    # TODO: Needs to be modified to handle historical data when it becomes available. 
    
    Parameters
    ----------
	    var : str, The CMIP5 variable name to be loaded. File names match variable names. 
	    rcp : str, "45" or "85", refers to representative concentration pathway. 
	    model : The name of the model that created the var
	    minX : float, min longitude (0-360) of the data to return. 
	    maxX : float, max longitude (0-360) of the data to return. 
	    minY : float, min latitude of the data to return. 
	    maxY : float, max latitude of the data to return. 
	    spatial_mean : Boolean, if True (default) a spatial mean of the era-interim data
	                   is not taken and the data are returned on a t,lon,lat grid. 
    
    Return
    ------
		valsCut : array[month, lat, lon] for the chosen spatial extent or 
		          array[month] spatial mean for the chosen spatial extent. 
		t_mon : array of months as pd.date_range object. 
		lonCut : array of longitude values that were used aftering trimming the global data
		latCut : array of latitude values that were used after trimming the global data
    
    """
        
    # Create link to the monthly file in this projects directory structure. 
    dataDir = os.path.join(".." ,"Data" ,"CMIP5" ,'r1i1p1_rcp45_rcp85_merged_t_COMMON_GRID')
    f = var + "_" + "Amon_" + model + "_rcp" + rcp + "_r1i1p1_200601-210012.nc"
    loadFile = os.path.join(dataDir, f)
    
    # Check to see if the file exists! Data logs indicate that not all requests exist
    if(not os.path.isfile(loadFile)):
        raise ValueError(f+ " File does not exist")
    
    # Load the nc data
    nc = Dataset(loadFile)
    vals = nc.variables[var][:]
    lon = nc.variables["lon"][:]
    lat = nc.variables["lat"][:]
    
    # Pandas handling of time so all models have the exact same origin and such. 
    # TODO: handle required changes to this section for when historical CMIP5 data
    # TODO: are also used. 
    t = nc.variables["time"]
    if(len(t) == 1140):
        # Convert to pandas time array, on the assumption t[0]=2006-01 & t[-1]=2100-12
        t_mon = pd.date_range("2006-01-01", periods=len(t), freq="M")
    else:
        raise ValueError('Error in number of months for file: '+ f + " 1140 expected.")
        
    # Now subset the data based on the passed max and min values for lon and lat
    lonIndex = np.where( ((lon >= minX) & (lon <= maxX)) )[0]
    latIndex = np.where( ( (lat >= minY) & (lat <= maxY) ) )[0]
    timeIndex = range(len(t)) # because we want all months, for now

    # Subset the 2D field
    lonCut = lon[lonIndex]
    latCut = lat[latIndex]

    # Subset the 3D field
    valsCut = vals[np.ix_(timeIndex, latIndex, lonIndex)]
    
    # Now, take the mean value in this spatial domain. Note, the 
    # same is done to the era-interim data using the method:
    # get_era_nc_vals()
    if spatial_mean :
        spatial_mean_values = np.mean(valsCut, axis=(1,2))
        valsCut = spatial_mean_values
    
    return valsCut, t_mon, lonCut, latCut


def get_era_nc_vals(var="t2m", minX=0., maxX=360., minY=-90., maxY=90., spatialMean=False, startYear=1992):
    """
    This function will be for loading a particular nc file, which will be spatially subset.
    The data loaded are from the merged_t_COMMON_GRID directory. These data have been
    regridded from thier native resolution using cdo remapbil. 
    
    TODO: Add "endyear" as an argument in addition to "startyear" as new MTBS fire data is 
    TODO: going to require going further back in time to train the model. 
    
    Parameters
    ----------
	    var : str, The variable (and file name) of the ECMWF era-interim data to be 
	          loaded. 
	    minX : float, min longitude (0-360) of the data to return. 
	    maxX : float, max longitude (0-360) of the data to return. 
	    minY : float, min latitude of the data to return. 
	    maxY : float, max latitude of the data to return. 
	    
	    spatialMean : Boolean, if False (default) a spatial mean of the era-interim data
	                  is not taken and the data are returned on a t,lon,lat grid. 
                  
    Return
    ------
	    valsCut : The chosen "var" as a ndarray(t, lat, lon) or if spatialMean = True
	              ndarray(t). 
	    t_monCut : pd.date_range describing the t axis of valsCut. 
	    lonCut : longitude ndarray
	    latCut : latitude ndarray
    
    """
    
    # Create link to the monthly file
    dataDir = os.path.join("..", "Data", "ERA-INTERIM", "merged_t_COMMON_GRID")
    f = var + "_1990-2015.nc" 

    loadFile = os.path.join(dataDir, f)
    
    # Check to see if the file exists! 
    if(not os.path.isfile(loadFile)):
        raise ValueError(f + " File does not exist")
    
    # Load the nc data
    nc = Dataset(loadFile)
    vals = nc.variables[var][:]
    lon = nc.variables["lon"][:]
    lat = nc.variables["lat"][:]
    
    # Pandas handling of time so all models have the exact same origin and such. 
    t = nc.variables["time"]
    if(len(t) == 312):
        # Convert to pandas time array, on the assumption t[0]=1990-01-01 and t[-1]=2015-12-01
        # 26 years times 12 months per year makes for 312 months expected in these arrays
        t_mon = pd.date_range("1990-01-01", periods=len(t), freq="M")
    else:
        # there should always be 312 months in these data (26 years of data)
        raise ValueError('Error in number of months for file: '+ f + " 312 expected.")
        
    # Now subset the data based on the passed max and min values for lon and lat
    lonIndex = np.where( ((lon >= minX) & (lon <= maxX)) )[0]
    latIndex = np.where( ( (lat >= minY) & (lat <= maxY) ) )[0]
    
    # Where is the first month of startYear
    t_first = np.where(t_mon.year == startYear)[0][0]
    
    # because we want all months, for now, but starting in startYear
    timeIndex = range(t_first, len(t)) 

    # Subset the 2D field
    lonCut = lon[lonIndex]
    latCut = lat[latIndex]
    t_monCut = t_mon[timeIndex]

    # Subset the 3D field
    valsCut = vals[np.ix_(timeIndex, latIndex, lonIndex)]

    # take a spatial mean? 
    if spatialMean:
        valsCut = np.mean(valsCut, axis=(1,2))
    
    return valsCut, t_monCut, lonCut, latCut