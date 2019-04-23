import numpy as np
from netCDF4 import Dataset
import sys
import matplotlib.pyplot as plt
import os
import pandas as pd # TODO: Lots of training
import glob as glob


def get_all_model_names(dataDir="../Data/CMIP5/r1i1p1_rcp_COMMON_GRID") :
    """
    Gets all the unique CMIP5 model names available in the dataDir 
    directory. Returned as an array.  
    """
    all_model_names = []
    all_files = glob.glob(dataDir+"/*.nc")#+"/tas_Amon*")
    #print(all_files)
    for f in all_files :
        all_model_names.append(f.split("/")[4].split("_")[2])
    return np.unique(all_model_names)

# Get CMIP5 estimates for the same time period to make comparisons
def get_cmip5_nc(var, rcp, model, spatial_mask, spatial_mean=False, inspect=False):
    """
    This function will be for loading a particular CMIP5 NetCDF file, which will be spatially subset 
    by this function. These CMIP5 model output to be loaded have been regridded to the "Common grid" 
    using $ cdo remapbil. Defualts to taking the spatial mean of all global data for the var argument. 
    Change area loaded using the minX, maxX, etc. arguments. 
    
    Parameters
    ----------
        var : str, The CMIP5 variable name to be loaded. File names match variable names. 
        rcp : str, "45" or "85", refers to representative concentration pathway. 
        model : The name of the model that created the var, to be loaded from. 
        spatial_mask : numpy array, where equal to 1 (True) are locations to mask, meaning
                       they will NOT be used in calculations, spatial means, etc. 
        spatial_mean : Boolean, if False (default) the data  will be returned on a 
                       [t, lon, lat] grid. 
    
    Return
    ------
        valsCut : array[month, lat, lon] for the chosen spatial extent or 
                  array[month] spatial mean for the chosen spatial extent. 
        t_mon : array of months as pd.date_range object. 
        lon : array of longitude values that were used aftering trimming the global data
        lat : array of latitude values that were used after trimming the global data
    
    """
        
    if (var == 'mrso') or () or (var =='mrlsl.integrated') or (var == "lai"):
        domain = 'Lmon'
    else :
        domain = 'Amon'
        
    # Create link to the monthly file in this projects directory structure. 
    dataDir = os.path.join(".." ,"Data" ,"CMIP5" ,'r1i1p1_rcp_COMMON_GRID')
    f = var + "_" + domain + "_" + model + "_rcp" + rcp + "_r1i1p1_198301-210012.nc"
    
    if inspect :
        print(f)
        
    loadFile = os.path.join(dataDir, f)
    
    # Check to see if the file exists! Data logs indicate that not all requests exist
    if(not os.path.isfile(loadFile)):
        raise ValueError(f+ " File does not exist")
    
    # Load the nc data
    nc = Dataset(loadFile)
    vals = nc.variables[var][:]
    lon = nc.variables["lon"][:]
    lat = nc.variables["lat"][:]
    
    # Create a new array, where each time slice will be a masked version 
    # of the data loaded from the nc file 
    vals_masked = np.ma.empty(shape=vals.shape)
    
    # Pandas handling of time so all models have the exact same origin and such. 
    # TODO: handle required changes to this section for when historical CMIP5 data
    # TODO: are also used. 
    t = nc.variables["time"]
    if(len(t) == 1416) :
        # Convert to pandas time array, on the assumption t[0]=2006-01 & t[-1]=2100-12
        t_mon = pd.date_range("1983-01-01", periods=len(t), freq="M")
    else:
        raise ValueError('Error in number of months for file: '+ f + " 1416 expected.")
        
    # Now subset the data based on the passed spatial_mask
    for i in range(len(t)) : 
        vals_masked[i,:,:] = np.ma.masked_array(vals[i,:,:], mask=spatial_mask, copy=True)
    
    # Now, take the mean value in this spatial domain. Note, the 
    # same is done to the era-interim data using the method:
    # get_era_nc_vals()
    if spatial_mean :
        
        # There is some extra stuff in here to make sure that
        # values of np.nan for soil moisture also get masked. 
        # i.e., in short, this masked all np.nan values too 
        # when time means of the division masked values are taken. 
        
        spatial_mask_3d = vals_masked.mask
        nan_mask = np.isnan(vals_masked.data)
        combined_mask = spatial_mask | nan_mask
        vals_masked_new = np.ma.array(vals_masked, mask=combined_mask)
        
        spatial_mean_values = np.ma.mean(vals_masked_new, axis=(1,2))
        vals_masked = spatial_mean_values.data # if you do not do this, forces dtype object in dataframe. 
                    
        if False :
            
            print("Spatial Mask types")
            print(type(spatial_mask_3d))
            print(spatial_mask_3d.shape)
            print("nan_mask")
            print(type(nan_mask))
            print(nan_mask.shape)
            print("combined_mask")
            print(type(combined_mask))
            print(combined_mask.shape)
            print("vals_masked_new")
            print(vals_masked_new.shape)
        
    nc.close()
    
    return vals_masked, t_mon, lon, lat

def get_era_nc_variable(var, spatial_mask, spatialMean=False, startYear=1983):
    """
    This function will be for loading a particular nc file, which will be spatially subset.
    The data loaded are from the merged_t_COMMON_GRID directory. These data have been
    regridded from thier native resolution using cdo remapbil. 
        
    Parameters
    ----------
        var : str, The variable (and file name) of the ECMWF era-interim data to be 
              loaded. 
        spatial_mask : numpy array, where equal to 0 are locations to keep data for and 
                       where equal to 1 (True) are grid points to mask (ignore data). 
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
    f = var + "_1983-2017.nc" 
    loadFile = os.path.join(dataDir, f)
    
    # Check to see if the file exists! 
    if(not os.path.isfile(loadFile)):
        raise ValueError(f + " File does not exist")
    
    # Load the nc data
    nc = Dataset(loadFile)
    vals = nc.variables[var][:]
    lon = nc.variables["lon"][:]
    lat = nc.variables["lat"][:]
    
    # Create a new array, where each time slice will be a masked version 
    # of the data loaded from the nc file 
    vals_masked = np.ma.empty(shape=vals.shape)
    
    # Pandas handling of time so all models have the exact same origin and such. 
    t = nc.variables["time"]
    if(len(t) == 420):
        # Convert to pandas time array, on the assumption t[0]=startYear-01-01 the data
        # being monthly. Defualt files are 1983-2017.  
        t_mon = pd.date_range(str(startYear)+"-01-01", periods=len(t), freq="M")
    else:
        # there should always be 420 months in these data (35 years of data times 12)
        raise ValueError('Error in number of months for file: '+ f + " 420 expected.")
        
    # Now subset the data based on the passed spatial_mask
    for i in range(len(t)) : 
        vals_masked[i,:,:] = np.ma.masked_array(vals[i,:,:], mask=spatial_mask, copy=True)
    
    # take a spatial mean? 
    if spatialMean:
        # Masked entries are ignored
        vals_masked = np.ma.mean(vals_masked, axis=(1,2))
    
    return vals_masked, t_mon, lon, lat


# The next three functions were developed and copied from make_era_and_CMIP5_units_identical.ipynb
# on Feb 14 2019. 
def m_to_mass_flux(tp_m) :
    """
    This function takes mean total precipitation in meters [m] and converts to a mass flux. 
    https://confluence.ecmwf.int/pages/viewpage.action?pageId=65218804
    Incoming tp_m argument must have units of mean m/day. Multiplication by the density
    of water, 1 day per 86400 seconds are used to obtain units
    of kg m**-2 s**-1, the so-called tp_prime value that is returned. 
    Will work anytime meters of water per month needs to be converted to kg/m**2/s
    """
    density = 1000. # kg m**-3
    seconds_per_day = 86400. 
    tp_prime = tp_m * density / seconds_per_day 

    return tp_prime

def e_to_mass_flux(e_m) :
    """
    Converts evaporation from m/day to monthly mean mass flux in units of kg m**-2 s**-1. 
    Multiplies era-interim values by -1 to account for difference in sign convention for ERA 
    and CMIP5. 
    """
    density = 1000. # kg m**-3
    seconds_per_day = 86400. 
    e_prime = e_m * density / seconds_per_day * -1
    
    return e_prime

def J_flux_to_W_flux(slhf) :
    """
    https://confluence.ecmwf.int/pages/viewpage.action?pageId=65218804
    Converts a flux of J m**-2 month**-1 to W m**-2 month**-1
    """
    seconds_per_day = 86400. 
    hfls_watts = slhf / seconds_per_day  * -1 # temp factor of dividing by 2. 
    
    return hfls_watts


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