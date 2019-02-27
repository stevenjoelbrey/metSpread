#!/usr/bin/env python

#------------------------------------------------------------------------------
# Description
#------------------------------------------------------------------------------
# The purpose of this script is to combine era-interim soil layers 1-4 
# into a variable the describes the total mass of water per square meter of soil. 
# The yearly soil files were downloaded with get_ecmwf.py. They are labeled 
# soil_YYYY.nc. They contain four soil layers, swvl1, swvl2, swvl3, swvl4. We 
# want to combine these layers because CMIP5 tracks soil moisture in all layers (mrso).
# ERA soil units are m**3/m**3, this needs to be converted to kg/m**2 in order
# to match th units of CMIP5 soil moisture. The sketch below is to aid the 
# unit conversion performed by this script.  
#
#    ({ }  Atmosphere ({ })
#        1m  
#     *-------*         Z0 = 0.0 m 
# 1m /|      /|         Z1 = 0.00 - 0.07 m 
#   / |     / | 2.89m   Z2 = 0.07 - 0.28 m
#  *--|----*  |         Z3 = 0.28 - 1.00 m 
#  |  *----|--*         Z4 = 1.00 - 2.89 m 
#  | /     | /      
#  *-------*      
# 
# Soil depths for layers 1-4 are Balsamo et al. 2009 but seen in Albergel et al 2012. 
# https://journals.ametsoc.org/doi/10.1175/2008JHM1068.1
#
# ERA gives us the fraction of a volume of soil of depth Z that is water in some
# phase (liquid or ice). The depth of this volume is 2.89 m, and the unit area of 
# interest is 1 m**2. 
# The total volume of this volume is 2.89 m * 1 m**2 = 2.89 m**3. To get the mass 
# in this cube multiply by the density of liquid water, which will have some error, 
# since some of (but not a known amount) this volume could be ice. 
 
import sys
import os
import glob
import numpy as np
from netCDF4 import Dataset
import cdo as cdo

# Read command line arguments
args  = sys.argv
if len(args) > 1 :
	year1 = int(args[1])
	year2 = int(args[2])
else :
	year1 = 1983
	year2 = 2017

# The soil moisture data live here: 
dataDir = "/Users/sbrey/GoogleDrive/sharedProjects/metSpreadData/ERA-Interim"
time_merge_out = "/Users/sbrey/GoogleDrive/sharedProjects/metSpreadData/ERA-Interim/merged_time"
common_grid_dir = "/Users/sbrey/GoogleDrive/sharedProjects/metSpreadData/ERA-Interim/merged_t_COMMON_GRID"
common_grid_txt = os.path.join(dataDir,"COMMON_GRID.txt")

# For access to handy $ cdo commands via python  
cdo = cdo.Cdo()

# Handle all years 
years = np.arange(year1, year2+1)

# from 1st to 4th (meters)
layer_depths = [0.07, 0.21, 0.72, 1.89]
print("Total depth of soil layers %f" % np.sum(layer_depths) )

# Make each soil_YYYY.nc file into a mrso like file 
for y in years : 

	#--------------------------------------------------------------------------
	# Open the nc file, scale each layers volume of water by the depth of that
	# layer
	#--------------------------------------------------------------------------
	print("Making mrso file for %i year soil file" %y)

	nc_file = os.path.join(dataDir, "soil_" + str(y) + ".nc")
	nc = Dataset(nc_file, "r") 
	
	# Get nc attributes to pass along to a new file written here. 

	# Scale volume by the depth of the layer. This gives the data a dimension of
	# depth
	swvl1 = nc.variables['swvl1'][:] * layer_depths[0]
	swvl2 = nc.variables['swvl2'][:] * layer_depths[1]
	swvl3 = nc.variables['swvl3'][:] * layer_depths[2]
	swvl4 = nc.variables['swvl4'][:] * layer_depths[3]

	# new variable, that is the sum of the total volume depths of water
	total_water_volume = swvl1 + swvl2 + swvl3 + swvl4 # (m water)
	water_kg_per_area = total_water_volume * 1000. # (m water) * (1000 kg/m**3) = (kg/m**2)

	#--------------------------------------------------------------------------
	# Write this as a new nc file
	#--------------------------------------------------------------------------
	outputFile = os.path.join(dataDir, "mrlsl.integrated_" + str(y) + ".nc")

	nc_out = Dataset(outputFile, 'w', format='NETCDF4')
	nc_out.description = 'Soil moisture (water+ice) per unit area'
	nc_out.location = 'Global'
	nc_out.createDimension('time',  len(nc.variables["time"]) )
	nc_out.createDimension('latitude', len(nc.variables["latitude"]) )
	nc_out.createDimension('longitude', len(nc.variables["longitude"]) )

	VAR_ = nc_out.createVariable("mrlsl.integrated", 'f4',('time', 'latitude','longitude'))
	VAR_.long_name = "Total Soil Moisture Content to 2.89 m"
	VAR_.units = "kg m-2"

	# Create time variable
	time_ = nc_out.createVariable('time', 'i4', ('time',))
	time_.units = nc.variables['time'].units
	time_.calendar = "gregorian"

	# create lat variable
	latitude_ = nc_out.createVariable('latitude', 'f4', ('latitude',))
	latitude_.units = nc.variables['latitude'].units

	# create longitude variable
	longitude_ = nc_out.createVariable('longitude', 'f4', ('longitude',))
	longitude_.units = nc.variables['longitude'].units

	# Write the actual data to these dimensions
	VAR_[:]       = water_kg_per_area
	time_[:]      = nc.variables["time"][:]
	latitude_[:]  = nc.variables["latitude"][:]
	longitude_[:] = nc.variables["longitude"][:]

	nc_out.close() # File we just wrote
	nc.close()     # File where the layers were stored separately 


#------------------------------------------------------------------------------
# cdo commands for merging yearly nc files and regridding to the common grid
#------------------------------------------------------------------------------
# Combine the yearly mrso files that we just wrote into a single merged time
# file
print("Combining and regridding files using cdo()")

f_in = glob.glob(os.path.join(dataDir, "mrlsl.integrated_*"))
f_out = os.path.join(time_merge_out, "mrlsl.integrated_"+str(year1)+"_"+str(year2)+".nc")
cdo.mergetime(input=" ".join(f_in), output=f_out, options="-b F64")

# Now, create a version of that file that lives on the common grid 
f_out_common_grid = os.path.join(common_grid_dir, "mrlsl.integrated_"+str(year1)+"-"+str(year2)+".nc")
cdo.remapbil(common_grid_txt, input=f_out, output=f_out_common_grid, options="-b F64")

print("Script ececuted without error")
print("Woot woot! Its the small things. There are enough tears in graduate school.")