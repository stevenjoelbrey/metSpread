#!/usr/bin/env python

#------------------------------------------------------------------------------
# Description
#------------------------------------------------------------------------------
# Format and combine raw monthly means of analysis and forecast fields from 
# era-interim. 
#
# NOTE: This will not include soil moisture, which is more tricky and is handled
# NOTE: with a script all of its own (Python/format_era_soil_moisture.py). 
#
# The two types of data that are downloaded are "Monthly Means of Daily Means" and 
# "Monthly Means of Daily Forecast Accumulations". 
#
# The difference between these data can be found here:
# http://apps.ecmwf.int/datasets/data/interim-full-moda/levtype=sfc/
#
# For forecast fields, they provide the forecast time period average for a variable, e.g.
# total precipitation (tp). "Step 0-12: averages of the precipitation produced 
# from the first 12 hours of each 00:00 forecast, plus the precipitation produced 
# from the first 12 hours of each 12:00 forecast. The resulting data is temporally 
# continuous."
# 
# https://confluence.ecmwf.int/pages/viewpage.action?pageId=65218804
#
# The yearly monthly mean files need to be combined into a single file
# for the years of interest. That can be done using
#
# $ cdo -b F64 mergtime e* e_year1_year2.nc 

import sys
import os
import glob
import numpy as np
import cdo as cdo

# Read command line arguments
args  = sys.argv
if len(args) > 1 :
	year1 = int(args[2])
	year2 = int(args[3])
else :
	print("Using hard-coded years")
	year1 = 1983
	year2 = 2017

# Execute required commands and operations for all years of interest. 
years = np.arange(year1, year2+1)
print("Combining data for years", years)

# Directories required for cdo commands and output files 
dataDir = "/Users/sbrey/GoogleDrive/sharedProjects/metSpreadData/ERA-Interim"
time_merge_out = "/Users/sbrey/GoogleDrive/sharedProjects/metSpreadData/ERA-Interim/merged_time"
common_grid_dir = "/Users/sbrey/GoogleDrive/sharedProjects/metSpreadData/ERA-Interim/merged_t_COMMON_GRID"
common_grid_txt = os.path.join(dataDir,"COMMON_GRID.txt")

# Gives access to CDO commands 
cdo = cdo.Cdo()

# Variables to combine data for
analysis_vars = ['t2m', 'si10', 'd2m', 'sp']
forecast_vars = ['e', 'tp', 'slhf', 'ro','uvb', 'strd', 'ssrd']


print("Handling the analysis fields.", analysis_vars)
# Combine yearly files for analysis fields, as there is nothing to do with different
# forecast accumulation periods for individual months within the yearly files. 
for an_var in analysis_vars :

	# Merge yearly files 
	f_in = glob.glob(os.path.join(dataDir, an_var+"*"))
	f_out = os.path.join(time_merge_out, an_var+"_"+str(year1)+"-"+str(year2)+".nc")
	cdo.mergetime(input=" ".join(f_in), output=f_out, options="-b F64")

	# Now regrid the merged t file to the common grid 
	f_in_common = f_out
	f_out_common = os.path.join(common_grid_dir, an_var+"_"+str(year1)+"-"+str(year2)+".nc")
	cdo.remapbil(common_grid_txt, input=f_in_common, output=f_out_common, options="-b F64")


print("Handling the forecasted fields.", forecast_vars)
# Combine yearly files for forecast fields, forecast accumulation periods for individual 
# months within the yearly files. 
for fc_var in forecast_vars :

	# Merge the yearly files into a single combined file
	f_in = glob.glob(os.path.join(dataDir, fc_var+"_*"))
	f_out_combine = os.path.join(time_merge_out, fc_var+"_"+str(year1)+"-"+str(year2)+".nc")
	cdo.mergetime(input=" ".join(f_in), output=f_out_combine, options="-b F64")

	# Now, regrid this merged file to the common grid 
	f_out_common_grid = os.path.join(common_grid_dir, fc_var+"_"+str(year1)+"-"+str(year2)+".nc")
	cdo.remapbil(common_grid_txt, input=f_out_combine, output=f_out_common_grid, options="-b F64")


print("Script ececuted without error")
print("Woot woot!")
