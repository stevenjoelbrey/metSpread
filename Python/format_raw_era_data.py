#!/usr/bin/env python

#------------------------------------------------------------------------------
# Description
#------------------------------------------------------------------------------
# Format and combine raw monthly means of analysis and forecast fields from 
# era-interim. 
# NOTE: This will not include soil moisture, which is more tricky and is handled
# NOTE: with a script all of its own. 
# TODO: Comment on the differences required for analysis and forecast fields. 


import sys
import os
import glob
import numpy as np
import cdo as cdo

# Read command line arguments
args  = sys.argv
if len(args) > 1 :
	var   = str(args[1])
	year1 = int(args[2])
	year2 = int(args[3])
else :
	year1 = 1983
	year2 = 2017

years = np.arange(year1, year2+1)

# Directories required for cdo commands and output files 
dataDir = "/Users/sbrey/GoogleDrive/sharedProjects/metSpreadData/ERA-Interim"
time_merge_out = "/Users/sbrey/GoogleDrive/sharedProjects/metSpreadData/ERA-Interim/merged_time"
common_grid_dir = "/Users/sbrey/GoogleDrive/sharedProjects/metSpreadData/ERA-Interim/merged_t_COMMON_GRID"
common_grid_txt = os.path.join(dataDir,"COMMON_GRID.txt")
# Gives access to CDO commands 
cdo = cdo.Cdo()

forecast_vars = ['e', 'tp', 'slhf']
analysis_vars = ['t2m', 'si10', 'd2m']

# Combine yearly files for forecast fields, as there is nothing to do with different
# forecast accumulation periods for individual months within the yearly files. 
for an_var in analysis_vars :

	# Merge yearly files 
	f_in = glob.glob(os.path.join(dataDir, an_var+"*"))
	f_out = os.path.join(time_merge_out, an_var+"_"+str(year1)+"_"+str(year2)+".nc")
	cdo.mergetime(input=" ".join(f_in), output=f_out, options="-b F64")

	# Now regrid the merged t file to the common grid 
	f_in_common = f_out
	f_out_common = os.path.join(common_grid_dir, an_var+"_"+str(year1)+"_"+str(year2)+".nc")
	cdo.remapbil(common_grid_txt, input=f_in_common, output=f_out_common, options="-b F64")


# TODO: Make common grid version of merged_time files 