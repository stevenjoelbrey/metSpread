#!/usr/bin/env python

#------------------------------------------------------------------------------
# Description
#------------------------------------------------------------------------------
# Format and combine CMIP5 model data for, historical data, rcp4.5, and rcp 8.5
# CMIP5 files are of the form:
#
# "variable_Amon_ModelName_scenario_ensembleMember_YYYYMM-YYYYMM.nc"
# 
# ModelName = e.g. [ACCESS1.0, ACCESS1.3]
# variables = [hfls, tp, hurs, huss, mrso, si10, t2m, e]
# scenarios = [historical, rcp45, rcp85]
# ensembleMember = r1i1p1
# 
# Goals/tasks of this script, in order
# 1) mergetime : for like files (many models store data in lots of small files) 
# 2) seldate : for the merged time files. Subset to dates of interest, 198301-210012
# 3) remapbil : for the 198301-210012 merged and cut files, regrid to a common grid 

import sys
import os
import glob
import numpy as np
import cdo as cdo
import shutil



cdo = cdo.Cdo()


# Read command line arguments
args  = sys.argv
if len(args) > 1 :
	var = args[1]
else :
	var = 'tas'

# Directories to read and write from 
raw_data_dir = "/Users/sbrey/GoogleDrive/sharedProjects/metSpreadData/getCMIP5/r1i1p1_raw_downloaded_files/"
merged_time_dir = "/Users/sbrey/GoogleDrive/sharedProjects/metSpreadData/getCMIP5/r1i1p1_merged_time/"

# Other variables for building file names 
ensemble = 'r1i1p1'
scenarios = ['historical', 'rcp45', 'rcp85']


def get_var_unique_models(var, scenario, raw_data_dir) :
	"""
	Function for getting the unique model names that output a specified
	CMIP5 variable for a given scenario. 

	Parameters
	----------
		var : str, the variable to list files for
		scenario : str, the scenario to list for a given var
		raw_data_dir : str, the directory where these queries will be
		               made

	return
	------
		A list of models that output desired variable. This can be of 
		length zero if nothing matches the query

	"""

	# Get unique model names based on the available data for this variable
	all_files = glob.glob(os.path.join(raw_data_dir, var+"*"))
	# Make into single string
	all_files_string = " ".join(all_files)
	# remove annoyung raw_data_dir, and seperate files
	all_var_files = all_files_string.replace(raw_data_dir, "").split()

	# Extract model name from file name 
	var_model_name = []
	for i in range(len(all_var_files)) : 
		
		if all_var_files[i].find(scenario)!= -1 :
			# found the desired scenario
			var_model_name.append(all_var_files[i].split('_')[2])

	# Get the unique model names from the long list
	unique_var_model_names = np.unique(var_model_name)

	return unique_var_model_names


#------------------------------------------------------------------------------
# For a given var, and scenario, cary out tasks 1-3 for all models 
#------------------------------------------------------------------------------
var = 'tas'
scenario = 'historical'
var_available_models = get_var_unique_models(var, scenario, raw_data_dir)

for model in var_available_models :

	# variable_Amon_ModelName_scenario_ensembleMember_YYYYMM-YYYYMM.nc
	s = var + '_Amon_' + model + '_' + scenario + '_' + ensemble + "*"
	l = glob.glob(os.path.join(raw_data_dir, s))

	# Write the name of the merged file, this will be a combo of s and
	# the span of the dates for the detected files 
	minDate = l[0][-16:-10] 
	maxDate = l[-1][-9:-3]
	f_merged_time_out = s.replace('*','') + '_' + minDate + '-' + maxDate + '.nc'

	if len(l) == 1 :
		# To time merge required, the dates are represented by a single
		# file. Copy this file and place in desired merged_time directory 
		f_in = os.path.join(raw_data_dir, f_merged_time_out)
		f_out = os.path.join(merged_time_dir, f_merged_time_out)
		# do the copying 
		shutil.copy(f_in, f_out)

	else : 
		# Make a single string of files to merge, to pass to cdo command 
		files_to_merge = " ".join(l)

		# Make merged file name based on span of dates 
		f_out = os.path.join(merged_time_dir, f_merged_time_out)

		cdo.mergetime(input=files_to_merge, output=f_out, options="-b F64")






