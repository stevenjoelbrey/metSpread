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

makeHistory = True

import sys
import os
import glob
import numpy as np
import cdo as cdo
import shutil
from netCDF4 import Dataset

cdo = cdo.Cdo()


# Directories to read and write from 
raw_data_dir = "/Users/sbrey/GoogleDrive/sharedProjects/metSpreadData/getCMIP5/r1i1p1_raw_downloaded_files/"
merged_time_dir = "/Users/sbrey/GoogleDrive/sharedProjects/metSpreadData/getCMIP5/r1i1p1_merged_time/"
cut_time_dir = "/Users/sbrey/GoogleDrive/sharedProjects/metSpreadData/getCMIP5/r1i1p1_merged_time_cut/"


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


def	make_var_files(var, scenario, raw_data_dir) : 
	"""
	Merges individual files for a given variable, and scenario, for all models
	that have output for those. 

	Parameters
	----------
		var : str, the variable to list files for
		scenario : str, the scenario to list for a given var
		raw_data_dir : str, the directory where these queries will be
		               made

	return
	------
		None, desired output are written as netCDF files. 

	"""

	# Get the names of the models available for this query
	var_available_models = get_var_unique_models(var, scenario, raw_data_dir)
	
	# loop through these models, merging and moving files as needed 
	for model in var_available_models :

		# e.g. variable_Amon_ModelName_scenario_ensembleMember_YYYYMM-YYYYMM.nc
		if var == "mrso" : 
			# Land parameter
			time_span = '_Lmon_'
		else :
			# atmosphere parameter
			time_span = '_Amon_'

		# Create a string that will list disired files
		s = var + time_span + model + '_' + scenario + '_' + ensemble + "*"
		l = glob.glob(os.path.join(raw_data_dir, s))

		# Write the name of the merged file, this will be a combo of s and
		# the span of the dates for the detected files 
		minDate = l[0][-16:-10] # All nc files end with dates spanned. Can be indexed
		maxDate = l[-1][-9:-3]  # from the end of the charactor string. 
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

		# Cut the newly merged files to the dates of interest only
		if scenario == 'historical' : 
			# Handle the desired cut dates for historical files as well as special
			# errors messages relevant to these cutoffs
			if int(maxDate) < int(200512) :
				print("------------------------------------------------------------------")
				print("ERROR- model last date not late enough")
				print("WRONG DATES for " + s)
				print("The max date in the file was " + maxDate)
				print("This file requires a manual addition of first month of an RCP file")
				print("------------------------------------------------------------------")

			if int(minDate) > int(198301) :
				print("------------------------------------------------------------------")
				print("ERROR- model file data do not have earlier enough start date")
				print("WRONG DATES for " + s)
				print("The min date in the file was " + minDate)
				print("This file requires getting more historical data for this variable")
				print("------------------------------------------------------------------")

			f_seldate = s.replace('*','') + '_198301_200512.nc'
			f_seldate_out = os.path.join(cut_time_dir, f_seldate) 
			cdo.seldate('1983-01-01,2005-12-31', input=f_out, output=f_seldate_out, options="-b F64")

			# Check for the length of this cut historical file written above.
			# 1983 - 2005 is 23 years, times 12 months per year, the time dimension
			# of the cut files should be length 276
			nc = Dataset(f_seldate_out, 'r')
			t = nc.variables['time'][:]
			nc.close()
			if len(t) != 276 :
				print("------------------------------------------------------------------")
				print("ERROR- merged and cut date file does not have correct # of months")
				print("The correct number of months should be 276, got %i " % len(t) )
				print(s)
				print("------------------------------------------------------------------")

		# TODO: RCP 8.5 statements
		# TODO: RCP 4.5 statements  

#------------------------------------------------------------------------------
# For a given var, and scenario, cary out tasks 1-3 for all models 
#------------------------------------------------------------------------------
if makeHistory == True :
	make_var_files('sfcWind', 'historical', raw_data_dir); print("sfcWind complete")
	make_var_files('tas', 'historical', raw_data_dir); print("tas complete")
	make_var_files('mrso', 'historical', raw_data_dir); print("mrso complete")
	make_var_files('huss', 'historical', raw_data_dir); print("huss complete")
	make_var_files('pr', 'historical', raw_data_dir); print("pr complete")
	make_var_files('hurs', 'historical', raw_data_dir); print("hurs complete")
	make_var_files('hfls', 'historical', raw_data_dir); print("hfls complete")

print("Script executed without fail")
print("yay")
