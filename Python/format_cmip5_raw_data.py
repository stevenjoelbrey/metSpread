#!/usr/bin/env python

#------------------------------------------------------------------------------
# Description
#------------------------------------------------------------------------------
# Format and combine CMIP5 model data for, historical data, rcp4.5, and rcp 8.5
# CMIP5 files are of the form:
#
# "variable_Amon_ModelName_scenario_ensembleMember_YYYYMM-YYYYMM.nc"
# 
# 
# Goals/tasks of this script, in order
# 1) mergetime : for like files (many models store data in lots of small files) 
#    it makes sense to merge variables into a single file. In this script this task
#    is completed by the make_var_files() method. These files are output into the 
#    directory r1i1p1_merged_time.
# 2) Cut merged historical files and rcp files to the desired date range. This uses
#    the data in r1i1p1_merged_time and cuts historical files to 1983-2005, puts them
#    in r1i1p1_history_cut and rcp files to 2006-2100 and puts them in r1i1p1_rcp_cut. 
#    This is performed by the make_var_files() method. 
# 3) Finally, we want to pair history files to rcp files, because this is one continious
#    set of model output for a given model, variable, scenario. This is done with the
#    pair_history_to_rcp() method. This file lists the history file for which a merged
#    rcp file match could not be found. It is worth making sure that these data do
#    not exist on https://esgf-node.llnl.gov/, or it may otherwise mean that the data
#    were never downloaded, or there is something wrong in steps 1 &| 2.
# 4) remapbil : for the 198301-210012 merged and cut files, they are regridded to a 
#    common grid using cdo.remapbil. These are the output to actually be used in this 
#    work. This regridding is only done with the 198301-210012 merged files have the
#    correct number of months. 
#
# NOTE: For optimal mrlsl_integrated availability, run interpolate_mrlsl first. 
#
#
# Issues raised by this code, box fixed and unresolved are documented in the file
# "format_cmip5_raw_data_print_statements.txt"
# created via $ python format_cmip5_raw_data.py > ../Data/CMIP5/format_cmip5_raw_data_print_statements.txt


# TODO: Write pandas dataframes to indicate what data are missing, 
# TODO: or where operations failed. 

# NOTE: HADGEM-ES has a very strange (late) history period. 

makeHistory = False # Runs all history file manipulations
makeNew     = False # Runs all file initial trimming and merging
makePairs   = False # Pairs available cut files to make single 1986-2100 file and common grid

import sys
import os
import glob
import numpy as np
import cdo as cdo
import shutil
from netCDF4 import Dataset
import pandas as pd

# Set up text file to save all print statements to
#orig_stdout = sys.stdout
#f = open('out.txt', 'w')
#sys.stdout = f

cdo = cdo.Cdo()

# Directories to read and write from 
base_dir = "../../metSpreadData/getCMIP5"
raw_data_dir = os.path.join(base_dir, "r1i1p1_raw_downloaded_files/")
merged_time_dir = os.path.join(base_dir, "r1i1p1_merged_time/")
common_grid_file = os.path.join(base_dir, 'COMMON_GRID.txt')

# Other variables for building file names 
ensemble = 'r1i1p1'

def get_var_unique_models(var, scenario, search_data_dir) :
	"""
	Function for getting the unique model names that output a specified
	CMIP5 variable for a given scenario. 

	Parameters
	----------
		var : str, the variable to list files for
		scenario : str, the scenario to list for a given var
		search_data_dir : str, the directory where these queries will be
		                  made

	return
	------
		A list of models that output desired variable. This can be of 
		length zero if nothing matches the query

	"""

	# Get unique model names based on the available data for this variable
	all_files = glob.glob(os.path.join(search_data_dir, var + "_*"))
	# Make into single string
	all_files_string = " ".join(all_files)
	# remove annoyung search_data_dir, and seperate files
	all_var_files = all_files_string.replace(search_data_dir, "").split()

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
	that have output for those. Files of a like type are merged with cdo.mergetime.
	That file is written. That file is then read and cut to the desired date range,
	which is 1983-01-01 to 2005-12-31 for history and 2006-01-01,2100-12-31 for 
	RCP 4.5 and 8.5. Wouldn't it be cool to live in an RCP4.5 world? 
	Why am I looking at these?

	Parameters
	----------
		var : str, the variable to list files for
		scenario : str, the scenario (historical, rcp45, rcp85) to list for a given var
		raw_data_dir : str, the directory where these queries will be
		               made

	return
	------
		None, desired output are written as netCDF files. 
		They should be a file of all dates merged and a file cut to 
		dates specified above. 

	"""

	# Get the names of the models available for this query
	var_available_models = get_var_unique_models(var, scenario, raw_data_dir)
	
	# loop through these models, merging and moving files as needed 
	for model in var_available_models :

		# e.g. variable_Amon_ModelName_scenario_ensembleMember_YYYYMM-YYYYMM.nc
		# TODO: Handle name of top soil variable here. 
		if var == "mrso" or var == "mrlsl.integrated" : 
			# Land parameter
			time_span = '_Lmon_'
		else :
			# atmosphere parameter
			time_span = '_Amon_'

		# Create a string that will list disired files
		s = var + time_span + model + '_' + scenario + '_' + ensemble + "*" # YYYYMM-YYYYMM.nc
		l = glob.glob(os.path.join(raw_data_dir, s))

		print("start--------------------------------------------")
		print("Working on query:")
		print(model)
		print(s)
		print("There are %i files to merge for query: " % len(l) )

		if len(l) == 0 :

			print("***************************************")
			print("No files found for inquery related to:")
			print(s)
			print("***************************************")

		else : 

			# There are data, proceed. ---------------------------------------------

			# Write the name of the merged file, this will be a combo of s and
			# the span of the dates for the detected files 
			minDate = l[0][-16:-10] # All nc files end with dates spanned. Can be indexed
			maxDate = l[-1][-9:-3]  # from the end of the charactor string. 
			f_merged_time_out = s.replace('*','') + '_' + minDate + '-' + maxDate + '.nc'

			if len(l) == 1 :

				# No time merge required, the dates are represented by a single
				# file. Copy this file and place in desired merged_time directory 
				print("Single file being coppied")
				f_in = os.path.join(raw_data_dir, f_merged_time_out)
				f_out = os.path.join(merged_time_dir, f_merged_time_out)
				# do the copying 
				shutil.copy(f_in, f_out)

			else : 

				print("Merging the %i files associated with query" % len(l))
				# Make a single string of files to merge, to pass to cdo command 
				files_to_merge = " ".join(l)

				# Make merged file name based on span of dates 
				# NOTE: 'f_merged_time_out' has actual span seen by file names 
				f_out = os.path.join(merged_time_dir, f_merged_time_out)

				cdo.mergetime(input=files_to_merge, output=f_out, options="-b F64")

			#######################################################################
			# Cut the newly merged files to the dates of interest only. scenarios
			# handled differently. 
			#######################################################################
			
			# TODO: Consider using CDO.mergetime on all history files and rcp of a 
			# TODO: type in one step. Would that solve the missing month problem? 

			# Historical ----------------------------------------------------------
			if scenario == 'historical' : 
				# Handle the desired cut dates for historical files as well as special
				# errors messages relevant to these cutoffs
				cut_time_dir = os.path.join(base_dir, "r1i1p1_history_cut/")

				# maxDate from r1i1p1_raw_downloaded_files for this type
				if int(maxDate) < int(200512) :
					# TODO: HadGem will often throw an error here. 
					# TODO: When this is the case, simply append the require set of months from
					# TODO: an RCP file. 
					print("***************************************")
					print("ERROR- model last date not late enough")
					print("WRONG DATES for " + s)
					print("The max date in the file was " + maxDate)
					print("This file requires a manual addition of first month of an RCP file")
					print("***************************************")
					hist_date_2 = maxDate
				else :
					# If the max date of the history files was bigger than date below, they 
					# will be subset so that the date below is the last date, so label the
					# file as such. 
					hist_date_2 = '200512'

				# maxDate from r1i1p1_raw_downloaded_files for this type
				if int(minDate) > int(198301) :
					print("***************************************")
					print("ERROR- model file data do not have early enough start date")
					print("WRONG DATES for " + s)
					print("The min date in the file was " + minDate)
					print("This file requires getting more historical data for this variable")
					print("***************************************")
					hist_date_1 = minDate

				else :
					# if the minDate is before the date required, those dates will be cut
					# so label the file accordingly. 
					hist_date_1 = '198301'

				# Write a new file name and write it, with descriptive dates
				f_seldate = s.replace('*','') + '_'+hist_date_1+'-'+hist_date_2+'.nc'
				f_seldate_out = os.path.join(cut_time_dir, f_seldate) 
				cdo.seldate('1983-01-01,2005-12-31', input=f_out, output=f_seldate_out, options="-b F64")

				# Check for the length of this cut historical file written above.
				# 1983 - 2005 is 23 years, times 12 months per year, the time dimension
				# of the cut files should be length 276
				nc = Dataset(f_seldate_out, 'r')
				t = nc.variables['time'][:]
				nc.close()
				if len(t) != 276 :
					print("***************************************")
					print("ERROR- merged and cut date file does not have correct # of months")
					print("The correct number of months should be 276, got %i " % len(t) )
					print(s)
					print("***************************************")
					# If the file is not corret, delete it! This will make it easier to see
					# what history files we are missing and then we can make changes to those
					# and fill in the gaps. 
					os.remove(f_seldate_out)

			# RCP 4.5 or 8.5 statements --------------------------------------------
			if scenario != "historical" :

				cut_time_dir = os.path.join(base_dir, "r1i1p1_rcp_cut/")

				if int(maxDate) < int(210012) :
					print("***************************************")
					print("ERROR- model last date not late enough")
					print("WRONG DATES for " + s)
					print("The max date in the file was " + maxDate)
					print("This file requires a manual addition of first month of an RCP file")
					print("***************************************")
					rcp_date_2 = maxDate
				else :
					rcp_date_2 = '210012'

				if int(minDate) > int(200601) :
					print("***************************************")
					print("ERROR- model file data do not have earlier enough start date")
					print("WRONG DATES for " + s)
					print("The min date in the file was " + minDate)
					print("***************************************")
					rcp_date_1 = minDate
				else :
					rcp_date_1 = "200601"

				# TODO: Consider cutting data to 2095 for creater participation among models,
				# TODO: e.g. evspsbl RCP85 CMCC-CESM. NOTE: I tried getting that file and failed. 
				f_seldate = s.replace('*','') + '_'+rcp_date_1+'-'+rcp_date_2+'.nc'
				f_seldate_out = os.path.join(cut_time_dir, f_seldate) 
				cdo.seldate('2006-01-01,2100-12-31', input=f_out, output=f_seldate_out, options="-b F64")	

				# Check for the length of this cut historical file written above.
				# 2006 - 2100 is 95 years, times 12 months per year, the time dimension
				# of the cut files should be length 1140
				nc = Dataset(f_seldate_out, 'r')
				t = nc.variables['time'][:]
				nc.close()

				if len(t) != 1140 :
					print("***************************************")
					print("ERROR- merged and cut date file does not have correct number of months")
					print("The correct number of months should be 1140, got %i " % len(t) )
					print(len(np.unique(t)))
					print(f_seldate_out)
					print("The file will be deleted. It is too dangerous to be left alive.")
					print("***************************************")
					# If the file is not corret, delete it! 
					os.remove(f_seldate_out)

		print("----------------------------------------------end")
		print(" ")		


# TODO: Pairing history with RCP files. 
# TODO: Regridding these merged files to the common grid for the ultimate final out! 
def pair_history_to_rcp(var, rcp, ensemble='r1i1p1') :
	"""
	This funtion pairs the cut historical files to RCP files. It can be used to 
	show what RCP files still need to be aquired. 

	Parameters
	----------
		var : str, the variable to find history and rcp pairs for
		rcp : str, the rcp to pair to the history file for a given var

	return
	------
		Writes a file concatenating RCP future files with their history
		files. written using cdo commands. A file in the model native
		resolution is written to r1i1p1_history_rcp_merged directory 
		and the same file with the new interpolated common grid is 
		written to r1i1p1_rcp_COMMON_GRID directory. 

	"""

		# e.g. variable_Amon_ModelName_scenario_ensembleMember_YYYYMM-YYYYMM.nc
	if (var == "mrso") or (var == "mrlsl.integrated") : 
		# Land parameter
		time_span = '_Lmon_'
	else :
		# atmosphere parameter
		time_span = '_Amon_'

	history_dir = os.path.join(base_dir, 'r1i1p1_history_cut', "")

	# This only works when history is available, everything builds on history files! 
	# get_var_unique_models(var, scenario, search_data_dir)
	history_models = get_var_unique_models(var, 'historical', history_dir)

	rcp_dir = os.path.join(base_dir, 'r1i1p1_rcp_cut', "")

	for model in history_models :

		print("start--------------------------------------------")

		# The dates of these files has been set and that was done by make_var_files() method
		history_file = var + time_span + model + "_historical_r1i1p1_198301-200512.nc"
		rcp_file = var + time_span + model + '_' + rcp + "_r1i1p1_200601-210012.nc"

		# Link file names to where they live on local machine 
		history_file_path = os.path.join(history_dir, history_file)
		rcp_file_path     = os.path.join(rcp_dir, rcp_file)

		print("Linking: " + history_file)
		print("To:      " + rcp_file)

		# Make sure the corresponding rcp file exists before trying to merge
		if os.path.exists(rcp_file_path) :

			# Hold onto RCP name, since that is what matters past 2005
			merged_name = var + time_span + model + "_" + rcp + "_r1i1p1_198301-210012.nc"
			history_rcp_out = os.path.join( base_dir, 'r1i1p1_history_rcp_merged', merged_name)
			
			# Make files to merge a string that cdo likes 
			files_to_merge = history_file_path + ' ' + rcp_file_path
			cdo.mergetime(input=files_to_merge, output=history_rcp_out, options="-b F64")

			# Test for correct number of months in the merged file
			# 1983-2100 is 118 years times 12 months = 1416
			nc = Dataset(history_rcp_out, 'r')
			t = nc.variables['time'][:]
			nc.close()
			# Remove and warn when the file is not the correct length. 
			if len(t) != 1416 :

				print("***************************************")
				print("ERROR- merged and cut date file does not have correct number of months")
				print(history_rcp_out)
				print("The correct number of months should be 1416, got %i " % len(t) )
				print("It has %i unique dates" %len(np.unique(t)))
				print("The file will be deleted. It is too dangerous to be left alive.")
				print("***************************************")
				os.remove(history_rcp_out)

			else :

				# Regrid this merged file to the common grid 
				f_out_common_grid = os.path.join(base_dir, 'r1i1p1_rcp_COMMON_GRID', merged_name)
				cdo.remapbil(common_grid_file, input=history_rcp_out, output=f_out_common_grid, options="-b F64")

		else :

			print('***************************************')
			print("The history file " + history_file + ' had no ' + rcp + ' pair.')
			print('Not available: ' + rcp_file)
			print('In ' + rcp_dir)
			print('***************************************')
			# TODO: Document where this happens! Then, once we know history is 
			# TODO: complete, we can easily figure out where the RCP data is missing. 

		print("----------------------------------------------end")
		print(" ")


#------------------------------------------------------------------------------
# For a given var, and scenario, cary out tasks 1-3 for all models 
#------------------------------------------------------------------------------
# NOTE: All lines to be uncommented when running all variables
if makeHistory :
	print("Working on making history files.")
	make_var_files('sfcWind', 'historical', raw_data_dir); print("sfcWind complete")
	make_var_files('tas', 'historical', raw_data_dir);     print("tas complete")
	#make_var_files('mrso', 'historical', raw_data_dir);    print("mrso complete")
	make_var_files('mrlsl.integrated', 'historical', raw_data_dir); print("mrlsl.integrated complete")
	#make_var_files('huss', 'historical', raw_data_dir);    print("huss complete")
	make_var_files('pr', 'historical', raw_data_dir);      print("pr complete")
	make_var_files('hurs', 'historical', raw_data_dir);    print("hurs complete")
	make_var_files('hfls', 'historical', raw_data_dir);    print("hfls complete")
	make_var_files('evspsbl', 'historical', raw_data_dir); print("evspsbl complete")
	print("----------------------------------")
	print("Made history files without error")
	print("----------------------------------")

if makeNew : 
	# Related to the newly downloaded data, when the project rebooted in concept
	# near the beginning of Nov 2018
	print("Make new RCP45")
	make_var_files('sfcWind', 'rcp45', raw_data_dir); print("sfcWind complete")
	make_var_files('tas', 'rcp45', raw_data_dir);     print("tas complete")
	#make_var_files('mrso', 'rcp45', raw_data_dir);    print("mrso complete")
	make_var_files('mrlsl.integrated', 'rcp45', raw_data_dir); print("mrlsl.integrated complete")
	#make_var_files('huss', 'rcp45', raw_data_dir);    print("huss complete")
	make_var_files('pr', 'rcp45', raw_data_dir);      print("pr complete")
	make_var_files('hurs', 'rcp45', raw_data_dir);    print("hurs complete")
	make_var_files('hfls', 'rcp45', raw_data_dir);    print("hfls complete")
	make_var_files('evspsbl', 'rcp45', raw_data_dir); print("evspsbl complete")

	print("Make new RCP85")
	make_var_files('sfcWind', 'rcp85', raw_data_dir); print("sfcWind complete")
	make_var_files('tas', 'rcp85', raw_data_dir);     print("tas complete")
	#make_var_files('mrso', 'rcp85', raw_data_dir);    print("mrso complete")
	make_var_files('mrlsl.integrated', 'rcp85', raw_data_dir); print("mrlsl.integrated complete")
	#make_var_files('huss', 'rcp85', raw_data_dir);    print("huss complete")
	make_var_files('pr', 'rcp85', raw_data_dir);      print("pr complete")
	make_var_files('hurs', 'rcp85', raw_data_dir);    print("hurs complete")
	make_var_files('hfls', 'rcp85', raw_data_dir);    print("hfls complete")
	make_var_files('evspsbl', 'rcp85', raw_data_dir); print("evspsbl complete")

	print("----------------------------------")
	print("Made RCP files without error")
	print("----------------------------------")

if makePairs : 

	#pair_history_to_rcp('tas', 'rcp85', ensemble='r1i1p1')
	for r in ['rcp85', 'rcp45'] :

		print("Working on pairing history to " + r)
		pair_history_to_rcp('tas', r);     print("completed tas")
		pair_history_to_rcp('sfcWind', r); print('completed sfcWind')
		#pair_history_to_rcp('mrso', r);    print('completed mrso')
		pair_history_to_rcp('mrlsl.integrated', r);    print('completed mrlsl.integrated')
		#pair_history_to_rcp('huss', r);    print('completed huss')
		pair_history_to_rcp('pr', r);      print('completed pr')
		pair_history_to_rcp('hurs', r);    print('completed hurs')
		pair_history_to_rcp('hfls', r);    print('completed hfls')
		pair_history_to_rcp('evspsbl', r); print('completed evspsbl')

	print("----------------------------------")
	print("Paired RCP to history without error")
	print("----------------------------------")

print("------------------------------------------------------------------")
print("Script executed without fail")
print("yay")
print("------------------------------------------------------------------")
