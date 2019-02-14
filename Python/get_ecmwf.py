#!/usr/bin/env python
import sys
from ecmwfapi import ECMWFDataServer
import os

# Read command line arguments
args  = sys.argv
if len(args) > 1 :
	var   = str(args[1])
	year1 = int(args[2])
	year2 = int(args[3])

else : 
	# Set the arguments manually
	print("Using hardcoded arguments, none passed at command line")
	var   = "tp"
	year1 = 1983
	year2 = 2017

#------------------------------------------------------------------------------
# DESCRIPTION:
#------------------------------------------------------------------------------
# This script was written to download data from ERA-interim, for various 
# forecast (fc) and analysis fields, for monthly values of surface parameters. 
# The dictionary feeding the server.retrieve() method in the for loop below
# format was determined by copying mars requests made by the apps.ecmwf.int.
# Resources describing thed data downloaded : 
# http://apps.ecmwf.int/datasets/data/interim-mdfa/levtype=sfc/
# http://apps.ecmwf.int/datasets/data/interim-full-moda/levtype=sfc/
# https://www.ecmwf.int/en/faq/what-are-definitions-radiation-fields
# 
# Detials on analysis fields and forecast fields can be seen in thier respective
# sections comments below. 

# unless otherwise specified, this is where these raw datafiles should go
saveDir = "/Users/sbrey/GoogleDrive/sharedProjects/metSpreadData/ERA-Interim"

# If that prefered saveDir does not exist, that means this work is being done
# on a different machine. Save to backup location.

if not os.path.isdir(saveDir) :
	# Save where the code is running or otherwise specified
	# SJB 01/14/2019
	#saveDir = os.getcwd()
	saveDir = '/Users/sbrey/Desktop/ecmwf'

print("Writing files to:")
print(saveDir)

# Instance of the server class. 
server = ECMWFDataServer()

# Var must be in paramDict below. 
# Dictionary created based on mars request samples created at the following URL:
# http://apps.ecmwf.int/datasets/data/interim-full-moda/levtype=sfc/?month_years=1990&param=167.128
years=range(year1, (year2+1) ) # span of year1 through year2
dateBase="YYYY0101/YYYY0201/YYYY0301/YYYY0401/YYYY0501/YYYY0601/YYYY0701/YYYY0801/YYYY0901/YYYY1001/YYYY1101/YYYY1201"



print("----------------------------------------------------------")
print("Getting variable: " + var)
print("For years: " + str(years))
print("----------------------------------------------------------")

# The param dictionary includes both analysis and forecast fields, and that
# is indicated by the second value in the list associated with each parameter. 
paramDict = {"t2m":["167.128", "an"],  # 2-meter temperature 
			 "si10":["207.128", "an"], # 10-meter wind speed
			 "d2m":["168.128", "an"],  # 2 meter dew point (needed to get hurs and huss)
			 "sp":["134.128", "an"],   # Surface Pressure
			 "e":["182.128","fc"],     # Evaporation
			 "tp":["228.128", "fc"],   # Total precipitation
			 "uvb":["57.128", "fc"],   # Surface UV radiation (not needed)
			 "ssrd":["169.128", "fc"], # Surface solar radiation downwards
			 "strd":["175.128", "fc"], # Surface thermal radiation downwards
			 "slhf":["147.128", "fc"], # surface latent heat flux
			 "ro":["205.128", "fc"],          # Runoff 
			 "soil":["39.128/40.128/41.128/42.128", "an"] # soil moisture, levels 1-4, to be combined. 
			 }

# Loop through the desired years to download. 
for year in years:
	
	# Replace "YYYY" with the year that is being downloaded 
	dates = dateBase.replace("YYYY", str(year))

	file_name = var + "_" + str(year) + ".nc"
	out = os.path.join(saveDir , file_name)

	# Forecast fields 
	# http://apps.ecmwf.int/datasets/data/interim-mdfa/levtype=sfc/
	#
	# Notes on how to get monthly totals from monthly forecast fields
	# https://confluence.ecmwf.int/pages/viewpage.action?pageId=65218804
	# 
	# "Step 0-12: averages of the precipitation produced from the first 
	# 12 hours of each 00:00 forecast, plus the precipitation produced from 
	# the first 12 hours of each 12:00 forecast. The resulting data is temporally 
	# continuous."
	
	if paramDict[var][1] == "fc" :
		
		server.retrieve({
		    "class": "ei",
		    "dataset": "interim",
		    "date": dates,
		    "expver": "1",
		    "grid": "0.75/0.75",
		    "levtype": "sfc",
		    "param": paramDict[var][0],
		    "step": "0-12",   # See note above. 
		    "stream": "mdfa", #(Monthly Means of Daily Forecast Accumulations)
		    "type": "fc", # forecast field
		    "target": out,
		    "format":"netcdf"
		})

	# Anlysis fields
	# http://apps.ecmwf.int/datasets/data/interim-full-moda/levtype=sfc/
	# https://apps.ecmwf.int/datasets/data/interim-full-moda/levtype=sfc/ for soil moisture also. ERA-Land only goes to 2010.
	elif paramDict[var][1] == "an" :
		server.retrieve({
		    "class": "ei",
		    "dataset": "interim",
		    "date": dates,
		    "expver": "1",
		    "grid": "0.75/0.75",
		    "levtype": "sfc",
		    "param": paramDict[var][0],
		    "stream": "moda",
		    "type": "an", # analysis field
		    "target": out,
		    "format":"netcdf"
		})



print("----------------------------------------------------------")
print("Script executed without fail")
print("----------------------------------------------------------")