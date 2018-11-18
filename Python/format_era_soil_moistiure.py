#!/usr/bin/env python

#------------------------------------------------------------------------------
# Description
#------------------------------------------------------------------------------
# The purpose of this script is to combine era-interim soil layers 1-4 
# into a single total volume of soil moisture. layer 1 from 0cm to 7cm, layer 
# 2 from 7cm to 28cm, layer 3 from 28cm to 1m and layer 4 from 1m to 2.89m. 
# The yearly soil files were downloaded with get_ecmwf.py. They are labeled 
# soil_YYYY.nc. They contain four soil layers, swvl1, swvl2, swvl3, swvl4. We 
# want to combine these layers because CMIP5 tracks soil moisture in all layers (mrso).
# ERA soil units are m**3/m**3, this needs to be converted to kg/m**2 in order
# to match th units of CMIP5 soil moisture. The sketch below is to aid this 
# unit conversion.  
#
#    *-------*
#   /|      /|
#  / |     / |
# *--|----*  |
# |  *----|--*
# | /     | /
# *-------*
#
# ERA gives us the fraction of a volume of soil of depth Z that is water in some
# phase. The depth of this volume is 2.89 m, and the unit area of interest is 1 m**2. 
# The total volume of this volume is 2.89 m * 1 m**2 = 2.89 m**3. ERA tells us what 
# fraction is water. total_water_vol = 2.89 m**3 * frac. To get the mass in this cube
# multiply by the density of liquid water, which will have some error, since some of
# this volume could be ice. 
# total_water_mass_per_area = water_vol_frac * 1000 kg/m**3 * 2.89 m = kg/m**2 


import sys
import os
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


dataDir = "/Users/sbrey/GoogleDrive/sharedProjects/metSpreadData/ERA-Interim"

cdo = cdo.Cdo()

years = np.arange(year1, year2+1)

for y in years : 
	f_in = os.path.join(dataDir, 'soil_' + str(y) + '.nc')
	f_out = os.path.join(dataDir, 'soil_vol_all_layers' + str(y) + '.nc')
	# TODO: This addition needs to be updated. Not all soil moisture levels 
	# TODO: are created equal. For example, 0.5 of the total volume for a layer
	# TODO: 7 cm deep is not the same volume as 0.5 of the total volume for a layer
	# TODO: that is 1 m deep. Grrrr. 
	cdo.expr('soil_vol_all_layers=swvl1+swvl2+swvl3+swvl4', input=f_in, output=f_out)

# Those created files need to be converted from units of m**3/m**3 to 
# kg/m**2
for y in years : 
	f_in = os.path.join(dataDir, 'soil_vol_all_layers' + str(y) + '.nc')
	f_out = os.path.join(dataDir, 'mrso' + str(y) + '.nc')
	cdo.expr('mrso=soil_vol_all_layers*2890', input=f_in, output=f_out)
	# TODO: Change units label here and then finally call it mrso 
