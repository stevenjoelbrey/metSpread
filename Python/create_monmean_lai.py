import sys
import os
import cdo as cdo
import numpy as np


y1=1985
y2=1989

cdo = cdo.Cdo()

os.chdir("/Users/sbrey/GoogleDrive/sharedProjects/metSpreadData/AVHRR")


for y in np.arange(y1, (y2+1) ) : 
	
	print("Working of year %i" %y)
	f_in = 'daily_merged_files/daily_merged_' + str(y) + ".nc"
	f_out = 'monthly_mean_files/monmeans_' + str(y) + ".nc"
	print("f_in : %s" %f_in)
	print("f_out: %s" %f_out)
	print("")
	cdo.monmean(input=f_in, output=f_out)