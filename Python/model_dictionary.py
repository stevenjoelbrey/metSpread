from netCDF4 import Dataset

def monthly_file_name(var, model, rcp):
	"""Function for creating file connections for different variables
	scenarios and models. Preasently, ensemble member r1i1p1 is assumed, as
	well as 200601-210012 dime span for the monthly data."""
	
	# e.g. hfls_Amon_CNRM-CM5_rcp85_r1i1p1_200601-210012.nc
	f = var + "_" + "Amon_" + model + "_rcp" + rcp + "_r1i1p1_200601-210012.nc"
	return f

def get_nc_vals(f, var="tas"): #minLon, maxLon, minLat, maxLat):
	"""This function will be for loading a particular nc file, which will be spatially subset"""
	nc = Dataset(f)
	vals = nc.variables[var][:]
	t = nc.variables["time"][:]
	lon = nc.variables["lon"][:]
	lat = nc.variables["lat"][:]

	return vals, t, lon, lat


model_name = {
		"ACCESS1-0":"ACCESS1-0",
		"ACCESS1-3":"ACCESS1-3",
		"CanESM2":"CanESM2",
		"CCSM4":"CCSM4",
		"CMCC-CESM":"CMCC-CESM",
		"CMCC-CM":"CMCC-CM",
		"CMCC-CMS":"CMCC-CMS",
		"CNRM-CM5":"CNRM-CM5",
		"FGOALS-g2":"FGOALS-g2",
		"FGOALS-s2":"FGOALS-s2",
		"CSIRO-Mk3-6-0":"CSIRO-Mk3-6-0",
		"GFDL-CM3":"GFDL-CM3",
		"GFDL-ESM2G":"GFDL-ESM2G",
		"GFDL-ESM2M":"GFDL-ESM2M",
		"GISS-E2-H":"GISS-E2-H",
		"GISS-E2-H-CC":"GISS-E2-H-CC",
		"GISS-E2-R":"GISS-E2-R",
		"GISS-E2-R-CC":"GISS-E2-R-CC",
		"HadGEM2-CC":"HadGEM2-CC",
		"HadGEM2-ES":"HadGEM2-ES",
		"inmcm4":"inmcm4",
		"IPSL-CM5A-LR":"IPSL-CM5A-LR",
		"IPSL-CM5A-MR":"IPSL-CM5A-MR",
		"IPSL-CM5B-LR":"IPSL-CM5B-LR",
		"MIROC-ESM":"MIROC-ESM",
		"MIROC-ESM-CHEM":"MIROC-ESM-CHEM",
		"MIROC4h":"MIROC4h",
		"MIROC5":"MIROC5",
		"MPI-ESM-LR":"MPI-ESM-LR",
		"MPI-ESM-MR":"MPI-ESM-MR",
		"MPI-ESM-P":"MPI-ESM-P", 
		"MRI-CGCM3":"MRI-CGCM3",
		"MRI-ESM1":"MRI-ESM1",
		"NorESM1-M":"NorESM1-M",
		"NorESM1-ME":"NorESM1-ME"
}

# A dictionary that uses the model name and shows the modeling center
# Full names of institutions and model names available here:
# https://cmip.llnl.gov/cmip5/availability.html
model_institution = {
		"ACCESS1-0":"CSIRO-BOM",
		"ACCESS1-3":"CSIRO-BOM",
		"CanESM2":"CCCma",
		"CCSM4":"NCAR",
		"CMCC-CESM":"CMCC",
		"CMCC-CM":"CMCC",
		"CMCC-CMS":"CMCC",
		"CNRM-CM5":"CNRM-CERFACS",
		"CSIRO-Mk3-6-0":"CSIRO-QCCCE",
		"GFDL-CM3":"NOAA-GFDL",
		"GFDL-ESM2G":"NOAA-GFDL",
		"GFDL-ESM2M":"NOAA-GFDL",
		"GISS-E2-H":"NASA GISS",
		"GISS-E2-H-CC":"NASA GISS",
		"GISS-E2-R":"NASA GISS",
		"GISS-E2-R-CC":"NASA GISS",
		"HadGEM2-CC":"MOHC",
		"HadGEM2-ES":"MOHC",
		"inmcm4":"RINM",
		"IPSL-CM5A-LR":"IPSL",
		"IPSL-CM5A-MR":"IPSL",
		"IPSL-CM5B-LR":"IPSL",
		"MIROC-ESM":"MIROC",
		"MIROC-ESM-CHEM":"MIROC",
		"MIROC5":"MIROC",
		"MPI-ESM-LR":"MPI-M",
		"MPI-ESM-MR":"MPI-M",
		"MRI-CGCM3":"MRI",
		"MRI-ESM1":"MRI",
		"NorESM1-M":"NCC",
		"NorESM1-ME":"NCC"
}

