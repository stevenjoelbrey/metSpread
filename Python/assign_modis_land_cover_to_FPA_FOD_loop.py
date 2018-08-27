import gdal # For working with the Geotiff.
import os
import numpy as np
from shapely.geometry import Polygon, Point
import geopandas
import pandas as pd
import sys


# List of arguments to be passed from command line.
if len(sys.argv) == 1:
	year = 2014
	print("NO ARGUMENTS PASSED. PROCEEDING IN DEVELPMENT MODE For year", year)
else:
	print("Arguments passed")
	print 'Number of arguments:', len(sys.argv), 'arguments.'
	print 'Argument List:', str(sys.argv)
	year = np.int(sys.argv[1])
	#min_month = np.int(sys.argv[2])
	#max_month = np.int(sys.argv[3])

	print "year: " + str(year)
	#print "Min Month: " + str(min_month)
	#print "Max Month: " + str(max_month)


# Make that nice table into a dictionary that can be called.
LC_dict = {0:"Water",
           1:"Evergreen Needleleaf forest",
           2:"Evergreen Broadleaf forest",
           3:"Deciduous Needleleaf forest",
           4:"Deciduous Broadleaf forest",
           5:"Mixed forest",
           6:"Closed shrublands",
           7:"Open shrublands",
           8:"Woody savannas",
           9:"Savannas",
           10:"Grasslands",
           11:"Permanent wetlands",
           12:"Croplands",
           13:"Urban and built-up",
           14:"Cropland/Natural vegetation mosaic",
           15:"Snow and ice",
           16:"Barren or sparsely vegetated"
          }


def get_land_cover_grid():

	# Get Land cover as array using gdal. This is a huge file.
	dataDir = os.path.join("..", "..", "metSpreadData", "GIS")
	f = os.path.join(dataDir, "LCType.tif")
	gtif = gdal.Open( f )
	LC_data = np.array(gtif.GetRasterBand(1).ReadAsArray())

	print "Shape of the LC data:"
	print LC_data.shape

	# Get the coordinates
	width  = gtif.RasterXSize
	height = gtif.RasterYSize
	gt     = gtif.GetGeoTransform()

	print "Width: " + str(width)
	print "Height: " + str(height)
	print "gtiff.GetGeoTransform() returned attributes"
	print gt

	minx = gt[0]
	miny = gt[3] + width*gt[4] + height*gt[5]
	maxx = gt[0] + width*gt[1] + height*gt[2]
	maxy = gt[3]

	print "Printing the coordinate limits as defined here:"
	print "minx " + str(minx)
	print "maxx " + str(maxx)
	print "miny " + str(miny)
	print "maxy " + str(maxy)

	# Logically speaking, the data should span the length of the min and max coord values
	# by the number of dim size. Becuase these go to the limits of the data, they
	# represent upper left hand corners for these data.
	LC_x = np.linspace(minx, maxx, LC_data.shape[1])
	LC_y = np.linspace(maxy, miny, LC_data.shape[0])


	# Subset these data, as we only have wildfires in North America. This will
	# save memory! We only need data in northern hemisphere for sure.
	x_mask = np.where( (LC_x >= -169.) & (LC_x <= -60.) )[0]
	y_mask = np.where( (LC_y >= 15.) & (LC_y <= 71.5) )[0]

	LC_x_return = LC_x[x_mask]
	LC_y_return = LC_y[y_mask]
	LC_return   = LC_data[np.ix_(y_mask, x_mask)]

	print("---------------------------------------")
	print("Dimensions of returned subset LC data")

	print "Shape of the LC data:"
	print LC_return.shape

	print "Width: " + str(len(LC_x_return))
	print "Height: " + str(len(LC_y_return))
	print("---------------------------------------")

	return LC_return, LC_x_return, LC_y_return


LC_data, LC_x, LC_y = get_land_cover_grid()


################################################################################
# Perform one time operations that should not be done in the loop
################################################################################

# In order to subset the land cover grid around a fire's location, need to
# know the width of the grid in deg units. Use max observed in the data, as
# it is ok to make a slightly larger than needed area to check fire overlap.
LC_grid_dx = np.max(np.diff(LC_x))
LC_grid_dy = np.max(np.abs(np.diff(LC_y)))

# Unit conversion
m2_per_acre = 4046.86

# CRS projection to assign to grid boxes when made into polygons
crs_assign = {'init':'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'}

################################################################################
# Read in the FPA FOD data
################################################################################
print("Reading in the wildfire dataframe")
fire_dir = os.path.join("..","Data","Fire","FPA_FOD")
fpa_fod_file = os.path.join(fire_dir, "FPA_FOD_1992_2015_eco.csv")
FPA_FOD_all = pd.read_csv(fpa_fod_file, low_memory=False)

# Subset FPA_FOD_all by the selected year.
year_mask  = (FPA_FOD_all["FIRE_YEAR"] == year)

#month_mask = (FPA_FOD_all["START_MONTH"] >= min_month) & (FPA_FOD_all["START_MONTH"] <= max_month)
time_mask  = year_mask #& month_mask
FPA_FOD    = FPA_FOD_all.iloc[time_mask.values]
nFires     = FPA_FOD.shape[0] # Each row is a wildfire

# Add land cover types as a new column in the FPA_FOD DataFrame
# Initialize the value with -1 to know where values have not been updated.
for new_column_name in LC_dict.values():
  FPA_FOD = FPA_FOD.assign(new_column = np.zeros(nFires))
  FPA_FOD = FPA_FOD.rename(index=str, columns={"new_column": new_column_name})

# Add columns that will store information about the execution of the overlap 
# calculations
FPA_FOD = FPA_FOD.assign(dArea = np.zeros(nFires))
FPA_FOD = FPA_FOD.assign(percent_retained = np.full((nFires), -1))

# I need to know where the LC_columns are in the FPA_FOD DataFrame
LC_columns = np.where(pd.Series(FPA_FOD.columns).isin(LC_dict.values()))[0]

print("All required columns have been added to the wildfire DataFrame.")
print(FPA_FOD.columns)

################################################################################
# LOOP through each fire
################################################################################
print("Working on the large for loop")
for fire_index in range(1000):#range(nFires):

    if (fire_index%1000==0):
        print str(fire_index*1.0/nFires*100) + " percent complete"

    # FPA FOD FIRE_SIZE is in acres, lat and lon in decimal degrees.
    # NOTE: These methods assume that the lat lon provided are centriods.
    fire_size_acres = FPA_FOD["FIRE_SIZE"].iloc[fire_index]
    fire_lon = FPA_FOD["LONGITUDE"].iloc[fire_index]
    fire_lat = FPA_FOD["LATITUDE"].iloc[fire_index]

    # Need to get from acres to SI units so we can think about what this
    # size fire means in terms of a lat lon grid
    fire_area_m2 = fire_size_acres * m2_per_acre

    # Assuming the fire area is a circle, calculate radius from area
    # NOTE: Area needs to be in m**2
    fire_rad_m = np.sqrt(fire_area_m2/np.pi)

    # Subset the LC grid needed to encompas this circle. Create spatial extents
    # using local meters per degree values in the x and y directions.
    # NOTE: The coords of LC_grid units are degrees lon and lat.
    m_per_deg_lat = (111.*1000.) # ~111 km per deg and 1000 meter per km
    m_per_deg_lon = (111.*1000.) * np.cos(fire_lat * np.pi/180.) # 0 at pole, 90 deg N.

    # Fire radius needs to be expressed in terms of degrees for polygon creation
    fire_rad_deg =  fire_rad_m / m_per_deg_lat

    # Establish LC_grid extent where overlap calculations are needed
    # Add one LC_grid box width to make sure the circle is covered for very
    # small wildfires.
    minLat = fire_lat - (fire_rad_deg + LC_grid_dy)
    maxLat = fire_lat + (fire_rad_deg + LC_grid_dy)
    minLon = fire_lon - (fire_rad_deg + LC_grid_dy)
    maxLon = fire_lon + (fire_rad_deg + LC_grid_dy)

    # Mask these extends on the LC coordinates and grid
    lonIndex = np.where( ( (LC_x >= minLon) & (LC_x <= maxLon)) )[0]
    latIndex = np.where( ( (LC_y >= minLat) & (LC_y <= maxLat) ) )[0]
    LC_x_subset = LC_x[lonIndex]
    LC_y_subset = LC_y[latIndex]
    LC_subset   = LC_data[np.ix_(latIndex, lonIndex)]

    # Get a local dx dy for this subset of grid boxes
    LC_subset_dx = np.max(np.diff(LC_x_subset))
    LC_subset_dy = np.max(np.abs(np.diff(LC_y_subset)))

    # Make fire a circle shapely polygon object
    fire_poly = Point(fire_lon, fire_lat).buffer(fire_rad_deg)

    # NOTE: The definition of the polygon corners depends on whether the LC_grid_coords
    # NOTE: represent grid corners or grid centers.
    # TODO: Retain land cover information here.
    ncolumn = len(LC_x_subset)
    nrow    = len(LC_y_subset)

    # NOTE: This will break near edges of LC_subset.
    ij_count = 0
    LC_grid_polygon_list = [None] * (ncolumn * nrow)
    LC_cover = [None] * (ncolumn * nrow)
    for i in range(ncolumn):
        for j in range(nrow):
            # Get land cover for this grid box
            LC_cover[ij_count] = LC_subset[j, i]

            # where coords are top left of pixels.
            x0 = LC_x_subset[i]
            x1 = LC_x_subset[i] + LC_subset_dx
            y0 = LC_y_subset[j]
            y1 = LC_y_subset[j] - LC_subset_dy
            LC_grid_polygon_list[ij_count]  = Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
            ij_count = ij_count + 1 # Advance the count of created polygons

    # Create GeoSeries of local grid and fire
    grid_polys = geopandas.GeoSeries(LC_grid_polygon_list)
    fire_polys = geopandas.GeoSeries(fire_poly)

    # Setup geographic coordinate system for fire and grid boxes to live on.
    grid_polys.crs = crs_assign
    fire_polys.crs = crs_assign

    # Create GeoDataFrame, which have desirable properties and methods
    grid_polys_df = geopandas.GeoDataFrame({'geometry': grid_polys, 'grid_LC_val':LC_cover})
    fire_polys_df = geopandas.GeoDataFrame({'geometry': fire_polys, 'Fire':"FireName"})

    # Returns a DataFrame with only the geometries that are contained by both
    # GeoDataFrames. This is an expensive operation.
    res_intersection = geopandas.overlay(grid_polys_df, fire_polys_df, how='intersection')

    # Add the overlapping area to this intersection dataframe
    res_intersection['grid_overlap_fraction'] = pd.Series(res_intersection.area/fire_poly.area,
                                                          index=res_intersection.index)

    # Append Percent Overlap information to the intersection DataFrame
    res_intersection['percent_of_fire_area'] = pd.Series(res_intersection.area/fire_poly.area*100.,
                                                         index=res_intersection.index)
    # nRows is equal to the total number of land cover grid boxes that the fire circle area
    # overlaps. 
    nRows = res_intersection.shape[0]
    
    # Attributes to append to FPA FOD dataframe
    FPA_FOD["percent_retained"].iloc[fire_index] = np.sum(res_intersection['percent_of_fire_area'])
    for res_row in range(nRows) :
       
        # Get the overlap area for this overlapping grid box in res_intersection
        area_to_add = res_intersection.at[res_row, "grid_overlap_fraction"] * fire_size_acres
        
        # Translate land cover numeric to name for column indexing. i.e. where does that
        # overlapping area go? 
        LC_column_name = LC_dict[res_intersection.at[res_row, "grid_LC_val"]]
        
        # Because grid boxes in res_intersection can be of the same LC_type, need to add to existing. 
        existing_area = FPA_FOD[LC_column_name].iloc[fire_index]
        FPA_FOD[LC_column_name].iloc[fire_index] = (existing_area + area_to_add)
        
    # Sanity check, make sure the sum of the assigned land cover areas is very close to equal the area
    # of the fire. 
    dArea = np.sum(FPA_FOD.iloc[fire_index, LC_columns]) - fire_size_acres
    FPA_FOD["dArea"].iloc[fire_index] = dArea
    
print("For loop complete. Writing land cover appended dataframe to disk.")
    
# When the loop is complete, perform sanity checks on the data and append warnings
# to the file name. 
if np.max(np.abs(FPA_FOD["dArea"])) > -1e6:
    warnString = "_double_check_max_dArea"
elif (np.max(FPA_FOD["percent_of_fire_area_overlaped"]) < 99.):
    warnString = "_check_fire_overlap_percentages"
else:
    warnString = ""

# Done with looping through fires, write the large dataframe
date_attributes = str(year) + "_months=" #+ str(min_month) + "-" + str(max_month)
file_name = "FPA_FOD_with_MODIS_Land_Cover" + warnString + "_" + str(year) +".csv"
save_name = os.path.join("..", "Data","Fire", "FPA_FOD", file_name)
FPA_FOD.to_csv(save_name)