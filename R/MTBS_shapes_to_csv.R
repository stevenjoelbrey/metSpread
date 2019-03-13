# MTBS_shapes_to_csv.R

#------------------------------------------------------------------------------
# DESCRIPTION:
#------------------------------------------------------------------------------

# 1) Reads in MTBS data as points and assigns them to an ecoregion. CONUS ONLY!
# 2) Make monthly burn area totals by ecoregion and write them to csv files. 
# 3) Make modular so I can change my mind about what I want. 

# NOTE: There have been restrictions placed on latitude. This affects the
# NOTE: 'Marine Regime Mountain' Division. This is by design. 

# NOTE: Because of our interest in biophysical controls on fire area we exclude
# NOTE: RX Types from these analysis. Prescribed fires are labeled both "RX" and 
# NOTE: "Rx" in these data, so we make sure to catch both types. 

maxLat <- 51. # Wiggle room for fires with coords near Canadian border. 
minLon <- -134. # Wiggle room for wildfires near the coast in Washington. 
# Ecoregions used for classification:
# https://www.fs.fed.us/rm/ecoregions/products/map-ecoregions-united-states/#
# https://www.fs.fed.us/rm/ecoregions/downloads/ecoregions-united-states/data-sets/eco-us-shp.zip
# 
# Fire Data URL: https://www.mtbs.gov/direct-download
# Ecoregion data URL: 

library(sf)
library(lubridate)
library(sp)
library(dplyr)
library(maps)

# Read in the region shape data data 
mtbs <- sf::st_read(dsn="Data/Fire/MTBS/", layer="mtbs_fod_pts_DD")

print("Unique 'Fire_Type' labels")
print(unique(mtbs$Fire_Type))

# TRUE where we want to keep
lat_mask <- mtbs$Lat <= maxLat 
lon_mask <- mtbs$Long > minLon
type_mask <- (mtbs$Fire_Type != "RX") & (mtbs$Fire_Type != "Rx")
# Combination of the masks
m <- lat_mask & lon_mask & type_mask 

# Subset the data 
mtbs <- mtbs[m,]
rm(lat_mask, lon_mask, m)

# Sanity check
print("Unique 'Fire_Type' labels after masking")
print(unique(mtbs$Fire_Type))

# Get total number of wildfires in area of interest
nFires <- dim(mtbs)[1]

mtbs_features <- names(mtbs)
print(mtbs_features)
# mtbs_desired_features <- c("Fire_ID", "Fire_Name", "Asmnt_Type", "Pre_ID",     
#                            "Post_ID", "Fire_Type", "ND_T", "IG_T",        
#                            "Low_T", "Mod_T", "High_T", "Ig_Date",    
#                            "Lat", "Long" ,"Acres", "geometry")

# Sanity checks 
quartz()
plot(mtbs$Long, mtbs$Lat, pch=".", col="red")
map("state", add=T)
title("Alaska and HI excluded, fire locations to grid")

# Get ignition dates as Date object 
ignition_dates <- mtbs$Ig_Date
# Add ignition month and year columns 
mtbs$Ig_month  <- lubridate::month(ignition_dates)
mtbs$Ig_year   <- lubridate::year(ignition_dates)

# Now make a date as a POSIXct object (why?)
mtbs_fire_month <- as.POSIXct(paste(mtbs$Ig_year, mtbs$Ig_month, "15", sep="-"), tz="UTC")
mtbs$mtbs_fire_month <- mtbs_fire_month

# Read in baileys ecoregions. These are the shapes by which we separate and 
# count the wildfire burn area
bailys <- sf::st_read(dsn="Data/LandCover/eco-us-shp/", layer="eco_us")
bailys <- st_transform(bailys, crs=st_crs(mtbs))

pdf(file="Data/LandCover/eco-us-shp/PROVINCE.pdf", height=6, width = 10)
plot(bailys["PROVINCE"], xlim=c(-130, -60), ylim=c(23,51))
#map("state", add=T)
#points(mtbs$Long, mtbs$Lat, pch=".")
dev.off()

pdf(file="Data/LandCover/eco-us-shp/DIVISION.pdf", height=6, width = 10)
plot(bailys["DIVISION"], xlim=c(-130, -60), ylim=c(23,51))
#map("state", add=T)
#points(mtbs$Long, mtbs$Lat, pch=".")
dev.off()

pdf(file="Data/LandCover/eco-us-shp/baily_polygons_and_fires.pdf", height=6, width = 10)
plot(bailys$geometry, xlim=c(-130, -60), ylim=c(23,51))
map("state", add=T)
points(mtbs$Long, mtbs$Lat, pch=".", col="red")
title("Baily polyons (smallest level) + states + wildfire centriods(?)")
dev.off()

# TODO: Optional spatial attribute
# # Level two ecoregion
# na_cec_eco_l2 <- sf::st_read(dsn="Data/LandCover/na_cec_eco_l2/", layer="NA_CEC_Eco_Level2")
# na_cec_eco_l2 <- st_transform(na_cec_eco_l2, crs=st_crs(mtbs))

# Perform overlap calculations with sf, has limitations given planer assumption
# returns the row of dim(bailys) where mtbs point falls inside of in bailys
within_return <- sf::st_within(mtbs, bailys)

# Set up attributes to store, length of number of fires, because we want to know
# answer for every fire. 
BAILY_DOMAIN   <- rep("", nFires)
BAILY_DIVISION <- rep("", nFires) # Of primary interest!!! 
BAILY_PROVINCE <- rep("", nFires)
BAILY_SECTION  <- rep("", nFires)

# Place each fire where it belongs
within_row <- rep(NA, nFires)

# Keeps track of the fires that overlapped a polygon vs. those that did not
had_overlap <- rep(TRUE, nFires)

for (i in 1:nFires){
  
  if (length(within_return[[i]]) > 0){
    
    bailys_row <- within_return[[i]]
    
    # Make assignments 
    BAILY_DOMAIN[i]   <- as.character(bailys$DOMAIN[bailys_row])
    BAILY_DIVISION[i] <- as.character(bailys$DIVISION[bailys_row])
    BAILY_PROVINCE[i] <- as.character(bailys$PROVINCE[bailys_row])
    BAILY_SECTION[i]  <- as.character(bailys$SECTION[bailys_row])
    
  } else {
    # No assignment made, keep track of these, we definitly want to plot thier
    # locations and track down just how much area we lost in this process. 
    had_overlap[i] <- FALSE
  }
  
}

# Place the BAILY_ATTRIBUTE onto mtbs sf:dataframe as new columns 
mtbs$BAILY_DOMAIN   <- BAILY_DOMAIN
mtbs$BAILY_DIVISION <- BAILY_DIVISION
mtbs$BAILY_PROVINCE <- BAILY_PROVINCE
mtbs$BAILY_SECTION  <- BAILY_SECTION

########################################################
# Where are the NA (fires without assignment) located?
########################################################

quartz()
plot(mtbs[!had_overlap,]$geometry, col="red", pch=19)
title(paste(sum(!had_overlap), "Fires with no assignment totaling :" ,
            sum(mtbs$Acres[!had_overlap]) ,"acres"))
map("world", add=T)
map("state", add=T)

pdf(file="Data/LandCover/eco-us-shp/fires_without_within.pdf", height=6, width = 10)
plot(mtbs[!had_overlap,]$geometry, col="red", pch=19)
title(paste(sum(!had_overlap), "Fires with no assignment totaling :" ,
            sum(mtbs$Acres[!had_overlap]) ,"acres"))
map("world", add=T)
map("state", add=T)
dev.off()

################################################################################
# Plot the three DIVISIONS of interest to ensure separation of assignment
################################################################################
divisions <- c("Marine Regime Mountains", "Temperate Steppe Regime Mountains",
               "Mediterranean Regime Mountains")

pdf(file="Data/LandCover/eco-us-shp/mountain_fires.pdf", height=6, width = 10)

map("state")
m1 <- mtbs$BAILY_DIVISION == "Marine Regime Mountains"
plot(mtbs[m1,]$geometry, col="blue", pch=".", add=T)

m2 <- mtbs$BAILY_DIVISION == "Temperate Steppe Regime Mountains"
plot(mtbs[m2,]$geometry, col="red", pch=".", add=T)

m3 <- mtbs$BAILY_DIVISION == "Mediterranean Regime Mountains"
plot(mtbs[m3,]$geometry, col="green", pch=".", add=T)

dev.off()


################################################################################
# Write the mtbs data as a csv. 
################################################################################
write.csv(mtbs, file="Data/Fire/MTBS/mtbs_with_bailys.csv")
save(mtbs, file="Data/Fire/MTBS/mtbs_with_bailys.RData")

################################################################################
# Now, we need monthly burn area for each DIVISION (or other spatial subset)
# of interest, such that the data are organized in a table
# Year/month | DIVISION 1 | DIVISION 2 | DIVISION 3 | and so on. 
# Jan 1984   |     10     |    23      |    4       |   ...

# Do not hold onto the geometry attribute
mtbs_no_geo <- as.data.frame(mtbs)

# TODO: Possible ways to replace function below with dplyr functionality
# mtcars_summary <- mtcars %>% 
#   group_by(cyl) %>%   # multiple group columns
#   summarise(max_hp = max(hp))#, mean_mpg = mean(mpg))  # multiple summary columns
# 
# mtbs_summary <- mtbs_no_geo %>% 
#   group_by(mtbs_fire_month, BAILY_DIVISION) %>% # multiple group columns
#   summarise(total_area = sum(Acres))            # multiple summary columns

sum_acres_by_region <- function(mtbs_no_geo, region="BAILY_DIVISION"){

  print(paste("Using", region, "as geographic separator"))
  
  # Set columns equal to the unique categories in region argument
  column_names <- unique(mtbs_no_geo[[region]])
  
  # I want time array to include ALL months between 1984 and 2016. Max span of 
  # data, zeros where needed, or rather, where there is no burn area
  d1 <- as.POSIXct("1984-01-01", tz="UTC")
  d2 <- as.POSIXct("2016-12-01", tz="UTC")
  month_sequence <- seq(from=d1, to=d2, by="month")
  row_names    <- as.character(month_sequence)
  
  # Create an empty dataframe to store these data
  df <- data.frame(matrix(NA, 
                          nrow = length(row_names), 
                          ncol = length(column_names))
                   )
  
  # Give the data nice column and row names
  colnames(df) <- column_names
  row.names(df) <- row_names
  
  # Loop through unique regions
  for (j in 1:length(column_names)) {
    
    div_mask <- mtbs_no_geo[[region]] == column_names[j]
    
    # sum the monthly burn area total
    for (i in 1:length(month_sequence)){
      
      MONTH <- lubridate::month(month_sequence[i])
      YEAR  <- lubridate::year(month_sequence[i])
      
      # Mask this date in the data (all regions in data still)
      year_mask <- mtbs_no_geo$Ig_year == YEAR
      month_mask <- mtbs_no_geo$Ig_month == MONTH
      
      # Overall Mask, combines year, month, and div
      m <- (div_mask) & (year_mask) & (month_mask)
      
      # If there is a TRUE value in M, there are acres that need to be summed
      if(sum(m)>0){
        df[i, j] <- sum(mtbs_no_geo$Acres[m])  
      } else{
        df[i, j] <- 0
      }
      
    } # End of month loop
    
  } # end of region loop

  return(df)
  
}

Division_burn_area <- sum_acres_by_region(mtbs_no_geo, region="BAILY_DIVISION")
Province_burn_area <- sum_acres_by_region(mtbs_no_geo, region="BAILY_PROVINCE")

# Write these region summaries 
# THESE ARE THE ONES USED FOR ANALYSIS
write.csv(Division_burn_area, file="Data/Fire/MTBS/bailys_division_acres_burned.csv")
write.csv(Province_burn_area, file="Data/Fire/MTBS/bailys_province_acres_burned.csv")

print("Script executed without error.")
