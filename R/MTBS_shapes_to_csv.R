# MTBS_shapes_to_csv.R

#------------------------------------------------------------------------------
# DESCRIPTION:
#------------------------------------------------------------------------------

# 1) Reads in MTBS data as points and assigns them to an ecoregion. CONUS ONLY!
# 2) Make monthly burn area totals by ecoregion and write them to csv files. 
# 3) Make modular so I can change my mind about what I want. 

# NOTE: There have been restrictions placed on latitude. This affects the
# NOTE: 'Marine Regime Mountain' Division. 

maxLat <- 51. # Wiggle room for fires with coords near canadian border. 
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


# Read in the region shape data data 
mtbs <- sf::st_read(dsn="Data/Fire/MTBS/", layer="mtbs_fod_pts_DD")
lat_mask <- mtbs$Lat <= maxLat 
lon_mask <- mtbs$Long > minLon
m <- lat_mask & lon_mask
mtbs <- mtbs[m,]
rm(lat_mask, lon_mask, m)



nFires <- dim(mtbs)[1]

mtbs_features <- names(mtbs)
# mtbs_desired_features <- c("Fire_ID", "Fire_Name", "Asmnt_Type", "Pre_ID",     
#                            "Post_ID", "Fire_Type", "ND_T", "IG_T",        
#                            "Low_T", "Mod_T", "High_T", "Ig_Date",    
#                            "Lat", "Long" ,"Acres", "geometry")

# Get ignition dates as Date object 
ignition_dates <- mtbs$Ig_Date
mtbs$Ig_month  <- lubridate::month(ignition_dates)
mtbs$Ig_year   <- lubridate::year(ignition_dates)

# Now a date as a POSIXct object
mtbs_fire_month <- as.POSIXct(paste(mtbs$Ig_year, mtbs$Ig_month, "15", sep="-"), tz="UTC")
mtbs$mtbs_fire_month <- mtbs_fire_month


# Read in baily ecoregions
bailys <- sf::st_read(dsn="Data/LandCover/eco-us-shp/", layer="eco_us")
bailys <- st_transform(bailys, crs=st_crs(mtbs))

pdf(file="Data/LandCover/eco-us-shp/PROVINCE.pdf", height=6, width = 10)
plot(bailys["PROVINCE"], xlim=c(-130, -60), ylim=c(23,51))
dev.off()

pdf(file="Data/LandCover/eco-us-shp/DIVISION.pdf", height=6, width = 10)
plot(bailys["DIVISION"], xlim=c(-130, -60), ylim=c(23,51))
dev.off()

# TODO: Optional spatial attribute
# # Level two ecoregion
# na_cec_eco_l2 <- sf::st_read(dsn="Data/LandCover/na_cec_eco_l2/", layer="NA_CEC_Eco_Level2")
# na_cec_eco_l2 <- st_transform(na_cec_eco_l2, crs=st_crs(mtbs))

# Perform overlap calculations with sf, has limitations given planer assumption
within_return <- sf::st_within(mtbs, bailys)

# Set up attributes to store
BAILY_DOMAIN   <- rep("", nFires)
BAILY_DIVISION <- rep("", nFires)
BAILY_PROVINCE <- rep("", nFires)
BAILY_SECTION  <- rep("", nFires)

# Place each fire where it belongs
within_row <- rep(NA, nFires)
for (i in 1:nFires){
  
  if (length(within_return[[i]]) > 0){
    
    bailys_row <- within_return[[i]]
    
    # Make assignments 
    BAILY_DOMAIN[i]   <- as.character(bailys$DOMAIN[bailys_row])
    BAILY_DIVISION[i] <- as.character(bailys$DIVISION[bailys_row])
    BAILY_PROVINCE[i] <- as.character(bailys$PROVINCE[bailys_row])
    BAILY_SECTION[i]  <- as.character(bailys$SECTION[bailys_row])
    
  }
  
}

# Place the BAILY_ATTRIBUTE onto mtbs
mtbs$BAILY_DOMAIN   <- BAILY_DOMAIN
mtbs$BAILY_DIVISION <- BAILY_DIVISION
mtbs$BAILY_PROVINCE <- BAILY_PROVINCE
mtbs$BAILY_SECTION  <- BAILY_SECTION
  
# # Where are the NA (fires without assignment) located?
# na_mask <- is.na(within_row)
# 
# quartz()
# plot(mtbs[na_mask,]$geometry, col="red", pch=19)
# title(paste(sum(na_mask), "Fires with no assignment"))
# map("world", add=T)
# map("state", add=T)

write.csv(mtbs, file="Data/Fire/MTBS/mtbs_with_bailys.csv")
save(mtbs, file="Data/Fire/MTBS/mtbs_with_bailys.RData")

# Now, we need monthly burn area for each DIVISION (or other spatial subset)
# of interest, such that the data are organized in a table
# Year/month | DIVISION 1 | DIVISION 2 | DIVISION 3 | and so on. 
# Jan 1984   |     10     |    23      |    4       |   ...

# Do not hold onto the geometry attribute
mtbs_no_geo <- as.data.frame(mtbs)

# Possible ways to replace function below with dplyr functionality
# mtcars_summary <- mtcars %>% 
#   group_by(cyl) %>%   # multiple group columns
#   summarise(max_hp = max(hp))#, mean_mpg = mean(mpg))  # multiple summary columns
# 
# mtbs_summary <- mtbs_no_geo %>% 
#   group_by(mtbs_fire_month, BAILY_DIVISION) %>% # multiple group columns
#   summarise(total_area = sum(Acres))            # multiple summary columns

################################################################################
# Now, we need monthly burn area for each DIVISION (or other spatial subset)
# of interest, such that the data are organized in a table
# Year/month | DIVISION 1 | DIVISION 2 | DIVISION 3 | and so on. 
# Jan 1984   |     10     |    23      |    4       |   ...
################################################################################

sum_acres_by_region <- function(mtbs_no_geo, region="BAILY_DIVISION"){

  # Set columns equal to the unique categories in region argument
  column_names <- unique(mtbs_no_geo[[region]])
  
  # I want time array to include ALL months between 1984 and 2016. Max span of 
  # data, zeros where needed.
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
      
      # Mask this date in the data
      year_mask <- mtbs_no_geo$Ig_year == YEAR
      month_mask <- mtbs_no_geo$Ig_month == MONTH
      
      # Overall Mask
      m <- div_mask & year_mask & month_mask
      
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
write.csv(Division_burn_area, file="Data/Fire/MTBS/bailys_division_acres_burned.csv")
write.csv(Province_burn_area, file="Data/Fire/MTBS/bailys_province_acres_burned.csv")

print("Script executed without error.")
