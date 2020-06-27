rm(list=ls()) # start with a clean R environment

my_dir <- "D:/Postdoc_Uppsala" # path of the folder where your data are stored
setwd(my_dir)

# Required libraries ------------------------------------------------------------------
library(rgbif)
library(raster)
library(CoordinateCleaner)
library(biomod2)
library(plyr)
library(usdm)
library(rgdal)
library(dismo)

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ------------------------ MATERIAL & METHODS --------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

# Previous data & information needed for analysis ----------------------------
## World raster:
# Empty raster created based on the land-use data, with the same dimensions and resolution
ref_raster <- raster("ref_raster_30k.grd")

## Europe raster:
# Subset from "ref_raster" above, created in ArcMap 10.5, by selecting the countries included in the Appendix S2 of the manuscript from a freely available shapefile containing the World land surface.
eu_raster <- raster("eu_ref_raster.tif")

## List of species to use:
list_sp <- read.csv("0Tables/Invasive_spp_to_use.csv",sep=";",header=F)
list_sp <- list_sp[c(1,4,2,5:7,3,8:15),1] # select the column with the name consistent throughout the dataset, and in alphabetical order

# Change "Mustela vison" by "Neovison vison", as it is named in new taxonomies
list_sp2 <- revalue(list_sp, c("Mustela vison"="Neovison vison"))


# DATA COMPILATION -----------------------------------------------------------

## SPECIES PRESENCES ###############

### GBIF ###########################
# Download OCURRENCES (all observations available in GBIF) ################
for (i in 1:length(list_sp)) {
  
  # Species to model:
  sp <- list_sp[i]
  sp2 <- sub(" ","_",sp)
  
  # Getting GBIF entries with this name: 
  all_occ <- name_backbone(sp)
  # Getting the species key, unique for the species and its synonimous:
  key_no <- all_occ$speciesKey
  # Get the number of occurrences in gbif:
  num_occ <- occ_count(taxonKey=key_no, georeferenced=TRUE)
  
  # If the species has less than 200,000 ocurrences, I can use R directly to   download the data:
  if (num_occ<=200000) {
    # Search for GBIF occurrences according to 'key_no'; add necesary filters:
    occ_search <- occ_search(taxonKey = key_no, return='data', limit = num_occ, hasCoordinate=TRUE, fields=c("name","key","decimalLatitude","decimalLongitude","basisOfRecord","year","collectionCode","institutionCode","coordinateUncertaintyInMeters"))
    occsub <- subset(occ_search, basisOfRecord=="HUMAN_OBSERVATION" | basisOfRecord=="MACHINE_OBSERVATION" | basisOfRecord=="OBSERVATION")
  } 
  else  { next }
  
  # Saving as .csv:
  file1 <- paste("Occ_vs_pres/",sp2,"_occ.csv",sep="")
  write.table(occsub,file1,sep=",",row.names=FALSE)
  
  # Converting to a spatial object
  data.sp <- SpatialPointsDataFrame(coords=cbind(occsub$decimalLongitude, occsub$decimalLatitude),proj4string = CRS(projection(ref_raster)),data=as.data.frame(occsub))
  
  # Saving as .shp:
  my_dsn <- "Occ_vs_pres"
  my_layer <- paste(sp2,"_occ",sep="")
  writeOGR(data.sp, dsn= my_dsn,layer=my_layer, driver="ESRI Shapefile", overwrite_layer=TRUE)
}

#**** For species with more than 200,000 (Branta canadiensis & Oxyura jamacensis), download directly from the GBIF web and store in your computer. These data should be filtered by human observation/machine observation/observation directly in the website. In this particular case, data were downloaded on 13th November 2018. Note that new data should be downloaded to be coherent with data downloaded in the previous step (direct download from R): ----
oxy_data <- fread("Occ_vs_pres/Oxyura_jamaicensis_occ.csv", select=c("name","key","decimalLatitude","decimalLongitude","basisOfRecord","year","collectionCode","institutionCode","coordinateUncertaintyInMeters"),quote="")
# Just in case, remove those entries without coordinates
oxy_data <- oxy_data[!is.na(oxy_data$decimalLongitude),]

branta_data <- fread("Occ_vs_pres/Branta_canadensis_occ.csv", select=c("name","key","decimalLatitude","decimalLongitude","basisOfRecord","year","collectionCode","institutionCode","coordinateUncertaintyInMeters"),quote="")
# Just in case, remove those entries without coordinates
branta_data <- branta_data[!is.na(branta_data$decimalLongitude),]

list1 <- list(oxy_data,branta_data)
list_sp3 <- c("Oxyura_jamaicensis", "Branta_canadensis")

for (i in 1:2) {
  
  ## Species to model:
  sp2 <- list_sp3[i]
  
  ## Data:
  occ1 <- list1[[i]]
  
  ## Converting to a spatial object
  data.sp <- SpatialPointsDataFrame(coords=cbind(occ1$decimalLongitude, occ1$decimalLatitude),proj4string = CRS(projection(ref_raster)),data=as.data.frame(occ1))
  data.sp$sp <- sp2
  
  ## Saving as shapefile:
  my_dsn <- "Occ_vs_pres"
  my_layer <- paste(sp2,"_occ",sep="")
  writeOGR(data.sp, dsn= my_dsn,layer=my_layer, driver="ESRI Shapefile", overwrite_layer=TRUE)
  
}

# Obtain PRESENCES (Filters: uncertainty of the observation & CoordinateCleaner. Then simplify to one observation per grid cell). Two options here: "cert" and "cert_na" sets of data (see Supplementary Information for a complete explanation) ################

radius_uncert <- 15000 # Grid size resolution (side), expressed in meters, divided by two, as uncertainty is expressed in radius around the occurrence data. It will be used to exclude occurrences with an uncertainty larger than this value.

for (j in 1:length(list_sp)) {
  
  # Species to model:
  sp <- list_sp[j]
  sp2 <- sub(" ","_",sp)
  
  # Getting occurrence data, downloaded directly from GBIF:
  name_dest <- paste("Occ_vs_pres/", sp2, "_occ.shp",sep="") # shape to allow reading also the Branta and Oxyurus 
  data.sp <- readOGR(name_dest) # spatial data (shapefile format)
  # data.sp <- data.sp@data     # only data.frame, excluding the spatial part
  
  ## Removing those with uncertainty > radius_uncert and NAs
  data_pres_cert <- subset(data.sp,crdnUIM<=radius_uncert)
  
  # Changing the column name where the species name is, so next function works
  names(data_pres_cert)[names(data_pres_cert) == 'sp'] <- 'name'
  
  ## Filtering with "coordinate cleaner"
  data_pres_cert2 <- clean_coordinates(data_pres_cert, lon="dcmlLng", lat="dcmlLtt",species="name",tests = c("centroids","equal","gbif", "institutions","zeros"))
  
  # Exclude problematic records, base on previous function:
  data_pres_cert3 <- data_pres_cert[data_pres_cert2$.summary,]
  
  ## Getting unique values per grid cell
  coord_data <- data_pres_cert3[,c("dcmlLng","dcmlLtt")]
  extract1 <- raster::extract(ref_raster,coord_data,cellnumbers = T)
  coord_raster_ref_sp <- as.data.frame(coordinates(ref_raster)[extract1[,1],])
  df1 <- as.data.frame(coord_raster_ref_sp) # OCURRENCES
  df2 <- df1[!duplicated(df1[,c('x','y')]),]
  df2 <- na.omit(df2) # PRESENCES
  
  # Converting to a spatial object - OCURRENCES
  occ.sp <- SpatialPointsDataFrame(coords=cbind(df1$x, df1$y),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),data=as.data.frame(cbind(df1$x, df1$y)))
  names(occ.sp) <- c("x","y")
  
  ## Writting points.shp
  my_dsn <- "Occ_vs_pres"
  my_layer4 <- paste(sp2,"_occ_cert_30k",sep="")
  writeOGR(occ.sp, dsn= my_dsn,layer=my_layer4, driver="ESRI Shapefile", overwrite_layer=TRUE)
  
  # Converting to a spatial object - PRESENCES
  pres.sp <- SpatialPointsDataFrame(coords=cbind(df2$x, df2$y),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),data=as.data.frame(cbind(df2$x, df2$y)))
  names(pres.sp) <- c("x","y")
  
  ## Writting points.shp
  my_dsn <- "Occ_vs_pres"
  my_layer2 <- paste(sp2,"_presences_cert_30k",sep="")
  writeOGR(pres.sp, dsn= my_dsn,layer=my_layer2, driver="ESRI Shapefile", overwrite_layer=TRUE)
  
  ## --------------------------------------
  ### Keeping the NAs as well as those with certain data: -------------------------------------
  data_pres_na <- data.sp[is.na(data.sp$crdnUIM),]
  
  # Changing the column name where the species name is, so next function works
  names(data_pres_na)[names(data_pres_na) == 'sp'] <- 'name'
  
  ## Filtering with "coordinate cleaner"
  data_pres_na2 <- clean_coordinates(data_pres_na, lon="dcmlLng", lat="dcmlLtt",species="name",tests = c("centroids","equal","gbif", "institutions","zeros"))
  
  # Exclude problematic records according to the previous function:
  data_pres_na3 <- data_pres_na[data_pres_na2$.summary,]
  
  # Joining to the previous set (cert) to save it:
  data_pres2 <- rbind(data_pres_cert3,data_pres_na3)
  
  ## Getting unique values per grid cell
  coord_data <- data_pres2[,c("dcmlLng","dcmlLtt")]
  extract1 <- raster::extract(ref_raster,coord_data,cellnumbers = T)
  coord_raster_ref_sp <- as.data.frame(coordinates(ref_raster)[extract1[,1],])
  df1 <- as.data.frame(coord_raster_ref_sp) # OCCURRENCES
  df2 <- df1[!duplicated(df1[,c('x','y')]),]
  df2 <- na.omit(df2) ## PRESENCES
  
  # Converting to a spatial object - OCURRENCES
  occ.sp2 <- SpatialPointsDataFrame(coords=cbind(df1$x, df1$y),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),data=as.data.frame(cbind(df1$x, df1$y)))
  names(occ.sp2) <- c("x","y")
  
  ## Writting points.shp
  my_dsn <- "Occ_vs_pres"
  my_layer4 <- paste(sp2,"_occ_cert_na_30k",sep="")
  writeOGR(occ.sp2, dsn= my_dsn,layer=my_layer4, driver="ESRI Shapefile", overwrite_layer=TRUE)
  
  # Converting to a spatial object - PRESENCES
  pres.sp2 <- SpatialPointsDataFrame(coords=cbind(df2$x, df2$y),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),data=as.data.frame(cbind(df2$x, df2$y)))
  names(pres.sp2) <- c("x","y")
  
  ## Writting points.shp
  my_dsn <- "Occ_vs_pres"
  my_layer2 <- paste(sp2,"_presences_cert_na_30k",sep="")
  writeOGR(pres.sp2, dsn= my_dsn,layer=my_layer2, driver="ESRI Shapefile", overwrite_layer=TRUE)
}


## PREDICTORS #####################

### CLIMATE #################
#### Downloaded from: <http://chelsa-climate.org/downloads/>, on 24th April 2017
#### Original resolution: 30 arc-sec (~1 km at the Equator) 
#### Original projection: WGS84
#### World

## Decide if resampling or averaging: 
 
clim_1k <- raster::stack(list.files(path = "CHELSA_clim/Present_original", pattern = "^CHELSA.*.tif$", full.names=TRUE))

## Resample/average to 0.25 degrees (~30 km at the Equator) to test which option is more accurate to up-scale data:

# Resampling. Without projection, for correlations of random samples
for (k in 1:dim(clim_1k)[3]) {
  name <-  names(clim_1k)[k]
  # Resample to the 5-min reference layer (built0). Non-projected, to avoid too many interpolations
  layer <- resample(clim_1k[[k]],ref_raster,method='bilinear') # 
  # Mask
  layer <- raster::mask(layer,ref_raster)
  # Save in the folder "resampled10k":
  file <- paste("CHELSA_clim/resampled_30k/",name,"_noproj_r.grd",sep="")
  writeRaster(layer,file,format="raster",prj=TRUE,overwrite=TRUE)
}

# Averaging. Without projection, for correlations of random samples
for (k in 1:dim(clim_1k)[3]) {
  name <- names(clim_1k)[k]
  # Averaging within the reference layer. Non-projected, to avoid too many interpolations
  layer <- aggregate(clim_1k[[k]],fact=30, fun=mean, expand=TRUE) # 
  # Allineate
  layer <- resample(layer,ref_raster,method='bilinear') # 
  # Mask
  layer <- raster::mask(layer,ref_raster)
  # Save in the folder "resampled10k":
  file <- paste("CHELSA_clim/resampled_30k/",name,"_noproj_av.grd",sep="")
  writeRaster(layer,file,format="raster",prj=TRUE,overwrite=TRUE)
}

# Reading results withouth projection (resample): [RUN ONLY IF YOU DON'T WANT TO RE-RUN THE LOOPS ABOVE AND YOU HAVE ALREADY YOUR DATA READY]
clim_10k_r <- stack(list.files(path = "CHELSA_clim/resampled_30k/", pattern = "noproj_r.grd$", full.names=TRUE))
# Reading results withouth projection (average):
clim_10k_a <- stack(list.files(path = "CHELSA_clim/resampled_30k/", pattern = "noproj_av.grd$", full.names=TRUE))

## What is better? 
# Extracting values from here and from the original data to check for correlations: 

# Random coordinates to extract from all data:
clim_10k_r_sample <- as.data.frame(sampleRandom(clim_10k_r,size=1000,xy=TRUE))
names(clim_10k_r_sample) 

coord_sample <- clim_10k_r_sample[,c(1:2)]

clim_10k_a_sample <- raster::extract(clim_10k_a,coord_sample,df=TRUE)
colnames(clim_10k_a_sample) <- paste(colnames(clim_10k_a_sample),"av",sep="_")

clim_1k_sample <- raster::extract(clim_1k,coord_sample,df=TRUE)
colnames(clim_1k_sample) <- paste(colnames(clim_1k_sample),"1k",sep="_")

table_clim1 <- cbind(clim_10k_r_sample,clim_1k_sample)
table_clim2 <- cbind(clim_10k_a_sample,clim_1k_sample)

# Correlations:
variable <- c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12", "bio13","bio14","bio15","bio16","bio17","bio18","bio19")

spearman_corr_r <- c(cor(table_clim1$CHELSA_bio10_1,table_clim1$CHELSA_bio10_1_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim1$CHELSA_bio10_2,table_clim1$CHELSA_bio10_2_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim1$CHELSA_bio10_3,table_clim1$CHELSA_bio10_3_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim1$CHELSA_bio10_4,table_clim1$CHELSA_bio10_4_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim1$CHELSA_bio10_5,table_clim1$CHELSA_bio10_5_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim1$CHELSA_bio10_6,table_clim1$CHELSA_bio10_6_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim1$CHELSA_bio10_7,table_clim1$CHELSA_bio10_7_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim1$CHELSA_bio10_8,table_clim1$CHELSA_bio10_8_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim1$CHELSA_bio10_9,table_clim1$CHELSA_bio10_9_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim1$CHELSA_bio10_10,table_clim1$CHELSA_bio10_10_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim1$CHELSA_bio10_11,table_clim1$CHELSA_bio10_11_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim1$CHELSA_bio10_12,table_clim1$CHELSA_bio10_12_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim1$CHELSA_bio10_13,table_clim1$CHELSA_bio10_13_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim1$CHELSA_bio10_14,table_clim1$CHELSA_bio10_14_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim1$CHELSA_bio10_15,table_clim1$CHELSA_bio10_15_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim1$CHELSA_bio10_16,table_clim1$CHELSA_bio10_16_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim1$CHELSA_bio10_17,table_clim1$CHELSA_bio10_17_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim1$CHELSA_bio10_18,table_clim1$CHELSA_bio10_18_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim1$CHELSA_bio10_19,table_clim1$CHELSA_bio10_19_1k,use="pairwise.complete.obs",method="spearman"))

corr_clim_r <- cbind(variable,spearman_corr_r)

write.table(corr_clim_r,"CHELSA_clim/resampled_30k/corr_clim_r.txt")

spearman_corr_a <- c(cor(table_clim2$CHELSA_bio10_1_av,table_clim2$CHELSA_bio10_1_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim2$CHELSA_bio10_2_av,table_clim2$CHELSA_bio10_2_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim2$CHELSA_bio10_3_av,table_clim2$CHELSA_bio10_3_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim2$CHELSA_bio10_4_av,table_clim2$CHELSA_bio10_4_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim2$CHELSA_bio10_5_av,table_clim2$CHELSA_bio10_5_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim2$CHELSA_bio10_6_av,table_clim2$CHELSA_bio10_6_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim2$CHELSA_bio10_7_av,table_clim2$CHELSA_bio10_7_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim2$CHELSA_bio10_8_av,table_clim2$CHELSA_bio10_8_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim2$CHELSA_bio10_9_av,table_clim2$CHELSA_bio10_9_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim2$CHELSA_bio10_10_av,table_clim2$CHELSA_bio10_10_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim2$CHELSA_bio10_11_av,table_clim2$CHELSA_bio10_11_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim2$CHELSA_bio10_12_av,table_clim2$CHELSA_bio10_12_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim2$CHELSA_bio10_13_av,table_clim2$CHELSA_bio10_13_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim2$CHELSA_bio10_14_av,table_clim2$CHELSA_bio10_14_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim2$CHELSA_bio10_15_av,table_clim2$CHELSA_bio10_15_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim2$CHELSA_bio10_16_av,table_clim2$CHELSA_bio10_16_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim2$CHELSA_bio10_17_av,table_clim2$CHELSA_bio10_17_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim2$CHELSA_bio10_18_av,table_clim2$CHELSA_bio10_18_1k,use="pairwise.complete.obs",method="spearman"),
                     cor(table_clim2$CHELSA_bio10_19_av,table_clim2$CHELSA_bio10_19_1k,use="pairwise.complete.obs",method="spearman"))

corr_clim_a <- cbind(variable,spearman_corr_a)

write.table(corr_clim_a,"CHELSA_clim/resampled_30k/corr_clim_a.txt")   


### LAND USE #################
# Land use harmonization (Couple Model Intercomparison Project, to link with climate models) 
#### Downloaded from: <http://luh.umd.edu/>, on 21st October 2018
####  The original information was in proportion of grid cell allocated for each use
####  Original resolution: 0.25 degrees
####  Original projection: WGS84
####  World

# Original. 0.25 degrees (~30 km) 
# Read original data (2015):  
lu1 <- raster("LUH/primf.grd")
lu2 <- raster("LUH/primn.grd")
lu3 <- raster("LUH/secdf.grd")
lu4 <- raster("LUH/secdn.grd")
lu5 <- raster("LUH/urban.grd")
lu6 <- raster("LUH/pastr.grd")
lu7 <- raster("LUH/range.grd")
# c3ann <- raster("LUH/c3ann.grd")
# c3per <- raster("LUH/c3per.grd")
# c3nfx <- raster("LUH/c3nfx.grd")
# c4ann <- raster("LUH/c4ann.grd")
# c4per <- raster("LUH/c4per.grd")
## Sum all crop categories:
# lu8 <- sum(c3ann,c3per,c3nfx,c4ann,c4per)  # Adding up all the crop categories
# # Save the new layer
# writeRaster(lu8,"LUH/all_crops.grd",overwrite=TRUE)
lu8 <- raster("LUH/all_crops.grd")


### WATER AVAILABILITY #######
# Rivers density ****
# Natural Earth (https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-rivers-lake-centerlines/) 
# Original: shapefile
# Manipulated in ArcGIS 
# Europe
# Final units: km/grid-cell 
rivers_eu <- raster("Hydro/rivers.tif")

# Water bodies ****
# Corine 2012
# Original: shapefile
# Manipulated in ArcGIS 
# Europe
# Final units: ha/grid-cell
wb_eu <- raster("Hydro/waterbodies.tif")

# Distance to the coast ****
# Mean distance from centroid to coastline 
# (Nielsen et al 2009).

# Inverse mask (1's=oceans; NA's=lands)
raster_inv <- raster::mask(is.na(ref_raster),ref_raster,maskvalue=1,updatevalue=NA)
# Calculate distance to nearest non-NA pixel. The distance unit is in meters if the RasterLayer is not projected (+proj=longlat) and in map units (typically also meters) when it is projected.
d <- distance(raster_inv)
# Optionally set non-land pixels to NA (otherwise values are "distance to non-land")
d <- d*ref_raster
writeRaster(d,"Coastline/coastline_30k.grd",format="raster",overwrite=TRUE, proj=TRUE)


### TOPOGRAPHY ####################
# Elevation ****
# 30 arc-sec (~1 km at the Equator) 
# elev_proj <- raster("GTOPO30/gtopo_mrg_cea.tif")
elev <- raster("GTOPO30/gtopo_mrg.tif")

# Resampled/averaged to 0.25 degrees (~30 km at the Equator) 
## Resampling...
# elev1 <- resample(elev,ref_raster,method="bilinear")
# elev1 <- raster::mask(elev1,ref_raster)
# writeRaster(elev1,"GTOPO30/elev_30k_r.tif",format="GTiff",prj=TRUE,overwrite=TRUE)
elev1 <- raster("GTOPO30/elev_30k_r.tif")

# Averaging...
# # Aggregate
# elev2 <- aggregate(elev,fact=30, fun=mean, expand=TRUE)
# # Allineate
# elev2 <- resample(elev2,ref_raster,method='bilinear')
# # Mask
# elev2 <- raster::mask(elev2,ref_raster)
# # Save in the folder "resampled10k":
# writeRaster(elev2,"GTOPO30/elev2_30k_a.tif",format="GTiff",prj=TRUE,overwrite=TRUE)
elev2 <- raster("GTOPO30/elev2_30k_a.tif")

## Extracting values from here and from the original data to compare values in a simple regression: 

# Random coordinates to extract from all data:
elev_10k_sample <- as.data.frame(sampleRandom(elev1,size=1000,xy=TRUE))
names(elev_10k_sample) <- c("x","y","elev_10k_r")

coord_sample <- elev_10k_sample[,c(1:2)]

elev_1k_sample <- raster::extract(elev,coord_sample,df=TRUE)
names(elev_1k_sample) <- paste(colnames(elev_1k_sample),"1k",sep="_")
elev_av_sample <- raster::extract(elev2,coord_sample,df=TRUE)
names(elev_av_sample) <- c("ID","gtopo30k_av")

table_elev <- cbind(elev_10k_sample,elev_1k_sample)
table_elev2 <- cbind(elev_av_sample,elev_1k_sample)

# Correlations:
variable <- "variable"

spearman_corr <- cor(table_elev$elev_10k_r,table_elev$gtopo_mrg_1k,use="pairwise.complete.obs",method="spearman")
corr_elev <- cbind(variable,spearman_corr)
write.table(corr_elev,"GTOPO30/corr_elev_r.txt")

spearman_corr2 <- cor(table_elev2$gtopo30k_av,table_elev2$gtopo_mrg_1k,use="pairwise.complete.obs",method="spearman")
corr_elev2 <- cbind(variable,spearman_corr2)
write.table(corr_elev2,"GTOPO30/corr_elev_a.txt")


# SD of the elevation ****
elev_sd <- raster::aggregate(elev, fact=30,fun=sd,expand=TRUE)
# Align:
elev_sd <- resample(elev_sd, ref_raster,method="bilinear")
# Mask:
elev_sd <- raster::mask(elev_sd,ref_raster)
# Save:
writeRaster(elev_sd,"GTOPO30/elev_sd_30k.tif",format="GTiff",prj=TRUE,overwrite=TRUE)


### ACCESIBILITY ##################
# Travel time to a location of interest using land (road/off road) or water (navigable river, lake and ocean) based travel. city=50,000 people. Available at: <https://forobs.jrc.ec.europa.eu/products/gam/download.php>

# 30 arc-sec (~1 km at the Equator; original) 
acc <- raster("Accesibility/acc_50k.tif")

# 0.25 degrees (~30 km at the Equator):
# Resampling....
acc1 <- resample(acc,ref_raster, method='bilinear')
acc1 <- raster::mask(acc1,ref_raster)
writeRaster(acc1, "Accesibility/acc_30k_r.tif",format="GTiff",prj=TRUE,overwrite=TRUE)

# Averaging...
# Aggregate
acc2 <- aggregate(acc,fact=30, fun=mean, expand=TRUE)
# Allineate
acc2 <- resample(acc2,ref_raster,method='bilinear') 
# Mask
acc2 <- raster::mask(acc2,ref_raster)
# Save in the folder "resampled10k":
writeRaster(acc2,"Accesibility/acc_30k_a.tif",format="GTiff",prj=TRUE,overwrite=TRUE)

## Extracting values from here and from the original data to compare values in a simple regression: 

# Random coordinates to extract from all data:
acc_30k_sample_r <- as.data.frame(sampleRandom(acc1,size=1000,xy=TRUE))
names(acc_30k_sample_r) <- c("x","y","acc_30k_r")

coord_sample <- acc_30k_sample_r[,c(1:2)]

acc_1k_sample <- raster::extract(acc,coord_sample,df=TRUE)
names(acc_1k_sample)
acc_30k_sample_a <- raster::extract(acc2,coord_sample,df=TRUE)
names(acc_30k_sample_a) <- c("ID","acc_30k_a")

table_acc <- cbind(acc_30k_sample_r,acc_1k_sample)
table_acc2 <- cbind(acc_30k_sample_a,acc_1k_sample)

# Correlations:
variable <- "variable"

spearman_corr <- cor(table_acc$acc_30k_r,table_acc$acc_50k,use="pairwise.complete.obs",method="spearman")
spearman_corr2 <- cor(table_acc2$acc_30k_a,table_acc2$acc_50k,use="pairwise.complete.obs",method="spearman")

corr_acc <- cbind(variable,spearman_corr)
corr_acc2 <- cbind(variable,spearman_corr2)

write.table(corr_acc,"Accesibility/corr_acc_r.txt")
write.table(corr_acc,"Accesibility/corr_acc_a.txt")


# MODELING ######################################################

## 1. Handling SDM pseudosabsences #############################
## Global models ####'''''''''''''''''''#########################

### Predictors: ----------------------------------------------------------------
# Check and remove correlation:
# nocorrvar <- vifstep(clim_30k,th=4)
# predictors <- as.character(nocorrvar@results$Variables)  # --> this gives me the names of the remaining variables
# 
# # Selection of only noncorrelated predictors from the full environmental set:
# clim_sub <- raster::subset(clim_30k,c(predictors))
# clim_sub <- stack(clim_sub)
# # 
# # Save to retrieve later:
# stackSave(clim_sub,"CHELSA_clim/resampled_30k/clim_sub")
clim_sub <- stackOpen("CHELSA_clim/resampled_30k/clim_sub")

### Set parameters for the model: -----------------------------------------------
name_model_folder <- "GLOBAL2_30k" # GLOBAL2_30k # GLOBAL1_30k 
model_class <- "global" 
dataset_type <- "_presences_cert_30k" # _presences_30k_cert_na ## from GBIF

### Fit & predict: ---------------------------------------------------------------

# Individual models: -------------------------------------------------------------
my_models <- c("GLM","GAM","MAXENT.Phillips","FDA","GBM") # Algorithms to fit the models   
my_pa_strategy <- c("random") # Strategy to select PAs   
my_pa_dist_min <- 0 # Distance where to start throwing PAs  
my_runs <- 3 # Number of PAs runs
my_n_pa <- 20000 # Number of pseudoabsences, when not dependent of no. presences
my_proj_name <- "global"  # Name for the output projections

for (i in 1:length(list_sp)) {
  
  setwd("D:/Postdoc_Uppsala")
  # Species name:
  sp <- list_sp[i]
  sp2 <- sub(" ","_",sp)
  sp3 <- sub("_",".",sp2)
  
  # Presences: ------------------------------------------------------------
  file_name1 <- paste("Occ_vs_pres/",sp2,dataset_type,".shp",sep="")
  sp_pres <- shapefile(file_name1)
  n_mod <- dim(sp_pres)[1]
  
  # Pseudoabsences: -------------------------------------------------------
  # Select number of pseudoabsences:
  n_pa <- my_n_pa
  
  # Format data for 'BIOMOD': ---------------------------------------------
  
  # Presences
  mypres <- rep(1,dim(sp_pres)[1])
  
  # Presences coordinates
  mycoord <- coordinates(sp_pres)
  
  # Clean data
  clean_data <- BIOMOD_FormatingData(resp.var=mypres, expl.var=clim_sub, resp.xy=mycoord,resp.name=sp2, PA.nb.rep=my_runs,PA.nb.absences=n_pa,PA.strategy = my_pa_strategy, PA.dist.min=my_pa_dist_min)
  
  # Settings for the model algorithms
  model_opt <- BIOMOD_ModelingOptions(GLM=list(type='quadratic', interaction.level=0),GBM=list(n.trees=1000),MAXENT.Phillips = list(path_to_maxent.jar = "D:/maxent.jar",product=FALSE),GAM=list(k=3))
  
  my_dir <- paste("D:/Postdoc_Uppsala/",name_model_folder, sep="")
  setwd(my_dir)
  
  # Run the models: ---------------------------------------------------------
  # Separate
  all_models <- BIOMOD_Modeling(data=clean_data, models=my_models, models.options=model_opt, NbRunEval = 4, DataSplit = 70, VarImport = 3, do.full.models = F, modeling.id = "global")
  
  # Project the models:---------------------------------------------------------
  # Separate
  models_proj_current <- BIOMOD_Projection(modeling.output = all_models,new.env = clim_sub,proj.name = my_proj_name ,binary.meth = "TSS",do.stack = FALSE )
}

# Ensemble model: -------------------------------------------------------------
for (i in 1:length(list_sp)) {
  
  setwd("D:/Postdoc_Uppsala")
  # Species name:
  sp <- list_sp[i]
  sp2 <- sub(" ","_",sp)
  sp3 <- sub("_",".",sp2)
  
  # Load separately fit models:
  model_folder <- paste("D:/Postdoc_Uppsala/",name_model_folder,sep="")
  setwd(model_folder)
  
  name_model <- paste(sp3,"/",sp3,".",model_class,".models.out",sep="")
  model_i <- load(name_model)
  
  # Load separate models projections:
  name_proj <- paste(sp3,"/proj_",model_class,"/",sp3,".",model_class,".projection.out",sep="")
  models_proj_current <- load(name_proj)
  
  # Fit ensemble model
  all_test <- get_evaluations(get(model_i))
  all_tss <- all_test[2, 1, 1:5, 1:4, 1:3]
  
  if (all_tss>=0.7) {
    all_ensemble_model <- BIOMOD_EnsembleModeling(modeling.output = get(model_i), em.by='all',eval.metric = 'TSS', eval.metric.quality.threshold = 0.7, models.eval.meth = c('KAPPA','TSS','ROC'), prob.mean=FALSE, prob.cv=TRUE, committee.averaging=TRUE, prob.mean.weight = TRUE, VarImport = 3)
  }
  
  else {
    val_tss <- stats::quantile(all_tss, probs = 0.9,names=FALSE,na.rm=TRUE)
    all_ensemble_model <- BIOMOD_EnsembleModeling(modeling.output = get(model_i), em.by='all',eval.metric = 'TSS', eval.metric.quality.threshold = val_tss, models.eval.meth = c('KAPPA','TSS','ROC'), prob.mean=FALSE, prob.cv=TRUE, committee.averaging=TRUE, prob.mean.weight = TRUE, VarImport = 3)
  }
  
  # Project ensemble model
  models_ensemble_proj_current <- BIOMOD_EnsembleForecasting(EM.output = all_ensemble_model,projection.output = get(models_proj_current),binary.meth = "TSS",do.stack = FALSE )
  
} 
  
### Accuracy ------------------------------------------------------------------
for (i in 1:length(list_sp)) {

  setwd("D:/Postdoc_Uppsala")
  # Species name:
  sp <- list_sp[i]
  sp2 <- sub(" ","_",sp)
  sp3 <- sub("_",".",sp2)
  
  # Load separated models:
  model_folder <- paste("D:/Postdoc_Uppsala/",name_model_folder,sep="")
  setwd(model_folder)
  
  name_model <- paste(sp3,"/",sp3,".",model_class,".models.out",sep="")
  model_i <- load(name_model)
  
  # Load unique ensemble model:
  name_ensemble <- paste(sp3,"/",sp3,".",model_class,"ensemble.models.out",sep="")
  ensemble_i <- load(name_ensemble)
  
  # Get evaluations of ensemble models to decide which one to choose:
  all_ensemble_models_scores <- get_evaluations(get(ensemble_i))
  file_name3 <- paste(sp3,"/",sp3,"_eval_ensemble.csv",sep="")
  write.table(all_ensemble_models_scores,file_name3,sep=";")
  
  # Get variables importance for all models
  all_models_var_import <- get_variables_importance(get(model_i))
  # Calculate the mean of variable importance by algorithm
  table_importance <- apply(all_models_var_import, c(1,2), mean, na.rm=TRUE)
  file_name2 <- paste(sp3,"/",sp3,"_importance_var.csv",sep="")
  write.table(table_importance,file_name2,sep=";")
  
}

### Spatial dissimilarity - Bhattacharyya distance (BD) - following Huang et al (2018) - Values of Bhattacharyya distance of 0 indicate identical predictions, values of 1 indicate the most different predictions ---------------------------

# Statistical measure of dissimilarity in a raster cell for the Bhattacharyya coefficient (eq. A.2 in the paper):
overlap <- function(A,B) { sqrt(A*B)+sqrt((1-A)*(1-B))}  # Eq. A.2

# Eq. A.5 without the summatory (in the paper)
Distance <- function(A,B) { -log(overlap(A,B)) }  # Eq. A.5

# BD per pixel (result --> one raster per species, and per global/regional) & table average per species
# Models to compare
model1 <- "GLOBAL1_30k"
model2 <- "GLOBAL2_30k"

my_proj_name <- "global"

table_disim_r <- c()

for (i in 1:length(list_sp)) {
  setwd("D:/Postdoc_Uppsala")
  
  ## Species: 
  sp <- list_sp[i]
  sp2 <- sub(" ","_",sp)
  sp3 <- sub(" ",".",sp)
  
  # Load models:
  file_mod1 <- paste(model1,"/",sp3,"/proj_", my_proj_name, "/individual_projections/",sp3,"_EMcaByTSS_mergedAlgo_mergedRun_mergedData.grd",sep="")
  raster_mod1 <- raster(file_mod1)
  raster_mod1_table <- raster::as.data.frame(raster_mod1/1000,xy=TRUE)
  
  # Continuous global cert+NA prediction
  file_mod2 <- paste(model2,"/",sp3,"/proj_", my_proj_name, "/individual_projections/",sp3,"_EMcaByTSS_mergedAlgo_mergedRun_mergedData.grd",sep="")
  raster_mod2 <- raster(file_mod2)
  raster_mod2_table <- raster::as.data.frame(raster_mod2/1000,xy=TRUE)
  
  # With this function I should get values between 0 and 1 (0=identical; 1=opposed)
  distance_preds <- Distance(raster_mod1_table$layer,raster_mod2_table$layer)
  distance_preds[is.infinite(distance_preds)] <- 1 
  
  # Table with coordinates
  table_dist <- as.data.frame(cbind(coordinates(raster_mod1),distance_preds))
  table_dist <- table_dist[!is.infinite(rowSums(table_dist)),]  # to remove infinite values
  mean_sp_r <- mean(table_dist$distance_preds,na.rm=TRUE)
  line2 <- c(sp2, mean_sp_r)
  
  table_disim_r <- rbind(table_disim_r,line2)
  
  # Create a raster with coordinates
  distance_preds_r <- rasterFromXYZ(table_dist)
  
  # Save:
  filey <- paste("D:/Postdoc_Uppsala/Dissimilarity/",sp2,"_",model1,model2,".grd",sep="")
  writeRaster(distance_preds_r,filey,format="raster",prj=TRUE,overwrite=TRUE)
  
}    

# Build table and save:
table_disim_r_clean <- as.data.frame(table_disim_r)
names(table_disim_r_clean) <- c("Species","Dissimilarity")

heading_mammals <- as.data.frame(cbind("Mammals",""))
names(heading_mammals) <- c("Species","Dissimilarity")

heading_birds <- as.data.frame(cbind("Birds",""))
names(heading_birds) <- c("Species","Dissimilarity")

heading_amph <- as.data.frame(cbind("Amphibians",""))
names(heading_amph) <- c("Species","Dissimilarity")

heading_reptiles <- as.data.frame(cbind("Reptiles",""))
names(heading_reptiles) <- c("Species","Dissimilarity")

table_disim_r_clean <- rbind(heading_mammals,table_disim_r_clean[1:9,],heading_birds,table_disim_r_clean[10:13,],heading_amph,table_disim_r_clean[14,],heading_reptiles,table_disim_r_clean[15,])

name1 <- paste("D:/Postdoc_Uppsala/Dissimilarity/",model1,model2,".txt",sep="")
write.table(table_disim_r_clean, name1,row.names=FALSE)

## 2. SDM ######################################################
## European models ####'''''''''''''''''''#########################

## Predictors: ----------------------------------------------------------------
# 1 # Climate: (averaged)
clim_30k <- stack(list.files(path = "CHELSA_clim/resampled_30k/", pattern = "_noproj_av.grd$", full.names=TRUE))
clim_30k <- raster::crop(clim_30k,eu_raster)
clim_30k <- raster::mask(clim_30k,eu_raster)

# 2 # Land use: (original)
lu1 <- raster("LUH/primf.grd")
lu2 <- raster("LUH/primn.grd")
lu3 <- raster("LUH/secdf.grd")
lu4 <- raster("LUH/secdn.grd")
lu5 <- raster("LUH/urban.grd")
lu6 <- raster("LUH/pastr.grd")
lu7 <- raster("LUH/range.grd")
lu8 <- raster("LUH/all_crops.grd")

lu <- stack(lu1,lu2,lu3,lu4,lu5,lu6,lu7,lu8)
lu <- raster::crop(lu,eu_raster)
lu <- raster::mask(lu,eu_raster)
names(lu) <- c("forest_prim","non_f_prim","forest_sec","non_f_sec","urban","pasture","rangeland","crops")

# 3 # Rivers:
rivers_eu <- raster("Hydro/rivers_europe.tif")
rivers_eu <- raster::crop(rivers_eu,eu_raster)
rivers_eu <- raster::mask(rivers_eu,eu_raster)

# 4 # Water bodies:
wb_eu <- raster("Hydro/waterbodies_europe.tif")
wb_eu <- raster::crop(wb_eu,eu_raster)
wb_eu <- raster::mask(wb_eu,eu_raster)

# 5 # Coastline:
coastline <- raster("Coastline/coastline_30k.grd") # meters
coastline <- coastline/1000 # kilometers
coastline <- raster::crop(coastline,eu_raster)
coastline <- raster::mask(coastline,eu_raster)
names(coastline) <- "coast"

# 6 # Elevation (average):
elev2 <- raster("GTOPO30/elev2_30k_a.tif")
elev2 <- raster::crop(elev2,eu_raster)
elev2 <- raster::mask(elev2,eu_raster)

# 7 # SD of the elevation:
elev_sd <- raster("GTOPO30/elev_sd_30k.tif")
elev_sd <- raster::crop(elev_sd,eu_raster)
elev_sd <- raster::mask(elev_sd,eu_raster)

# 8 # Accesibility: (average)
acc <- raster("Accesibility/acc_30k_a.tif")
acc <- raster::crop(acc,eu_raster)
acc <- raster::mask(acc,eu_raster)

# Stacking all:
all_pred1 <- stack(clim_30k,lu,rivers_eu, wb_eu,coastline,elev2,elev_sd,acc)

# # Check and remove correlation (avoid running every time, stochasticity may add/remove one variable in the limit):
# nocorrvar <- vifstep(all_pred1,th=4) 
# predictors <- as.character(nocorrvar@results$Variables)  # --> this gives me the names of the remaining variables
# save(predictors,file="D:/Postdoc_Uppsala/predictors.r")
load("D:/Postdoc_Uppsala/predictors.r")

# Selection of only noncorrelated predictors from the full environmental set:
pred_eu_sub <- raster::subset(all_pred1,c(predictors))


### Set parameters for the model: -----------------------------------------------
name_model_folder <- "REGIONAL2_30k" # REGIONAL1_30k
model_class <- "regional" 
dataset_type <- "_presences_cert_na_30k" # _presences_30k_cert_na ## from GBIF

### Fit & predict: ---------------------------------------------------------------

# Individual models: -------------------------------------------------------------
my_models <- c("GLM","GAM","MAXENT.Phillips","FDA","GBM") # Algorithms to fit the models   
my_pa_strategy <- c("random") # Strategy to select PAs   
my_pa_dist_min <- 0 # Distance where to start throwing PAs  
my_runs <- 3 # Number of PAs runs
my_n_pa <- 5000 # Number of pseudoabsences, when not dependent of no. presences
my_proj_name <- "regional"  # Name for the output projections

for (i in 1:length(list_sp)) {
  
  setwd("D:/Postdoc_Uppsala")
  # Species name:
  sp <- list_sp[i]
  sp2 <- sub(" ","_",sp)
  sp3 <- sub("_",".",sp2)
  
  # Presences: ------------------------------------------------------------
  file_name1 <- paste("Occ_vs_pres/",sp2,dataset_type,".shp",sep="")
  sp_pres <- shapefile(file_name1)
  # Subset of Europe:
  sp_pres_eu <- rasterize(sp_pres,eu_raster,1)
  n_mod <- summary(as.factor(sp_pres_eu@data@values))[1]
  
  # Pseudoabsences: -------------------------------------------------------
  # Select number of pseudoabsences:
  n_pa <- my_n_pa
  
  # Weights based on the global model:
  file_pa <- paste("D:/Postdoc_Uppsala/GLOBAL2_30k/",sp3,"/proj_global/individual_projections/",sp3,"_EMcaByTSS_mergedAlgo_mergedRun_mergedData.grd",sep="")
  all_pa <- raster(file_pa)
  pa_eu <- raster::crop(all_pa,eu_raster)
  pa_eu <- raster::mask(pa_eu,eu_raster)
  pa_eu_1 <- pa_eu/1000 # conversion from 0 to 1
  weights_pa <- 1/(1+((pa_eu_1/(pa_eu_1-1))^2)) # inverse logistic transformation
  
  # Format data for 'BIOMOD': ---------------------------------------------
  # Subset of only presences:
  subset_1 <- as.data.frame(rasterToPoints(sp_pres_eu,fun=function(x){x==1}))
  
  # Presences
  mypres <- rep(1,dim(subset_1)[1])
  
  # Presences coordinates
  mycoord <- subset_1[,c("x","y")]
 
  # Clean data
  clean_data <- BIOMOD_FormatingData(resp.var=mypres, expl.var=pred_eu_sub, resp.xy=mycoord,resp.name=sp2, PA.nb.rep=my_runs,PA.nb.absences=n_pa,PA.strategy = my_pa_strategy, PA.dist.min=my_pa_dist_min)
  
  # Adding "weights" to the data
  # Get the presences + PA dataset
  my_pres_PA_df <- data.frame(clean_data@coord, obs = clean_data@data.species, clean_data@PA)
  
  # Add the weigth vector to the presences
  pres <- subset(my_pres_PA_df, obs==1)
  pres$yweights <- 1
  # Weight of absences based on global model
  abs <- subset(my_pres_PA_df, is.na(obs)==TRUE)
  weight_values_ab <- raster::extract(weights_pa,cbind(abs$x,abs$y))
  abs$yweights <- weight_values_ab
  
  new_my_pres_PA_df <- rbind(pres,abs)
  
  # Settings for the model algorithms
  model_opt <- BIOMOD_ModelingOptions(GLM=list(type='quadratic', interaction.level=0),GBM=list(n.trees=1000),MAXENT.Phillips = list(path_to_maxent.jar = "C:/Users/espa0001/maxent.jar",product=FALSE),GAM=list(k=3))
  
  my_dir <- paste("D:/Postdoc_Uppsala/",name_model_folder, sep="")
  setwd(my_dir)
  
  # Run the models: ---------------------------------------------------------
  # Separate
  all_models <- BIOMOD_Modeling(data=clean_data, models=my_models, models.options=model_opt, NbRunEval = 4, DataSplit = 70, VarImport = 3, do.full.models = F, modeling.id = "regional", Yweights=new_my_pres_PA_df$yweights)
  
  # Project the models:---------------------------------------------------------
  # Separate
  models_proj_current <- BIOMOD_Projection(modeling.output = all_models,new.env = pred_eu_sub,proj.name = my_proj_name ,binary.meth = "TSS",do.stack = FALSE )
}

# Ensemble model: -------------------------------------------------------------
for (i in 1:length(list_sp)) {
  
  setwd("D:/Postdoc_Uppsala")
  # Species name:
  sp <- list_sp[i]
  sp2 <- sub(" ","_",sp)
  sp3 <- sub("_",".",sp2)
  
  # Load separately fit models:
  model_folder <- paste("D:/Postdoc_Uppsala/",name_model_folder,sep="")
  setwd(model_folder)
  
  name_model <- paste(sp3,"/",sp3,".",model_class,".models.out",sep="")
  model_i <- load(name_model)
  
  # Load separate models projections:
  name_proj <- paste(sp3,"/proj_",model_class,"/",sp3,".",model_class,".projection.out",sep="")
  models_proj_current <- load(name_proj)
  
  # Fit ensemble model
  all_test <- get_evaluations(get(model_i))
  all_tss <- all_test[2, 1, 1:5, 1:4, 1:3]
  
  if (all_tss>=0.7) {
    all_ensemble_model <- BIOMOD_EnsembleModeling(modeling.output = get(model_i), em.by='all',eval.metric = 'TSS', eval.metric.quality.threshold = 0.7, models.eval.meth = c('KAPPA','TSS','ROC'), prob.mean=FALSE, prob.cv=TRUE, committee.averaging=TRUE, prob.mean.weight = TRUE, VarImport = 3)
  }
  
  else {
    val_tss <- stats::quantile(all_tss, probs = 0.9,names=FALSE,na.rm=TRUE)
    all_ensemble_model <- BIOMOD_EnsembleModeling(modeling.output = get(model_i), em.by='all',eval.metric = 'TSS', eval.metric.quality.threshold = val_tss, models.eval.meth = c('KAPPA','TSS','ROC'), prob.mean=FALSE, prob.cv=TRUE, committee.averaging=TRUE, prob.mean.weight = TRUE, VarImport = 3)
  }
  
  # Project ensemble model
  models_ensemble_proj_current <- BIOMOD_EnsembleForecasting(EM.output = all_ensemble_model,projection.output = get(models_proj_current),binary.meth = "TSS",do.stack = FALSE )
  
} 

### Accuracy ------------------------------------------------------------------
for (i in 1:length(list_sp)) {
  
  setwd("D:/Postdoc_Uppsala")
  # Species name:
  sp <- list_sp[i]
  sp2 <- sub(" ","_",sp)
  sp3 <- sub("_",".",sp2)
  
  # Load separated models:
  model_folder <- paste("D:/Postdoc_Uppsala/",name_model_folder,sep="")
  setwd(model_folder)
  
  name_model <- paste(sp3,"/",sp3,".",model_class,".models.out",sep="")
  model_i <- load(name_model)
  
  # Load unique ensemble model:
  name_ensemble <- paste(sp3,"/",sp3,".",model_class,"ensemble.models.out",sep="")
  ensemble_i <- load(name_ensemble)
  
  # Get evaluations of ensemble models to decide which one to choose:
  all_ensemble_models_scores <- get_evaluations(get(ensemble_i))
  file_name3 <- paste(sp3,"/",sp3,"_eval_ensemble.csv",sep="")
  write.table(all_ensemble_models_scores,file_name3,sep=";")
  
  # Get variables importance for all models
  all_models_var_import <- get_variables_importance(get(model_i))
  # Calculate the mean of variable importance by algorithm
  table_importance <- apply(all_models_var_import, c(1,2), mean, na.rm=TRUE)
  file_name2 <- paste(sp3,"/",sp3,"_importance_var.csv",sep="")
  write.table(table_importance,file_name2,sep=";")
}

### Spatial dissimilarity - Bhattacharyya distance (BD) - following Huang et al (2018) - Values of Bhattacharyya distance of 0 indicate identical predictions, values of 1 indicate the most different predictions ---------------------------

# Statistical measure of dissimilarity in a raster cell for the Bhattacharyya coefficient (eq. A.2 in the paper):
overlap <- function(A,B) { sqrt(A*B)+sqrt((1-A)*(1-B))}  # Eq. A.2

# Eq. A.5 without the summatory (in the paper)
Distance <- function(A,B) { -log(overlap(A,B)) }  # Eq. A.5

# BD per pixel (result --> one raster per species, and per global/regional) & table average per species
# Models to compare
model1 <- "REGIONAL5_30k"
model2 <- "REGIONAL2_30k"

my_proj_name <- "regional"

table_disim_r <- c()

for (i in 1:length(list_sp)) {
  setwd("D:/Postdoc_Uppsala")
  
  ## Species: 
  sp <- list_sp[i]
  sp2 <- sub(" ","_",sp)
  sp3 <- sub(" ",".",sp)
  
  # Load models:
  file_mod1 <- paste(model1,"/",sp3,"/proj_", my_proj_name, "/individual_projections/",sp3,"_EMcaByTSS_mergedAlgo_mergedRun_mergedData.grd",sep="")
  raster_mod1 <- raster(file_mod1)
  raster_mod1_table <- raster::as.data.frame(raster_mod1/1000,xy=TRUE)
  
  # Continuous global cert+NA prediction
  file_mod2 <- paste(model2,"/",sp3,"/proj_", my_proj_name, "/individual_projections/",sp3,"_EMcaByTSS_mergedAlgo_mergedRun_mergedData.grd",sep="")
  raster_mod2 <- raster(file_mod2)
  raster_mod2_table <- raster::as.data.frame(raster_mod2/1000,xy=TRUE)
  
  # With this function I should get values between 0 and 1 (0=identical; 1=opposed)
  distance_preds <- Distance(raster_mod1_table$layer,raster_mod2_table$layer)
  distance_preds[is.infinite(distance_preds)] <- 1 
  
  # Table with coordinates
  table_dist <- as.data.frame(cbind(coordinates(raster_mod1),distance_preds))
  table_dist <- table_dist[!is.infinite(rowSums(table_dist)),]  # to remove infinite values
  mean_sp_r <- mean(table_dist$distance_preds,na.rm=TRUE)
  line2 <- c(sp2, mean_sp_r)
  
  table_disim_r <- rbind(table_disim_r,line2)
  
  # Create a raster with coordinates
  distance_preds_r <- rasterFromXYZ(table_dist)
  
  # Save:
  filey <- paste("D:/Postdoc_Uppsala/Dissimilarity/",sp2,"_",model1,model2,".grd",sep="")
  writeRaster(distance_preds_r,filey,format="raster",prj=TRUE,overwrite=TRUE)
  
}    

# Build table and save:
table_disim_r_clean <- as.data.frame(table_disim_r)
names(table_disim_r_clean) <- c("Species","Dissimilarity")

heading_mammals <- as.data.frame(cbind("Mammals",""))
names(heading_mammals) <- c("Species","Dissimilarity")

heading_birds <- as.data.frame(cbind("Birds",""))
names(heading_birds) <- c("Species","Dissimilarity")

heading_amph <- as.data.frame(cbind("Amphibians",""))
names(heading_amph) <- c("Species","Dissimilarity")

heading_reptiles <- as.data.frame(cbind("Reptiles",""))
names(heading_reptiles) <- c("Species","Dissimilarity")

table_disim_r_clean <- rbind(heading_mammals,table_disim_r_clean[1:9,],heading_birds,table_disim_r_clean[10:13,],heading_amph,table_disim_r_clean[14,],heading_reptiles,table_disim_r_clean[15,])

name1 <- paste("D:/Postdoc_Uppsala/Dissimilarity/",model1,model2,".txt",sep="")
write.table(table_disim_r_clean, name1,row.names=FALSE)


# MULTI-SPECIES SUMMARY #########################################

## Total Richness -----------------------------------------------------------
setwd("D:/Postdoc_Uppsala")

name_model_folder <- "REGIONAL2_30k" 
model_class <- "regional" # global / regional


all_bin <- stack()

for (i in 1:length(list_sp)) {
  # Species: 
  sp <- list_sp[i]
  sp2 <- sub(" ","_",sp)
  sp3 <- sub(" ",".",sp)
  
  # load binary projection
  file_bin <- paste(name_model_folder,"/",sp3,"/proj_",model_class,"/individual_projections/",sp3,"_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.grd",sep="")
  bin_sp <- raster(file_bin)
  
  all_bin <- stack(all_bin, bin_sp)
}

rich_sum <- stackApply(all_bin,indices=15,fun=sum)
rich_sum <- raster::crop(rich_sum,eu_raster)
rich_sum <- raster::mask(rich_sum,eu_raster)

name_rich <- paste("D:/Postdoc_Uppsala/Maps/Richness",name_model_folder,sep="_")
writeRaster(rich_sum,filename=name_rich,overwrite=TRUE)


## CV Europe predictions ----------------------------------------------------
name_model_folder <- "REGIONAL2_30k" 
model_class <- "regional" # global / regional
raster_to_use <- eu_raster  # ref_raster / eu_raster

# Rescaling & getting one map per species + calculate mean per species (table)
raster_list <- list()

for (i in 1:length(list_sp)) {
  
  sp <- list_sp[i]
  sp2 <- sub(" ","_",sp)
  sp3 <- sub("_",".",sp2)
  
  file_name <- paste(name_model_folder,"/",sp3,"/proj_", model_class,"/individual_projections/",sp3,"_EMcvByTSS_mergedAlgo_mergedRun_mergedData.grd",sep="")
  raster_sp <- raster(file_name)
  
  r.min <- cellStats(raster_sp, "min")
  r.max <- cellStats(raster_sp, "max")
  
  raster_sp_scale <- (raster_sp - r.min) / (r.max - r.min) 
  raster_sp_scale <- raster::mask(raster_sp_scale,raster_to_use)
  raster_list[[i]] <- raster_sp_scale
  
  # Saving the rescaled version for each species to map it:
  file_cv <- paste("Maps/",sp2,name_model_folder,"CV", sep="_")
  writeRaster(raster_sp_scale,file_cv,prj=TRUE,overwrite=TRUE)
}

# Calculating the mean per grid-cell
raster_stack <- stack(raster_list)
cv_prd <- stackApply(raster_stack,indices=15,fun=mean,na.rm=TRUE)
cv_prd <- raster::mask(cv_prd,raster_to_use)

# Scale again to obtain values between 0 and 1
r.min <- cellStats(cv_prd, "min")
r.max <- cellStats(cv_prd, "max")
cv_prd_scale <- (cv_prd - r.min) / (r.max - r.min)

# Writing resulting raster
name1 <- paste("Maps/Mean_cv",name_model_folder,"CV", sep="_")
writeRaster(cv_prd_scale,name1,prj=TRUE,overwrite=TRUE)


## Priority management strategy zones ---------------------------------------
## STEPS 
name_model_folder_r <- "REGIONAL2_30k"
name_model_folder_g <- "GLOBAL2_30k"

# (1) Raster layers to include in the index:
## Richness of IATV (EU)
name_raster_eu <- paste("D:/Postdoc_Uppsala/Maps/Richness",name_model_folder_r,sep="_")
rich_sum <- raster(name_raster_eu)

## Average CV of the predictions 
name_raster_cv_eu <- paste("D:/Postdoc_Uppsala/Maps/Mean_cv",name_model_folder_r,"CV",sep="_")
cv_prd_scale <- raster(name_raster_cv_eu)

## Richness of IATV (global)
name_raster_g <- paste("D:/Postdoc_Uppsala/Maps/Richness",name_model_folder_g,sep="_")
rich_g_sum <- raster(name_raster_g)

# (2) Reclassify each raster into 2 classes:
## Matrices of reclassification  -----    ----   -----   -----   ------  --- 1=low / 2=high
class_r <- c(0,8,1, 8,16,2)
rclmat_r <- matrix(class_r, ncol=3,byrow=TRUE)

class_cv <- c(0,0.5,1, 0.5,1,2)
rclmat_cv <- matrix(class_cv,ncol=3,byrow=TRUE)

## New binary rasters:
richness2 <- reclassify(rich_sum, rclmat_r,right=FALSE)
cv_prd2 <- reclassify(cv_prd_scale,rclmat_cv,right=FALSE)
richness_g2 <- reclassify(rich_g_sum, rclmat_r,right=FALSE)
# Writing resulting rasters (not compulsory):
#
# name3 <- paste("Maps/Tot_richBIN",model_name, sep="_")
# writeRaster(richness2,name3,format="raster",prj=TRUE,overwrite=TRUE)
# 
# name4 <- paste("Maps/Mean_CVBIN",model_name,sep="_")
# writeRaster(cv_prd2,name4,format="raster",prj=TRUE,overwrite=TRUE)
# 
# name6 <- paste("Maps/Tot_rich_g_BIN",model_name,sep="_")
# writeRaster(richness_g2,name6,format="raster",prj=TRUE,overwrite=TRUE)

## Convert to data.frame to selece according to criteria
richness2_table <- as.data.frame(richness2, xy = TRUE)
names(richness2_table) <- c("x","y","richness")

cv_prd2_table <- as.data.frame(cv_prd2, xy = TRUE)
names(cv_prd2_table) <- c("x","y","cv")

richness_g2_table <- as.data.frame(richness_g2, xy = TRUE)
names(richness_g2_table) <- c("x","y","richness_g")

# Joining all:
table1 <- merge(richness2_table,cv_prd2_table,by=c("x","y"))
table1 <- merge(table1,richness_g2_table)

# (3) Selecting combinations:
hotspot <- subset(table1,cv==1 & richness==2)
hotspot$zone <- 1 # A

coldspot <- subset(table1, cv==1 & richness==1)
coldspot$zone <- 2 # B

dd_highprio <- subset(table1, cv==2 & richness==2 & richness_g==2)
dd_highprio$zone <- 3 # C

dd_lowprio <- subset(table1, cv==2 & richness==1 & richness_g==1)
dd_lowprio$zone <- 4 # F

adapt <- subset(table1, cv==2 & richness==2 & richness_g==1)
adapt$zone <- 5 # D

colon <- subset(table1, cv==2 & richness==1 & richness_g==2)
colon$zone <- 6 # E

# One table
table_zones <- rbind(hotspot, coldspot, dd_highprio, dd_lowprio, adapt, colon)

## Convert to a raster
zones_raster <- rasterize(x=cbind(table_zones$x,table_zones$y),y=eu_raster,field=table_zones$zone,update=TRUE)

## Counting cells:
table_count <- summary(as.factor(zones_raster@data@values))
name7 <- paste("Maps/Zones",name_model_folder_r,name_model_folder_g,".txt",sep="")
write.table(table_count,name7)

## Save as raster: 
writeRaster(zones_raster,name7,format="raster",prj=TRUE,overwrite=TRUE)


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ---------------------------- RESULTS ---------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

# TABLE: MODEL ACCURACY ########################################################
setwd("D:/Postdoc_Uppsala")

# Here you decide which results you want the output from. First line is the name of the folder where the output of the model is, second is the quality of GBIF data from which it is fit, third indicates if it included easin data (only for European models), forth indicates if it is a global model or not.
name_model_folder <- "REGIONAL2_30k"
dataset_type <- "_presences_cert_na_30k" # Options: _presences_cert_30k / _presences_30k_cert_na ## from GBIF
global <- "no"
model_class <- "regional"
raster_to_use <- eu_raster # eu_raster # ref_raster

## Sample size (n)
# empty vector to store n for all species
col_n <- c()

for (i in 1:length(list_sp)) {
  
  # Species: 
  sp <- list_sp[i]
  sp2 <- sub(" ","_",sp)
  sp3 <- sub(" ",".",sp)
  
  if (global=="no"){raster1 <- eu_raster}
  else {raster1 <- ref_raster}
  
  # GBIF points
  file_name1 <- paste("Occ_vs_pres/",sp2,dataset_type,".shp",sep="")
  sp_pres <- shapefile(file_name1)
  # rasterize
  points_mod_r <- rasterize(sp_pres,raster1,1,background=0)
  
  n_mod <- summary(as.factor(points_mod_r@data@values))["1"]
  
  # Create a vector with all the species:
  col_n <- c(col_n,n_mod)
  
}  

## TSS, Sensitivity, Specificity, Cut-off (for binary transformation)
# empty vectors to store information for all species
col_tss <- c()
col_sens <- c()
col_spec <- c()
col_cut <- c()

for (i in 1:length(list_sp)) {
  
  setwd("D:/Postdoc_Uppsala")
  # Species name:
  sp <- list_sp[i]
  sp2 <- sub(" ","_",sp)
  sp3 <- sub("_",".",sp2)
  
  # File with all the info. per species:
  file_eval <- paste("D:/Postdoc_Uppsala/", name_model_folder,"/",sp3,"/",sp3,"_eval_ensemble.csv",sep="")
  model_eval <- read.csv2(file_eval)
  # selecting the parameters that I want
  tss_mod <- round(as.numeric(as.character(model_eval[2,5])),2)
  sens_mod <- round(as.numeric(as.character(model_eval[2,7])),2)
  spec_mod <- round(as.numeric(as.character(model_eval[2,8])),2)
  cut_mod <- round(as.numeric(as.character(model_eval[2,6]))/1000,2)
  # Create vectors with all the species:
  col_tss <- c(col_tss,tss_mod)
  col_sens <- c(col_sens,sens_mod)
  col_spec <- c(col_spec,spec_mod)
  col_cut <- c(col_cut, cut_mod)
  
}

## Mean CV within each ensemble model
sp_mean_cv <- c()

for (i in 1:length(list_sp)) {
  sp <- list_sp[i]
  sp2 <- sub(" ","_",sp)
  sp3 <- sub("_",".",sp2)
  
  file_name <- paste(name_model_folder,"/",sp3,"/proj_", model_class,"/individual_projections/",sp3,"_EMcvByTSS_mergedAlgo_mergedRun_mergedData.grd",sep="")
  raster_sp <- raster(file_name)
  
  r.min <- cellStats(raster_sp, "min")
  r.max <- cellStats(raster_sp, "max")
  
  raster_sp_scale <- (raster_sp - r.min) / (r.max - r.min) 
  raster_sp_scale <- raster::mask(raster_sp_scale,raster_to_use)
  raster_list[[i]] <- raster_sp_scale
  
  # Calculate the mean CV for each species:
  mean_sp <- cellStats(raster_sp_scale, stat='mean', na.rm=TRUE)
  sp_mean_cv <- c(sp_mean_cv,mean_sp)
}

## Final table with info for all: 
table_a <- as.data.frame(cbind(as.character(list_sp2),col_n,col_tss,col_sens,col_spec,col_cut,sp_mean_cv))
names(table_a) <- c("Species","n","TSS","Sensitivity","Specificity","Cut-off binary","Mean CV")

heading_mammals <- as.data.frame(cbind("Mammals","","","","","",""))
names(heading_mammals) <- c("Species","n","TSS","Sensitivity","Specificity","Cut-off binary","Mean CV")

heading_birds <- as.data.frame(cbind("Birds","","","","","",""))
names(heading_birds) <- c("Species","n","TSS","Sensitivity","Specificity","Cut-off binary","Mean CV")

heading_amph <- as.data.frame(cbind("Amphibians","","","","","",""))
names(heading_amph) <- c("Species","n","TSS","Sensitivity","Specificity","Cut-off binary","Mean CV")

heading_reptiles <- as.data.frame(cbind("Reptiles","","","","","",""))
names(heading_reptiles) <- c("Species","n","TSS","Sensitivity","Specificity","Cut-off binary","Mean CV")

table_a <- rbind(heading_mammals,table_a[1:9,],heading_birds,table_a[10:13,],heading_amph,table_a[14,],heading_reptiles,table_a[15,])

# Write in the directory:
dir_table <- paste("0Tables/",name_model_folder,"_accuracy.txt",sep="")
write.table(table_a,dir_table,row.names = FALSE)


# TABLE: RESAMPLING vs. AVERAGING ############################################################
setwd("D:/Postdoc_Uppsala")
# Loading resampled Spearman correlation coefficients
clim_r <- read.table("CHELSA_clim/resampled_30k/corr_clim_r.txt")
topo_r <- read.table("GTOPO30/corr_elev_r.txt")
acc_r <- read.table("Accesibility/corr_acc_r.txt")

# Loading resampled Spearman correlation coefficients
clim_a <- read.table("CHELSA_clim/resampled_30k/corr_clim_a.txt")   
topo_a <- read.table("GTOPO30/corr_elev_a.txt")
acc_a <- read.table("Accesibility/corr_acc_a.txt")

# Build table:
heading <- as.data.frame(cbind("Predictor","Resampled","Average"))
names(heading) <- c("","Spearman","")
subh_1 <- as.data.frame(cbind("Climatic","",""))
names(subh_1) <- c("","Spearman","")
subh_2 <- as.data.frame(cbind("Topography","",""))
names(subh_2) <- c("","Spearman","")
subh_3 <- as.data.frame(cbind("Accesibility","",""))
names(subh_3) <- c("","Spearman","")

chunk_clim <- as.data.frame(cbind(c("Annual Mean Temperature","Mean Diurnal Range","Isothermality","Temperature Seasonality","Max Temperature of Warmest Month","Min Temperature of Coldest Month","Temperature Annual Range","Mean Temperature of Wettest Quarter","Mean Temperature of Driest Quarter","Mean Temperature of Warmest Quarter","Mean Temperature of Coldest Quarter","Annual Precipitation","Precipitation of Wettest Month","Precipitation of Direst Month","Precipitation Seasonality","Precipitation of Wettest Quarter","Precipitation of Driest Quarter","Precipitation of Warmest Quarter","Precipitation of Coldest Quarter"),clim_r[,2],clim_a[,2]))
chunk_clim$V2 <- as.character(chunk_clim$V2)
chunk_clim$V3 <- as.character(chunk_clim$V3)
names(chunk_clim) <- c("","Spearman","")

chunk_topo <- as.data.frame(cbind("Elevation",topo_r[2],topo_a[,2]))
chunk_topo[,2] <- as.character(chunk_topo[,2])
chunk_topo[,3] <- as.character(chunk_topo[,3])
names(chunk_topo) <- c("","Spearman","")

chunk_acc <- as.data.frame(cbind("Distance to the closest major city (>50,000)",acc_r[2],acc_a[,2]))
chunk_acc[,2] <- as.character(chunk_acc[,2])
chunk_acc[,3] <- as.character(chunk_acc[,3])
names(chunk_acc) <- c("","Spearman","")

table_corr <- rbind(heading,subh_1,chunk_clim,subh_2,chunk_topo,subh_3,chunk_acc)

# Write in the directory:
write.table(table_corr,"D:/Postdoc_Uppsala/0Tables/corr_july19.txt",row.names = FALSE)


# TABLE: DATA COUNT ############################################################

# empty vectors to store information for all species
table_count <- c()

for (i in 1:length(list_sp)) {
  
  # Species to model:
  sp <- list_sp[i]
  sp2 <- sub(" ","_",sp)
  
  # GLOBAL
  # Ocurrences (gross):
  file1 <- paste("Occ_vs_pres/",sp2,"_occ.shp",sep="")
  glob_occ <- shapefile(file1)
  
  # Presences:
  ## Cert dataset
  file2 <- paste("Occ_vs_pres/",sp2,"_presences_cert_30k.shp",sep="")
  glob_cert <- shapefile(file2)
  
  ## Cert+NA dataset
  file3 <- paste("Occ_vs_pres/",sp2,"_presences_30k_cert_na.shp",sep="")
  glob_cert_na <- shapefile(file3)
  
  # EUROPEAN
  # Ocurrences:
  eu_occ <- crop(glob_occ,eu_raster)
  
  # Presences cert:
  eu_cert <- crop(glob_cert,eu_raster)
  
  # Presentes cert_NA:
  eu_cert_na <- crop(glob_cert_na,eu_raster)
  
  # Compiling for each sp
  line_sp <- c(dim(glob_occ)[1],dim(glob_cert)[1],dim(glob_cert_na)[1],dim(eu_occ)[1],dim(eu_cert)[1],dim(eu_cert_na)[1])
  
  # Joining all spp
  table_count <- rbind(table_count,line_sp)
}

## Final table with info for all: 
table_count2 <- as.data.frame(cbind(as.character(list_sp2),table_count[,c(1:6)]))
names(table_count2) <- c("Species","Occurrences","Certain dataset","Certain+NA dataset","Occurrences","Certain dataset","Certain+NA dataset")

heading_up <- as.data.frame(cbind("","Global","","","European","",""))
names(heading_up) <- c("Species","Occurrences","Certain dataset","Certain+NA dataset","Occurrences","Certain dataset","Certain+NA dataset")

heading_up2 <- as.data.frame(cbind("Species","Occurrences","Certain dataset","Certain+NA dataset","Occurrences","Certain dataset","Certain+NA dataset"))
names(heading_up2) <- c("Species","Occurrences","Certain dataset","Certain+NA dataset","Occurrences","Certain dataset","Certain+NA dataset")

heading_mammals <- as.data.frame(cbind("Mammals","","","","","",""))
names(heading_mammals) <- c("Species","Occurrences","Certain dataset","Certain+NA dataset","Occurrences","Certain dataset","Certain+NA dataset")

heading_birds <- as.data.frame(cbind("Birds","","","","","",""))
names(heading_birds) <- c("Species","Occurrences","Certain dataset","Certain+NA dataset","Occurrences","Certain dataset","Certain+NA dataset")

heading_amph <- as.data.frame(cbind("Amphibians","","","","","",""))
names(heading_amph) <- c("Species","Occurrences","Certain dataset","Certain+NA dataset","Occurrences","Certain dataset","Certain+NA dataset")

heading_reptiles <- as.data.frame(cbind("Reptiles","","","","","",""))
names(heading_reptiles) <- c("Species","Occurrences","Certain dataset","Certain+NA dataset","Occurrences","Certain dataset","Certain+NA dataset")

table_count2 <- rbind(heading_up,heading_up2,heading_mammals,table_count2[1:9,],heading_birds,table_count2[10:13,],heading_amph,table_count2[14,],heading_reptiles,table_count2[15,])

# Write in the directory:
write.table(table_count2,"D:/Postdoc_Uppsala/0Tables/count_global_july19.txt",row.names = FALSE)


# TABLE : RANGE FILLING & ALL-MODELS-OVERLAPED RASTER PER-SPECIES ######################
name_model_folder_g <- "GLOBAL2_30k"
name_model_folder_r <- "REGIONAL2_30k"
dataset_type <- "_presences_cert_na_30k"  # Options: _presences_cert_30k / _presences_30k_cert_na ## from GBIF

vector1 <- c()

for (i in 1:length(list_sp)) {
  
  setwd("D:/Postdoc_Uppsala")
  # Species: 
  sp <- list_sp[i]
  sp2 <- sub(" ","_",sp)
  sp3 <- sub(" ",".",sp)
  
  # Binary global:
  file1 <- paste(name_model_folder_g,"/",sp3,"/proj_global/individual_projections/",sp3,"_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.grd",sep="")
  bin_global <- raster(file1)
  
  # Binary regional:
  file2 <- paste(name_model_folder_r,"/",sp3,"/proj_regional/individual_projections/",sp3,"_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.grd",sep="")
  bin_regional <- raster(file2)
  
  # Presences:
  file_name1 <- paste("Occ_vs_pres/",sp2,dataset_type,".shp",sep="")
  sp_pres <- shapefile(file_name1)
  # Subset of Europe:
  sp_pres_eu <- rasterize(sp_pres,eu_raster,2, background=0)
  
  ## Range filling:
  # Presences against regional binary within EU 
  sum_pres_reg <- overlay(bin_regional,sp_pres_eu,fun=sum)
  # Cells that sum 3 (pres =2 and bin_reg =1) / total bin reg=1
  niche_filling_reg <- summary(as.factor(sum_pres_reg@data@values))["3"] / as.numeric(summary(getValues(bin_regional==1))["TRUE"])
  # Adding up all species:
  vector1 <- c(vector1,niche_filling_reg)
  
  ## Sum of global bin, regional bin and presences, to build summary maps:
  # Multiply global *5 to get unique values afterwards
  bin_global <- bin_global*5
  # Mask to EU
  bin_global <- raster::crop(bin_global,eu_raster)
  bin_global <- raster::mask(bin_global,eu_raster)
  # Overlay of the three
  sum_pres_both <- overlay(sum_pres_reg,bin_global,fun=sum)
  # Save raster for each species:
  file_range <- paste("Maps/",sp2,name_model_folder_g,name_model_folder_r,"_Overlay",sep="")
  writeRaster(sum_pres_both,file_range,format="raster",prj=TRUE,overwrite=TRUE)
}

## Final table with info for all: 
table_range <- as.data.frame(cbind(as.character(list_sp2),vector1))
names(table_range) <- c("Species","Range filling")

heading_mammals <- as.data.frame(cbind("Mammals",""))
names(heading_mammals) <- c("Species","Range filling")

heading_birds <- as.data.frame(cbind("Birds",""))
names(heading_birds) <- c("Species","Range filling")

heading_amphibians <- as.data.frame(cbind("Amphibians",""))
names(heading_amphibians) <- c("Species","Range filling")

heading_reptiles <- as.data.frame(cbind("Reptiles",""))
names(heading_reptiles) <- c("Species","Range filling")

table_range <- rbind(heading_mammals,table_range[c(1:9),],heading_birds,table_range[c(10:13),],heading_amphibians,table_range[14,],heading_reptiles,table_range[15,])

# Write in the directory:
dir_table <- paste("0Tables/",name_model_folder_r,"_range_filling.txt",sep="")
write.table(table_range,dir_table,row.names = FALSE)

name_model_folder <- "REGIONAL5_30k" 
model_class <- "regional"
dataset_type <- "_presences_cert_30k" # presences_cert_na_30k # _presences_cert_30k

# FIGURE: MESS ------------

for (i in 1:length(list_sp)) {
  
  setwd("D:/Postdoc_Uppsala")
  # Species name:
  sp <- list_sp[i]
  sp2 <- sub(" ","_",sp)
  sp3 <- sub("_",".",sp2)
  
  # Presences: 
  file_name1 <- paste("Occ_vs_pres/",sp2,dataset_type,".shp",sep="")
  sp_pres <- shapefile(file_name1)
  
  # Extract reference values:
  ref_pres <- raster::extract(pred_eu_sub,sp_pres)
  
  # Calculate mess, respect to present data:
  mess_current <- mess(pred_eu_sub, ref_pres,full=TRUE)  # full=TRUE provides one raster for each variable (n), and the selected minimum value (n+1), which is the so-called "mess"
  names(mess_current) <- c(names(pred_eu_sub), "rmess")
  
  # Save the multivariate environmental similarity surface:
  file_mess <- paste("D:/Postdoc_Uppsala/MESS/Current_conditions_",sp2,dataset_type,".grd",sep="")
  writeRaster(mess_current[[19]],file_mess,format="raster",prj=TRUE,overwrite=TRUE)
}

my_pattern <- paste("^Current_conditions_.*.",dataset_type,".grd$",sep="")

mess_stack <- stack(list.files(path = "MESS/", pattern = my_pattern, full.names=TRUE))

names(mess_stack) <- c(as.character(list_sp))

mess_table <- as.data.frame(mess_stack,xy=TRUE)

# Count the number of negative values per pixel: 
count_neg <- function(x) {sum(x < 0, na.rm=TRUE)}

mess_table$mess_neg <- apply(mess_table[,c(3:17)],1,count_neg)

# Count the number of species present in each grid-cell to calculate the proportion instead of the absolute count of species with negative values
count <- function(x) {sum(!is.na(x))}

mess_table$tot_spp <- apply(mess_table[,c(3:17)],1,count)

mess_table$prop_mess <-  mess_table$mess_neg / mess_table$tot_spp

# Create a raster with coordinates
count_neg_r <- rasterFromXYZ(mess_table[,c("x","y","mess_neg")])
count_neg_r <- raster::mask(count_neg_r,eu_raster)

my_raster <- paste("MESS/Count_neg_mess",dataset_type,sep="")
writeRaster(count_neg_r,my_raster,format="raster",prj=TRUE,overwrite=TRUE)

# Create a raster with coordinates
prop_mess <- rasterFromXYZ(mess_table[,c("x","y","prop_mess")])
prop_mess <- raster::mask(prop_mess,eu_raster)

my_raster <- paste("MESS/Prop_neg_mess",dataset_type,sep="")
writeRaster(prop_mess,my_raster,format="raster",prj=TRUE,overwrite=TRUE)

# FIGURE: IGNORANCE MAPS ----------

# Europe ----
for (i in 1:length(list_refs)) {
  
  setwd("D:/Postdoc_Uppsala")
  # Family name:
  fam <- list_refs[i]
  
  # Ni:
  name_occ <- paste("GBIF_all/",fam,".csv",sep="")
  spp_occ <- fread(name_occ, select=c("family","species","decimalLatitude","decimalLongitude"),quote="")
  names(spp_occ) <- c("family","species","y","x")
  spp_occ$pres <- 1
  
  # Assign these occurrences to our grid 
  coord_data <- spp_occ[,c("x","y")]
  extract1 <- raster::extract(eu_raster,coord_data,cellnumbers = T)
  coord_raster_ref_sp <- as.data.frame(coordinates(eu_raster)[extract1[,1],])
  df1 <- as.data.frame(coord_raster_ref_sp) # OCURRENCES
  df1$pres <- 1
  
  # Add all cells of Europe to count them as 0:
  eu_table <- as.data.frame(eu_raster,xy=TRUE)
  eu_table <- na.omit(eu_table)
  
  joint_eu <- merge(eu_table[,c("x","y")],df1,by=c("x","y"),all.x=TRUE)
  joint_eu$pres[is.na(joint_eu$pres)] <- 0
  joint_eu <- joint_eu[,c("x","y","pres")]
  
  sum_gridcell <- aggregate(joint_eu$pres,list(joint_eu$x,joint_eu$y),sum,na.rm=TRUE)
  names(sum_gridcell) <- c("x","y","Sum_occ")
  
  ## Ign3: I3 = O0.5 / (Ni + o0.5). setting O0.5 = 1 setting the reference number O = 1 means that one observation is enough to consider that the absence of reports of a target species from any grid cell is 50% due to true absence from the site and 50% due to failure to detect the species.
  sum_gridcell$Ign3 <- 1 /(sum_gridcell$Sum_occ + 1)
  
  # Save this in raster to visualize:
  ign3 <- SpatialPointsDataFrame(coords = cbind(sum_gridcell$x,sum_gridcell$y),proj4string = CRS(projection(eu_raster)), data=as.data.frame(sum_gridcell$Ign3))
  ign_raster <- rasterize(ign3,eu_raster,field="sum_gridcell$Ign3")
  ign_raster <- raster::crop(ign_raster,eu_raster)
  ign_raster <- raster::mask(ign_raster,eu_raster)
  # Save
  name_to_save4 <- paste("Maps/Fam_Ign3_eu_index_",fam,".grd",sep="")
  writeRaster(ign_raster,name_to_save4,format="raster",prj=TRUE,overwrite=TRUE)
}

# Calculate the average over all species ---
mess_stack <- stack(list.files(path = "Maps/", pattern = "^Fam.*.grd$", full.names=TRUE))
mean_ign <- calc(mess_stack, fun = mean, na.rm = T)
mean_ign <- raster::crop(mean_ign,eu_raster)
mean_ign <- raster::mask(mean_ign,eu_raster)

writeRaster(mean_ign,"Maps/Mean_fam_ign3.grd",format="raster",prj=TRUE,overwrite=TRUE)


# SUMMARY TABLE OF MESS, IGNORANCE VALUES AND CV ####################

## Zones map
zones_map <- raster("Maps/ZonesREGIONAL5_30kGLOBAL2_30k.grd")

## Mean CV map
cv_map <- raster("Maps/Mean_cv_REGIONAL5_30k_CV.grd")

## Ignorance map
ign_map <- raster("Maps/Mean_fam_ign3.grd")

## MESS
# Prop of spp per grid cell
mess_map <- raster("MESS/Prop_neg_mess_presences_cert_30k.grd") # Prop of negative spp

# All species 
mess_stack <- stack(list.files(path = "MESS/", pattern = my_pattern, full.names=TRUE))
names(mess_stack) <- c(as.character(list_sp))

all_mess_df <- as.data.frame(mess_stack, xy=TRUE)
all_mess_df$average <- rowMeans(all_mess_df[,c(3:17)])

mess_av <- rasterFromXYZ(all_mess_df[,c("x","y","average")])

# Stack and convert to table
stack_all <- stack(zones_map,ign_map,mess_map,cv_map)
stack_all_df <- as.data.frame(stack_all,xy=TRUE)
names(stack_all_df) <- c("x","y","zone","ign","mess","cv")

new_table <- aggregate(ign ~ zone, stack_all_df, mean)
new_table$ign_min <- aggregate(ign ~ zone, stack_all_df, min)
new_table$ign_max <- aggregate(ign ~ zone, stack_all_df, max)
new_table$mess <- aggregate(mess ~ zone, stack_all_df, mean)
new_table$mess_min <- aggregate(mess ~ zone, stack_all_df, min)
new_table$mess_max <- aggregate(mess ~ zone, stack_all_df, max)
new_table$cv <- aggregate(cv ~ zone, stack_all_df, mean)
new_table$cv_min <- aggregate(cv ~ zone, stack_all_df, min)
new_table$cv_max <- aggregate(cv ~ zone, stack_all_df, max)

write.table(new_table,"Maps/table_ign_mess.csv",row.names = FALSE)


## REPETIR LO DEL MESS CON EL NUMERO DE ESPECIES, PQ LA MEDIA SALE NEGATIVA Y NO ME DICE NADA (SI USO LA MEDIA DEBERIA ESTANDARIZAR O ALGO). PENSAR SI HACER MAPA O TABLA #####
zone_ign <- sp::merge(zone_shp,new_table,by="zone",incomparables=NA)
