## Luke Ozsanlav-Harris

## November 2021
## Uploaded to GitHub March 2023

## Download remote sensing data for Manx Shearwater colony on the Isle of Rum
## This will be Sentinel data (For NDVI) and a Digital elevation model (DEM)
## Correct spatial data for terrain and derive slope
## Make model that can be used to predict burrow denisty acroos the whole cony from the denisty plots


## Packages required
library(tidyverse)
library(raster)
library(sf)
library(rgdal)
library(rsq)
library(glmmTMB)
library(performance)
library(DHARMa)
library(sjPlot)
library(stargazer)

## 1. Download the DEM data using elevatr
## 2. Import Sentinel 2 data and create spectral indices
## 3. Align the Sentinel data and the DEM
## 4. Add slope/aspect and surface area to Raster stack
## 5. Create a Raster that contains colony identity
## 6. Format colony grid centers and extract the ones with density plot counts
## 7. Interpolate Raster stack values onto the centers of the density plots
## 8. 'Average' Raster values over the whole density plots
## 9. Model burrow Density against spatial covariates
## 10. Predict density over the colony extent
## 11. Calculate number of burrows and surface for each pixel within the colony extent
## SEND OUTPUT TO RICH TO COMNVERT FROM NUMBER OF BURROWS TO NUMBER OF 'ACTIVE' BURROWS

## **NOT NEEDED FOR FINAL ANALYSIS**
## 12. Calculate surface area of low, medium and high density
## 13. Plot relationship between burrow density and spatial covariates
## 14. Label Barkeval south slope pixels in predicted density list






##
####---- 1: Download the DEM data using elevatr ----####
##

## read in the colony bounding box
Bound_box <- st_read(dsn = "SpatialData", 
                     layer = "Rum colony bounding box")
st_crs(Bound_box) ## check the crs of the layer
prj_dd <- "EPSG:4326"
BB_Box <- as.data.frame(Bound_box[[2]][[1]][[1]]) ## create a bounding box of lat longs as the shapefile didn't work


## retrieve the raster elevation data
## units of the raster are lat longs
## MORE INFO ON DATA SET USED: SRTMGL1 (SRTM GL1 30m) https://portal.opentopography.org/apidocs/#/Public/getGlobalDem 
DEM <- elevatr::get_elev_raster(locations = BB_Box, prj = prj_dd,
                                src = "gl1", verbose = TRUE)
DEM <-raster("SpatialData/Original_Rum_DEM.grd")
DEM # get the raster info
plot(DEM) # plot the raster





##
####---- 2: Import Sentinel 2 data and create spectral indices ----####

## Read in the data down loaded from Landviewer
Sent_B02 <-raster("NDVI data/S2A_tile_20210530_B04-B03-B02-B8A-B80/S2A_tile_20210530_29VPD_0_R10B02.tif")
Sent_B03 <-raster("NDVI data/S2A_tile_20210530_B04-B03-B02-B8A-B80/S2A_tile_20210530_29VPD_0_R10B03.tif")
Sent_B04 <-raster("NDVI data/S2A_tile_20210530_B04-B03-B02-B8A-B80/S2A_tile_20210530_29VPD_0_R10B04.tif")
Sent_B08 <-raster("NDVI data/S2A_tile_20210530_B04-B03-B02-B8A-B80/S2A_tile_20210530_29VPD_0_R20B08.tif")
Sent_B8A <-raster("NDVI data/S2A_tile_20210530_B04-B03-B02-B8A-B80/S2A_tile_20210530_29VPD_0_R20B8A.tif")
# plot(Sent_B02)

## Create 20m x 20m versions of 10m rasters
Sent_R20_B02 <- projectRaster(from = Sent_B02, to =Sent_B8A)
Sent_R20_B03 <- projectRaster(from = Sent_B03, to =Sent_B8A)
Sent_R20_B04 <- projectRaster(from = Sent_B04, to =Sent_B8A)
Sent_R20_B08 <- projectRaster(from = Sent_B08, to =Sent_B8A)

## Calculate spectral indices from Sentinel ##

## Calculate NDVI from the two layers
Sent_NDVI <- (Sent_B8A-Sent_R20_B04)/(Sent_R20_B04+Sent_B8A)
Sent_NDVI; plot(Sent_NDVI)
names(Sent_NDVI) <- "NDVI"


## Calculate Green Coverage Index (GCI):
## GCI = (NIR) / (Green) - 1
Sent_GCI <- Sent_B08/(Sent_B03-1)
Sent_GCI; plot(Sent_GCI)
names(Sent_GCI) <- "GCI"


## Calculate Adjusted transformed soil-adjusted VI (SAVI)
## Found formula at this link: https://www.indexdatabase.de/db/si-single.php?sensor_id=96&rsindex_id=209
Sent_SAVI <- 1.22*(((Sent_B8A-1.22)*(Sent_R20_B04-0.03))/((1.22*Sent_B8A)+(Sent_R20_B04)-(1.22*0.03)+(1*(1.22^2))))
Sent_SAVI; plot(Sent_SAVI)
names(Sent_SAVI) <- "SAVI"

## Calculate Enhanced vegetation index (EVI)
## EVI (Sentinel 2) = 2.5 * ((B8 - B4) / (B8 + 6 * B4 - 7.5 * B2 + 1))
## https://www.indexdatabase.de/db/si-single.php?sensor_id=96&rsindex_id=16





##
####---- 3: Align the Sentinel data and the DEM ----####
##

## Transform the DEM CRS to that of Sentinel
DEM_reporj <- projectRaster(DEM, crs = crs(Sent_NDVI))
DEM_reporj
names(DEM_reporj) <- "Elevation"

## Change the resolution and extent of the Sentinel data to that of the DEM
values(Sent_SAVI)[values(Sent_SAVI)>2000] <- NA
Sent_stack <- stack(Sent_NDVI, Sent_SAVI)
Sent_stack_tr <- projectRaster(from = Sent_stack, to =DEM_reporj)
plot(Sent_stack_tr[[2]])

## write out the transformed remote sensing data
#writeRaster(Sent_stack_tr[[1]], filename = "SpatialData/For map 08-11/Sentinel_NDVI")
#writeRaster(Sent_stack_tr[[2]], filename = "SpatialData/For map 08-11/Sentinel_SAVI")
#writeRaster(DEM_reporj, filename = "Elevation")





##
####---- 4: Add slope/aspect and surface area to Raster stack ----####
##

## Calculate and save the slope of the raster
DEM_slope <- terrain(DEM_reporj, opt="slope", unit= 'degrees', neighbors=8)
# plot(DEM_slope)
#writeRaster(DEM_slope, filename = "Slope")

## Calculate and save the aspect of the raster
DEM_aspect <- terrain(DEM_reporj, opt="aspect", unit= 'degrees', neighbors=8)
# plot(DEM_aspect)
#writeRaster(DEM_aspect, filename = "Aspect")


## calculate the total surface area of each pixel
## have to convert to a SpatialPixelsDataFrame first to calculate surface area
DEM_SpatPix <- as(DEM_reporj, "SpatialPixelsDataFrame")
DEM_SurfaceArea <- surfaceArea(DEM_SpatPix, byCell = T)
DEM_SurfaceArea2 <- as(DEM_SurfaceArea, "RasterLayer")
plot(DEM_SurfaceArea2)
names(DEM_SurfaceArea2) <- "SurfaceArea"
#writeRaster(DEM_SurfaceArea2, filename = "SurfaceArea")


## Create a Raster stack of all the spatial predictor vairables
Rum_stack <- raster::stack(list(DEM_reporj, DEM_slope, DEM_aspect, Sent_stack_tr))
Rum_stack@layers





##
####---- 5: Create a Raster that contains colony identity ----####
##

## Read in the shape files with the colony shapefiles in
Colony <- st_read(dsn = "SpatialData/New colony shapefiles", 
                          layer = "temp2")

Additions <- st_read(dsn = "SpatialData/New colony shapefiles", 
                     layer = "Extra_shapes")

## make sure the shape files have the same crs
Add_tr <- st_transform(Additions, crs = st_crs(Colony))

## now merge the two shape files and remove internal boundaries
Colony_outline <- st_union(x=Colony, y=Add_tr)
Colony_outline$id <- 1
Colony_comb <- Colony_outline %>% st_buffer (1.5) %>% st_union() %>% st_sf()
Colony_comb; plot(Colony_comb)

## transform the CRS of the shape file to that of the Raster data
Colony_comb_tr <- st_transform(Colony_comb, crs = st_crs(DEM_reporj)) 

## create an empty raster with the same specification of the DEM data
Colony_extent_raster <- DEM_reporj
values(Colony_extent_raster) = 0

## extract the rater cells covered by each of the colony
colony_index <- extract(Colony_extent_raster, Colony_comb_tr, cellnumbers=TRUE)

## Assign pixels within the colony extent with a 1
values(Colony_extent_raster)[c(colony_index[[1]][,1])] = 1
plot(Colony_extent_raster)

## change the raster layer name
names(Colony_extent_raster) <- "colony_ext"

## write out the shape file of the extent of active burrows
plot(Colony_comb_tr)
#st_write(Colony_comb_tr, dsn = "SpatialData", 
#         layer = "Active_burrow_extent", driver = "ESRI Shapefile")





## Now Create colony shape file with colony ID ###

## read in the shape file of the colony boundaries
Colony_outline <- st_read(dsn = "SpatialData", 
                          layer = "colony_boundaries")

## transform the CRS of the shape file to that of the Raster data
Colony_outline_trans <- st_transform(Colony_outline, crs = st_crs(DEM_reporj))

## create an empty raster with the same specification of the DEM data
Colony_raster <- DEM_reporj
values(Colony_raster) = 0

## extract the rater cells covered by each of the three colonies
colony_index <- extract(Colony_raster, Colony_outline_trans, cellnumbers=TRUE)

## Assign colony IDs to the raster
## 1= Barkeval, 2= Askival and Hallival and 3 = Trollobal
values(Colony_raster)[c(colony_index[[1]][,1])] = 1
values(Colony_raster)[c(colony_index[[2]][,1])] = 2
values(Colony_raster)[c(colony_index[[3]][,1])] = 3
plot(Colony_raster)

## change the name of the layer
names(Colony_raster) <- "ColonyID"



## Save the colony data
#writeRaster(Colony_raster, filename = "SpatialData/For map 08-11ColonyID")





##
####---- 6. Format colony grid centers and extract the ones with density plot counts ----####
##

## Read in the centers of each density plot
Plot_centers <- st_read(dsn = "SpatialData", 
                        layer = "centroid_survey_point")

## read in the density plot surveys
Density_plots <- read.csv("TabularData/Denisty plot data.csv")
Density_plots <- Density_plots[1:(nrow(Density_plots)-1),] ## remove last row as it just summaries



## Change the CRS of the plot centers so that it matches with all of the Raster data
Plot_centers_reproj <- st_transform(Plot_centers, crs = st_crs(DEM_reporj))

## Extract just the grid numbers from the density plot data
## Streamline the Density plot data
Density_plots_slim <- subset(Density_plots, select=c("Grid_Centroid", "Number_Burrows", "Area_Surveyed"))
data.table::setnames(Density_plots_slim, old = "Grid_Centroid", new = "plot_id")

## Join the Density plot data and the grid data to keep just the centers for the Density plot 
Density_centers_reproj <- inner_join(Plot_centers_reproj, Density_plots_slim, by = "plot_id")
stopifnot(nrow(Density_centers_reproj)==nrow(Density_plots_slim))



## Create a column for the density plot radius
Density_centers_reproj$Radius <- Density_centers_reproj$Area_Surveyed/31.4

## Extract just the columns with data in
Density_centers_reproj <- Density_centers_reproj[colSums(!is.na(Density_centers_reproj)) > 0]





##
####---- 7. Interpolate Raster stack values onto the centers of the density plots ----####
##

## Simplify name of main data set
Dens_plots <- Density_centers_reproj

## Use nearest neighbor interpolation for each center
nn <- raster::extract(x= Rum_stack, y= Dens_plots, method = "simple")
Dens_plots <- cbind(Dens_plots, nn)
data.table::setnames(Dens_plots, old =c("Elevation", "slope", "aspect", "NDVI", "SAVI"), new = c("NN_Elevation", "NN_slope", "NN_aspect", "NN_NDVI", "NN_SAVI"))

## Use bi linear interpolation for each center
bi <- raster::extract(x= Rum_stack, y= Dens_plots, method = 'bilinear')
Dens_plots <- cbind(Dens_plots, bi)
data.table::setnames(Dens_plots, old =c("Elevation", "slope","aspect", "NDVI", "SAVI"), new = c("Bi_Elevation", "Bi_slope", "Bi_aspect", "Bi_NDVI", "Bi_SAVI"))

## Assign colony ID to each of the 
col <- raster::extract(x= Colony_raster, y= Dens_plots, method = 'simple')
Dens_plots <- cbind(Dens_plots, col)




##
####---- 8. 'Average' Raster values over the whole density plots ----####
##

## Calculate a new radius based off the slope for the plot center
## Use trigonometry to correct for the slope
Dens_plots$cor_Radius <- Dens_plots$Radius*cos((Dens_plots$NN_slope*(pi/180)))

## Add the corrected radius as a buffer to these points
Dens_Buf <- st_buffer(Dens_plots, dist = Dens_plots$cor_Radius)

## Average the values in the raster stack over the whole Density plot buffers
## This works by finding the pixels the buffer lies in or overlaps with
## It then calculates weights for each of these pixels based off the % of the buffer that overlaps with the pixel
## Then the function argument then calculates the mean for elevation, slope, NDVI using the %s as weights
Ave1 <- extract(x= Rum_stack, y= Dens_Buf, weights= TRUE, normalizeWeights= TRUE, fun= mean)
Dens_plots <- cbind(Dens_plots, Ave1)
Dens_plots$aspect <- NULL # remove aspect as you don't want to average it
data.table::setnames(Dens_plots, old =c("Elevation", "slope", "NDVI", "SAVI"), new = c("WeightedBuff_Elevation", "WeightedBuff_slope", "WeightedBuff_NDVI", "WeightedBuff_SAVI"))


## write out this data set now
#write_csv(Dens_plots, file = "Spatial_covariates_for_density_plots.csv")

## write out the Density plot buffers
#st_write(Dens_Buf, dsn = "SpatialData", 
#         layer = "Density_plots_Buff", driver = "ESRI Shapefile")





##
####---- 9. Model burrow Density against spatial covariates ----####
##

## Calculate burrow Density
Dens_plots$Density <- (Dens_plots$Number_Burrows/Dens_plots$Area_Surveyed)

## check correlation between predictors
cor_subset <- subset(as.data.frame(Dens_plots), select= c("WeightedBuff_slope", "WeightedBuff_Elevation", "WeightedBuff_NDVI"))
cor(cor_subset)

## Scale the explanatory
Dens_plots2 <- as.data.frame(Dens_plots)
Dens_plots2[, c("WeightedBuff_slope", "WeightedBuff_Elevation", "WeightedBuff_NDVI")] <- scale(
  Dens_plots2[, c("WeightedBuff_slope", "WeightedBuff_Elevation", "WeightedBuff_NDVI")])


## Model the response as a beta family using glmmTMB
Dens_plots$col <- as.factor(Dens_plots$col)
Dens_plots$WeightedBuff_Elevation_sc <- scale(Dens_plots$WeightedBuff_Elevation)
Dens_plots$WeightedBuff_NDVI_sc <- scale(Dens_plots$WeightedBuff_NDVI)
Dens_plots$WeightedBuff_slope_sc <- scale(Dens_plots$WeightedBuff_slope)

Model2 <- glmmTMB(Density~ col + WeightedBuff_slope_sc + WeightedBuff_Elevation_sc + WeightedBuff_NDVI_sc + col*WeightedBuff_NDVI_sc, 
              data = Dens_plots, 
              family = beta_family)
summary(Model2)

## test model assumptions using DHARMa, simulate randomized quantile residuals
simulationOutput <- simulateResiduals(fittedModel = Model2, plot = T) # not violated

## test for outliers, dispersion and conformity to expected distribution
testResiduals(simulationOutput)


## Run model selection
## change default "na.omit" to prevent models being fitted to different datasets
library(MuMIn)
options(na.action = "na.fail") 

## create all candidate models using dredge, specify any dependencies (trace shows progress bar)
ms2 <- MuMIn::dredge(Model2, trace = 2)
ms2_sub <- subset(ms2, !MuMIn::nested(.), recalc.weights=T)
ms2_sub6 <- subset(ms2_sub, delta<=6, recalc.weights=T)




##
####---- 9.2 Plot MuMin model selection table ----####
##

## create a table from the model selection
## run models from top set
Top1 <- glmmTMB(Density~ WeightedBuff_slope_sc + WeightedBuff_Elevation_sc + WeightedBuff_NDVI_sc, 
                  data = Dens_plots, 
                  family = beta_family)
Top2 <- glmmTMB(Density~ WeightedBuff_Elevation_sc + WeightedBuff_NDVI_sc, 
                  data = Dens_plots, 
                  family = beta_family)
## extract CIs
CIs1 <- as.data.frame(round(confint(Top1), digits = 2))
CIs1$report <- paste0(CIs1$Estimate, " [", CIs1$`2.5 %`, " ", CIs1$`97.5 %`, "]")
CIs2 <- as.data.frame(round(confint(Top2), digits = 2))
CIs2$report2 <- paste0(CIs2$Estimate, " [", CIs2$`2.5 %`, " ", CIs2$`97.5 %`, "]")

## re-organised the data so that i can add the CIs to the model selection table
CIs1$names <- rownames(CIs1)
CIs2$names <- rownames(CIs2)
CIs3 <- left_join(CIs1, CIs2[,4:5], by = "names")
rownames(CIs3) <- CIs3$names
CIs3 <- t(CIs3)
CIs3 <- as.data.frame(CIs3[c(4,6),])

## change column names from the model selection table
ms2_sub6$`disp((Int))` <- NULL
colnames(ms2_sub6)
data.table::setnames(ms2_sub6, old = c("cond((Int))", "cond(col)", "cond(WeightedBuff_Elevation_sc)","cond(WeightedBuff_NDVI_sc)",     
                           "cond(WeightedBuff_slope_sc)", "cond(col:WeightedBuff_NDVI_sc)"), 
                   new = c("Intercept", "Colony", "Elevation", "NDVI", "Slope", "Colony*NDVI"))
## add CIs to the modle selection table
ms2_sub6$Intercept <- CIs3$`(Intercept)`
ms2_sub6$NDVI <- CIs3$WeightedBuff_NDVI_sc
ms2_sub6$Elevation <- CIs3$WeightedBuff_Elevation_sc
ms2_sub6$Slope <- CIs3$WeightedBuff_slope_sc
ms2_sub6$weight <- NULL

## create and save the model selection table
#tab_df(ms2_sub6, digits = 2, alternate.rows = T, file = "Model_selection_table_files/Model_selection_table.doc")





##
####---- 10. Predict density over the colony extent ----####
##

## Create a raster stack with all the rasters i need in
Rum_stack2 <- raster::stack(list(DEM_reporj, DEM_slope, DEM_aspect, Sent_stack_tr[[1]], Colony_raster, DEM_SurfaceArea2, Colony_extent_raster))

## Now turn it into a data frame so I can filter out just the colony extent
Col_pixels <- as.data.frame(Rum_stack2)
Col_pixels2 <- filter(Col_pixels, colony_ext > 0)

## Change the names of the variables so that we can predict from a model
data.table::setnames(Col_pixels2, old =c("Elevation", "slope", "NDVI"), new = c("WeightedBuff_Elevation", "WeightedBuff_slope", "WeightedBuff_NDVI"))

## Run the model and then predict density
Pred_model <- glmmTMB(Density~ WeightedBuff_slope + WeightedBuff_Elevation + WeightedBuff_NDVI, 
                  data = Dens_plots, 
                  family = beta_family)
summary(Pred_model)
Pred_dens <- predict(Pred_model, newdata = Col_pixels2)
Pred_dens <- exp(Pred_dens)
Col_pixels2 <- cbind(Col_pixels2, Pred_dens)





##
####---- 11. Calculate number of burrows and surface for each pixel within the colony extent ----####
##


## Create raster of predicted density
Colony_PredDens <- Colony_extent_raster
Colony_PredDens[values(Colony_PredDens)>0] <- Col_pixels2$Pred_dens
Colony_PredDens[values(Colony_PredDens)==0] <- NA
plot(Colony_PredDens)
names(Colony_PredDens) <- "Predicted_Density"
#writeRaster(Colony_PredDens, file = "SpatialData/New colony shapefiles/PredDensityRaster", overwrite=TRUE)

## extract the predicted densities from the raster
colnames(Col_pixels2)
PredDens_values <- subset(Col_pixels2, select = c("SurfaceArea", "Pred_dens"))


## Create a raster of predicted number of burrows in the colony extent
## Calculate the number of burrows
Col_pixels2$NoBurrows <- Col_pixels2$Pred_dens*Col_pixels2$SurfaceArea

## Assign these values to the raster of the colony
Colony_NoBurrows <- Colony_extent_raster
Colony_NoBurrows[values(Colony_NoBurrows)>0] <- Col_pixels2$NoBurrows
Colony_NoBurrows[values(Colony_NoBurrows)==0] <- NA
plot(Colony_NoBurrows)
names(Colony_NoBurrows) <- "No_burrows"
#writeRaster(Colony_NoBurrows, file = "SpatialData/New colony shapefiles/DensityRaster")

## return the total number of burrows
sum(Col_pixels2$NoBurrows)
sum(Col_pixels2$NoBurrows[Col_pixels2$ColonyID ==1])
sum(Col_pixels2$NoBurrows[Col_pixels2$ColonyID ==2])
sum(Col_pixels2$NoBurrows[Col_pixels2$ColonyID ==3])





##
####---- 12. Calculate surface area of low, medium and high density ----####
##

## Create low, medium and high strata for predicted density and create raster of it
## First calculate the three Quantiles from the density data
Q33 <- quantile(Dens_plots$Density, probs = 1/3)
Q66 <- quantile(Dens_plots$Density, probs = 2/3)

## create raster of predicted density in the colony
plot(Colony_extent_raster)
Colony_dens <- Colony_extent_raster
Colony_dens[values(Colony_dens)>0] <- Col_pixels2$Pred_dens
Colony_dens[values(Colony_dens)==0] <- NA
plot(Colony_dens)

## Now reassign values in the raster for the three strata
Colony_strata <- Colony_dens
Colony_strata[values(Colony_strata)> Q66 & values(Colony_strata)< 2] = 3
Colony_strata[values(Colony_strata)> Q33 & values(Colony_strata)< Q66] = 2
Colony_strata[values(Colony_strata)< Q33] = 1
plot(Colony_strata)
names(Colony_strata) <- "strata"

## Save this raster of strata so i can plot it in QGIS
#writeRaster(Colony_strata, file = "SpatialData/New colony shapefiles/DensityStrata", overwrite=TRUE)

## Calculate the surface area of each strata
Strata_stack <- stack(Colony_strata, DEM_SurfaceArea2, Colony_PredDens, Colony_NoBurrows)
Strata_area <- as.data.frame(Strata_stack)
Strata_summary <- Strata_area %>% 
                  drop_na(strata) %>% 
                  group_by(strata) %>% 
                  summarise(total_area = sum(SurfaceArea), 
                            No_burrows = sum(No_burrows))




##
####---- 13. Plot relationship between burrow density and spatial covariates ----####
##



ModelTop <- glmmTMB(Density~ WeightedBuff_slope + WeightedBuff_Elevation + WeightedBuff_NDVI, 
                  data = Dens_plots, 
                  family = beta_family)

## extract the fitted values plus CIs using the effects package
effectz <- effects::effect(term= "WeightedBuff_NDVI", mod= ModelTop, xlevels= 100)
effectz2 <- as.data.frame(effectz)

ggplot() + 
  geom_point(data=Dens_plots, aes(x=WeightedBuff_NDVI, y=Density, colour = "#D55E00"), size =2) + 
  geom_line(data=effectz2, aes(x= WeightedBuff_NDVI, y = fit), size = 1.25)  +
  geom_ribbon(data = effectz2, aes(x=WeightedBuff_NDVI, ymin = lower, ymax = upper), alpha = 0.5, colour = NA, fill = "grey") + 
  ylab("Density(Burrows/m^2)") + xlab("Normalized vegetation index") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), legend.position = "none",
        axis.text=element_text(size=13, face = "bold"), axis.title=element_text(size=17, face = "bold"), 
        plot.title = element_text(size=14, face="bold"), legend.text=element_text(size=12), legend.title=element_text(size=12),
        panel.grid.minor.x = element_blank(), strip.text.x = element_text(size = 13, face = "bold"))






## extract the fitted values plus CIs using the effects package
effectz <- effects::effect(term= "WeightedBuff_slope", mod= ModelTop, xlevels= 100)
effectz3 <- as.data.frame(effectz)

ggplot() + 
  geom_point(data=Dens_plots, aes(x=WeightedBuff_slope, y=Density, colour = "#D55E00"), size =1.5) + 
  geom_line(data=effectz3, aes(x= WeightedBuff_slope, y = fit), size = 1.25)  +
  geom_ribbon(data = effectz3, aes(x=WeightedBuff_slope, ymin = lower, ymax = upper), alpha = 0.5, colour = NA, fill = "grey") + 
  ylab("Density(Burrows/m^2)") + xlab("Slope (degrees)") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), legend.position = "none",
        axis.text=element_text(size=13, face = "bold"), axis.title=element_text(size=17, face = "bold"), 
        plot.title = element_text(size=14, face="bold"), legend.text=element_text(size=12), legend.title=element_text(size=12),
        panel.grid.minor.x = element_blank(), strip.text.x = element_text(size = 13, face = "bold"))








## extract the fitted values plus CIs using the effects package
effectz <- effects::effect(term= "WeightedBuff_Elevation", mod= ModelTop, xlevels= 100)
effectz4 <- as.data.frame(effectz)

ggplot() + 
  geom_point(data=Dens_plots, aes(x=WeightedBuff_Elevation, y=Density, colour = "#D55E00"), size =1.5) + 
  geom_line(data=effectz4, aes(x= WeightedBuff_Elevation, y = fit), size = 1.25)  +
  geom_ribbon(data = effectz4, aes(x=WeightedBuff_Elevation, ymin = lower, ymax = upper), alpha = 0.5, colour = NA, fill = "grey") + 
  ylab("Density(Burrows/m^2)") + xlab("Elevation (m)") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), legend.position = "none", 
        axis.text=element_text(size=13, face = "bold"), axis.title=element_text(size=17, face = "bold"), 
        plot.title = element_text(size=14, face="bold"), legend.text=element_text(size=12), legend.title=element_text(size=12),
        panel.grid.minor.x = element_blank(), strip.text.x = element_text(size = 13, face = "bold"))







##
####---- 14. Label Barkeval south slope pixels in predicted density list ----####
##


## Read in a shapefile of the barkeval south slope section 
Bark_slope <- st_read(dsn = "SpatialData", 
                     layer = "Bark_south_slope")
plot(Bark_slope)


## Create a raster with Barkeval south slope pixels labelled as a 1 and all others as a 0
## transform the CRS of the shape file to that of the Raster data
Bark_slope_tr <- st_transform(Bark_slope, crs = st_crs(DEM_reporj))


## create an empty raster with the same specification of the DEM data
South_slope_raster <- DEM_reporj
values(South_slope_raster) = 0

## extract the rater cells covered by each of the colony
South_slope_index <- extract(South_slope_raster, Bark_slope_tr, cellnumbers=TRUE)

## Assign pixels within the colony extent with a 1
values(South_slope_raster)[c(South_slope_index[[1]][,1])] = 1
plot(South_slope_raster)

## change the raster layer name
names(South_slope_raster) <- "South_slope"



## Create a raster stack with all the raters I need in
Rum_stack3 <- raster::stack(list(DEM_reporj, DEM_slope, DEM_aspect, Sent_stack_tr[[1]], Colony_raster, DEM_SurfaceArea2, Colony_extent_raster, South_slope_raster))

## Now turn it into a data frame so I can filter out just the colony extent
All_pixels <- as.data.frame(Rum_stack3)
All_pixels2 <- filter(All_pixels, colony_ext > 0)

## Change the names of the variables so that we can predict from a model
data.table::setnames(All_pixels2, old =c("Elevation", "slope", "NDVI", "ColonyID"), new = c("WeightedBuff_Elevation", "WeightedBuff_slope", "WeightedBuff_NDVI", "col"))

## Run the model and then predict density
Pred_model <- glmmTMB(Density~ WeightedBuff_slope + WeightedBuff_Elevation + WeightedBuff_NDVI, 
                      data = Dens_plots, 
                      family = beta_family)
summary(Pred_model)
Pred_dens <- predict(Pred_model, newdata = All_pixels2)
Pred_dens <- exp(Pred_dens) # back transform the density estimates

## bind the predicted densities back to the main data set
All_pixels2 <- cbind(All_pixels2, Pred_dens)
summary(All_pixels2)

## write out the data set
PredDens_values2 <- subset(All_pixels2, select = c("SurfaceArea", "Pred_dens", "South_slope"))
#write_csv(PredDens_values2, file = "TabularData/Predicted_densities_with_southslope.csv")

## Check the number of burrows it is predicting in the colony and on that section of the south slope
southslope <- All_pixels2 %>% filter(South_slope == 1)
sum(southslope$Pred_dens*southslope$SurfaceArea)
sum(All_pixels2$Pred_dens*All_pixels2$SurfaceArea)


## Read in the data set with the predicted number of AOBs
Pred_AOB <- read_csv("TabularData/Predicted_AOBs_Rich.csv")

## Bind on the column to denote if a pixel is on the south slope of Barkeval
Pred_AOB$Bark_south <- All_pixels2$South_slope

## Now wirte this file back out 
write_csv(Pred_AOB, file = "TabularData/Predicted_AOBs_with_southslope.csv")

## check the number of AOBs in the south slpoe and in the whole colony
southslope2 <- Pred_AOB %>% filter(Bark_south == 1)
sum(southslope2$Pred_AOBs)
sum(Pred_AOB$Pred_AOBs)
