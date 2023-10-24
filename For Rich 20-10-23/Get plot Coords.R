## Get the coordinates for our density plots


## Read in the centers of each density plot
Plot_centers <- st_read(dsn = "SpatialData", 
                        layer = "centroid_survey_point")

## read in the density plot surveys
Density_plots <- read.csv("TabularData/Denisty plot data.csv")
Density_plots <- Density_plots[1:(nrow(Density_plots)-1),] ## remove last row as it just summaries

## get the lat longs and add them on
latlong <- st_transform(Plot_centers, crs = 4326) %>% st_coordinates()
Plot_centers$Lat <- latlong[,2]
Plot_centers$Long <- latlong[,1]

## change this column name for the join
data.table::setnames(Density_plots, old = "Grid_Centroid", new = "plot_id")

## DO the join to keep just the centorids that we surveyed
Plot_centers <- subset(Plot_centers, select=c("plot_id", "xcoord", "ycoord", "Lat", "Long"))
Density_centers_reproj <- inner_join(Plot_centers, Density_plots, by = "plot_id")
stopifnot(nrow(Density_centers_reproj)==nrow(Density_plots))


## read out the data as a csv and a shapefile
Density_centers_reproj <- Density_centers_reproj %>% dplyr::select(-c(X, Survey_Team))
Density_centers_txt <- Density_centers_reproj %>% st_drop_geometry()
write_csv(Density_centers_txt, "For Rich 20-10-23/Shearwater_Density_Plots_Rum21.csv")
st_write(Density_centers_reproj, "For Rich 20-10-23/Shearwater_Density_Plots_Rum21.shp")





