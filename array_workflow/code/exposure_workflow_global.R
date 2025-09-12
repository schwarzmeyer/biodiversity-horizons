
# This script demonstrates the workflow for estimating species exposure to climate change using refactored functions. 
# The workflow is divided into five modules: preparing climate data, preparing species data, estimating niche limits, computing exposure, and calculating exposure times.

library(dplyr)
library(tidyr)
library(stringr)
library(terra)
library(arrow)
library(glue)
library(here)
library(sf)
library(parallel)

# Load functions
source(here("code/functions/raster_to_tibble.R"))
source(here("code/functions/array_to_tibble.R"))
source(here("code/functions/get_range_cells.R"))
source(here("code/functions/climate_limits.R"))
source(here("code/functions/niche_limits.R"))
source(here("code/functions/exposure.R"))
source(here("code/functions/exposure_times.R"))


# Module 1: Prepare climate data
# Read climate data
path <- "/Users/andreasschwarzmeyer/Library/CloudStorage/Dropbox/Projects/species-overshoot/raw_data/climate_data/CMIP6/Overshoot/ssp585/CNRM-ESM2-1"

climate_data <-  
  rast(here(path, "tas_Amon_CNRM-ESM2-1_ssp585_r1i1p1f2_gr_201501-210012.nc"))

# Transforming the climate before running the function.
climate_data <- climate_data - 273.15 # Kelvin to Celsius
climate_data <- rotate(climate_data) # Rotate raster so longitudes vary between -180 and 180

# The data come in a monthly format. We can transform it to yearly.
climate_data_yearly <- tapp(climate_data, 
                            index = "years", 
                            fun = "mean")


# So now we have two climate rasters: a monthly and a yearly.

# To transform the raster into a tibble, we need a raster template.
raster_template <- rast(extent = c(-180, 180, -90, 90), 
                        resolution = 1, 
                        crs = "EPSG:4326")

values(raster_template) <- 1:ncell(raster_template)



# Now we use the `raster_to_tibble` function to transform the raster data into a tibble.
# Monthly data
climate_tbl_monthly <- raster_to_tibble(climate.data = climate_data,
                                        raster.template = raster_template,
                                        temporal.resolution = "monthly", 
                                        year.min = 2015,
                                        year.max = 2100,
                                        format = "long",
                                        data.as.integer = TRUE)

# Yearly data
climate_tbl_yearly <- raster_to_tibble(climate.data = climate_data_yearly,
                                       raster.template = raster_template,
                                       temporal.resolution = "yearly", 
                                       year.min = 2015,
                                       year.max = 2100,
                                       format = "long",
                                       data.as.integer = TRUE)




array_data <- readRDS("/Users/andreasschwarzmeyer/Downloads/DCPP/OBS_tas_ERA5.rds")


array_tbl <- array_to_tibble(array.data = array_data, 
                             raster.template = raster_template,
                             latitude = "lat",
                             longitude = "lon",
                             data.as.integer = TRUE)

## 2. Module 2: Prepare species data

path <- "/Users/andreasschwarzmeyer/Library/CloudStorage/Dropbox/Projects/biodiversity-horizons/data-raw/"

species_shapefiles <- readRDS(here(path, "primates_shapefiles.rds"))




species_names <- unique(species_shapefiles$sci_name)

n_cores <- 11
species_data <- mclapply(species_names, 
                         function(.x) get_range_cells(species.names = .x, 
                                                      species.ranges = species_shapefiles, 
                                                      raster.template = raster_template), 
                         mc.cores = n_cores) 

names(species_data) <- species_names


## 3. Module 3: Estimate niche limits
 
# Number of standard deviations to define outliers
sd.threshold <- 3


climate_limits_monthly <- climate_limits(climate_tbl_monthly, 
                                         sd.threshold = sd.threshold, 
                                         temporal.resolution = "monthly",
                                         year.max = 2040) 


climate_limits_yearly <- climate_limits(climate_tbl_yearly, 
                                        sd.threshold = sd.threshold, 
                                        temporal.resolution = "yearly",
                                        year.max = 2040) 






niche_lim_monthly <- mclapply(species_data,
                              function(.x) niche_limits(.x,
                                                        climate.data = climate_limits_monthly,
                                                        percentiles = c(0, 1),
                                                        sd.threshold = sd.threshold),
                              mc.cores = n_cores)


niche_lim_yearly <- mclapply(species_data,
                             function(.x) niche_limits(.x,
                                                       climate.data = climate_limits_yearly,
                                                       percentiles = c(0.05, 0.95),
                                                       sd.threshold = sd.threshold),
                             mc.cores = n_cores)


# Converting the results from list to data frame. 


niche_lim_monthly <- niche_lim_monthly |>
  bind_rows(.id = "species") |>
  mutate(species = factor(species))

niche_lim_yearly <- niche_lim_yearly |>
  bind_rows(.id = "species") |>
  mutate(species = factor(species))

## 4. Module 4: Estimating exposure

exposure_monthly <- mclapply(species_names, 
                             function(.x) exposure(.x, 
                                                   species.data = species_data,
                                                   climate.data = climate_tbl_monthly,
                                                   niche.data = niche_lim_monthly,
                                                   month.specific = TRUE,
                                                   return.magnitude = TRUE),
                             mc.cores = 5) 


exposure_monthly <- exposure_monthly |> 
  bind_rows() |> 
  as_tibble() 




# Running exposure using yearly data, and selecting the wide format output. 

exposure_yearly <- mclapply(species_names, 
                            function(.x) exposure(.x, 
                                                  species.data = species_data,
                                                  climate.data = climate_tbl_yearly,
                                                  niche.data = niche_lim_yearly,
                                                  long.format = TRUE,
                                                  return.magnitude = FALSE),
                            mc.cores = 5) 



# In this example, the number of species isn't large. 
# But when running analyses with many species, the output size can be huge. 
# To deal with this, I've been splitting the species data into chunks and saving the results as `.parquet` files. Below is an example.

# First, set a directory to save the results:
  
mydir <- getwd()

# Run

species_names <- names(species_data)
chunk_size <- 100
species_chunks <- split(species_names, ceiling(seq_along(species_names) / chunk_size))


for(i in seq_along(species_chunks)){
  
  exposure_results <- mclapply(species_chunks[[i]], 
                               function(.x) exposure(.x, 
                                                     species.data = species_data,
                                                     climate.data = climate_tbl_monthly,
                                                     niche.data = niche_lim_monthly,
                                                     month.specific = TRUE,
                                                     return.magnitude = TRUE),
                               mc.cores = 7) 
  
  exposure_results <- exposure_results |> 
    bind_rows() |> 
    as_tibble() 
  
  
  file_name <- glue("{mydir}/chunk_{sprintf('%03d', i)}.parquet")
  write_parquet(exposure_results, file_name)
  
}



  
# Module 5: Calculating exposure times
  
# This step estimates where and when exposure occurs by computing, for every species in every grid cell, the year exposure starts and ends. This step currently only makes sense when using yearly climate data. For monthly climate data, we will calculate indicators (see below). Therefore, I'll use yearly data here.


result <- mclapply(X = exposure_yearly, 
                   FUN = function(.x){
                     apply(X = .x, 
                           MARGIN = 1, 
                           FUN = exposure_times, 
                           consecutive.elements = 5, 
                           first.year = 2015)
                   },
                   mc.cores = n_cores)

result <- result |> 
  bind_rows() |> 
  drop_na(exposure) 

