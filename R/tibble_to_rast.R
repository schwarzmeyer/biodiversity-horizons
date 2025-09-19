tibble_to_rast <- function(data, var, mask = FALSE, robin = TRUE){
  
  raster.template <- rast(extent = c(-180, 180, -90, 90), res = 1, crs = "EPSG:4326")
  values(raster.template) <- 1:ncell(raster.template)
  names(raster.template) <- "world_id"

  if(mask){
    
    world <- ne_countries(scale = "large", returnclass = "sf") |> 
      st_transform(crs = st_crs(raster.template))
    
    raster.template <- raster.template |> 
      terra::mask(world, touches = TRUE) 
    
  }
  
  ids <- raster.template |>
    as.data.frame(xy = T) |>
    as_tibble()

  r_plot <- data |>
    left_join(ids, by = "world_id") |>
    select(-world_id) |>
    select(x, y, !!sym(var))

  r <- rast(r_plot, type = "xyz", crs = "EPSG:4326")
  
  if(robin) r <- project(r, "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs")
  
  return(r)
  
}
