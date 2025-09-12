species_thresholds <- function(
    climate_data,          # 5D array: [member, lon, lat, year, month]
    species_data,          # 3D array: [lat, lon, species]
    percentiles = c(0.95, 0.99, 1.0), 
    n_cores = 1, 
    verbose = TRUE
) {
  
  # load required packages
  require(matrixStats)
  require(parallel)
  
  species_mask <- species_data[[1]]

  # input validation -------------------------------------
  
  # check array dimensions
  if (length(dim(climate_data)) != 5)  stop("climate_data must have dimensions [member, lon, lat, year, month]")
  if (length(dim(species_mask)) != 3)  stop("species_data must have dimensions [lat, lon, species]")
  
  # check dimension names
  if(!all(names(dim(climate_data)) %in% c("member","lon","lat","month","year"))) stop("Dimension names should be 'member', 'lon', 'lat', 'year', 'month'")
  if(!all(names(dim(species_mask)) %in% c("lon","lat","species"))) stop("Dimension names should be 'lon', 'lat', 'species'")
  
  # check spatial alignment
  species_mask <- aperm(species_mask, c(
    which(names(dim(species_mask)) == "lon"),
    which(names(dim(species_mask)) == "lat"),
    which(names(dim(species_mask)) == "species")))
  
  
  if (any(dim(species_mask)[1:2] != dim(climate_data)[2:3])) stop("Spatial dimensions mismatch between climate and species data")
  
  
  
  
  # compute climate percentiles --------------------------
  
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl))
  
  clusterExport(cl, "percentiles", envir = environment())
  
  if (verbose) message("1/4 Computing climate percentiles...")
  
  climate_percentiles <- parApply(climate_data, c(1, 2, 3, 5), function(x) {
      q <- quantile(x, percentiles, na.rm = TRUE, names = FALSE, type = 7)
      return(q)
    }, cl = cl) 
  
  
  
  climate_percentiles <- array(climate_percentiles, 
                  dim = list(
                  percentile = length(percentiles),
                  member = dim(climate_percentiles)[names(dim(climate_percentiles)) == "member"],
                  lon = dim(climate_percentiles)[names(dim(climate_percentiles)) == "lon"],
                  lat = dim(climate_percentiles)[names(dim(climate_percentiles)) == "lat"],
                  month = dim(climate_percentiles)[names(dim(climate_percentiles)) == "month"]))


  
  if (verbose) message("2/4 Reshaping data matrices...")
  
  percent <- which(names(dim(climate_percentiles)) == "percentile")
  member <- which(names(dim(climate_percentiles)) == "member")
  month <- which(names(dim(climate_percentiles)) == "month")
  lon <- which(names(dim(climate_percentiles)) == "lon")
  lat <- which(names(dim(climate_percentiles)) == "lat")
  
  climate_reshape <-  aperm(climate_percentiles, c(percent, member, month, lon, lat)) 
  
  dim(climate_reshape) <- c(prod(dim(climate_reshape)[1],
                                 dim(climate_reshape)[2],
                                 dim(climate_reshape)[3]),
                            prod(dim(climate_reshape)[4],
                                 dim(climate_reshape)[5]))
  
  species <- which(names(dim(species_mask)) == "species")
  lon_species <- which(names(dim(species_mask)) == "lon")
  lat_species <- which(names(dim(species_mask)) == "lat")
  
  dim(species_mask) <- c(prod(dim(species_mask)[lon_species],
                              dim(species_mask)[lat_species]),
                         dim(species_mask)[species])

    
  # export required data to workers
  clusterExport(cl, c("climate_reshape", "species_mask"), envir = environment())
  clusterEvalQ(cl, library(matrixStats))
  
  n_species <- dim(species_mask)[2]
    
  if (verbose) message("3/4 Processing species (", n_species, " total)...")
  
  # process species in parallel --
  result <- parLapply(cl, 1:n_species, function(s) {
    
    species_presence <- which(species_mask[, s] == 1)
    res <- rowMaxs(climate_reshape[, species_presence, drop = FALSE], na.rm = TRUE)  
      
    return(res)
    }) 
    
  if (verbose) message("4/4 Finalizing output...")
  
    result <- simplify2array(result)  
  
    
    result_array <- array(result, dim = list(
            percentile = dim(climate_percentiles)[percent],
            member = dim(climate_percentiles)[member],
            month = dim(climate_percentiles)[month],
            species = n_species
          ))
    
    output <- list(
      data = result_array,
      percentiles = percentiles,
      species = species_data[[4]]
    )
    
    return(output)
    
    
}




historical <- array(climate$hist[,,,1:30,],
                    dim = list(
                      member = 1,
                      lon = 360,
                      lat = 180,
                      year = 30,
                      month = 12
                    ))

fut <- array(climate$hist[,,,31:54,],
                    dim = list(
                      member = 1,
                      lon = 360,
                      lat = 180,
                      year = 30,
                      month = 12
                    ))


dim(historical)


result <- species_thresholds(
  climate_data = historical,  # [member, lon, lat, year, month]
  species_data = species_arrays,  # [lat, lon, species]
  percentiles = c(0.95, 1),
  n_cores = parallel::detectCores() - 1,
  verbose = TRUE
)

names(result)


dim(result)


future_climate <- fut
thermal_data <- result
species_mask <- species_arrays$data
selected_percentile <- 1

compute_exceedances <- function(
    future_climate,        # 5D array [member, lon, lat, year, month]
    thermal_data,     # 4D array [percentile, member, month, species]
    species_mask,          # 3D array [lat, lon, species] (0/1)
    selected_percentile = NULL,  # NULL = all percentiles
    n_cores = 1
) {
  
  thermal_tolerance <- thermal_data[[1]]
  
  
  # ------------------------------------------
  # Step 1: Input Validation
  # ------------------------------------------
  # Check array dimensions
  if (length(dim(future_climate)) != 5) stop("future_climate must be 5D: [member, lon, lat, year, month]")
  if (length(dim(thermal_tolerance)) != 4) stop("thermal_tolerance must be 4D: [percentile, member, month, species]")
  if (length(dim(species_mask)) != 3) stop("species_mask must be 3D: [lat, lon, species]")
  
  # check dimension names
  if(!all(names(dim(future_climate)) %in% c("member","lon","lat","month","year"))) stop("Dimension names should be 'member', 'lon', 'lat', 'year', 'month'")
  if(!all(names(dim(thermal_tolerance)) %in% c("percentile","member","month","species"))) stop("Dimension names should be 'percentile', 'member', 'month', 'species'")
  if(!all(names(dim(species_mask)) %in% c("lon","lat","species"))) stop("Dimension names should be 'lon', 'lat', 'species'")
  
  
  # Check spatial alignment

  # Get dimensions
  c(n_member, n_lon, n_lat, n_year, n_month) %<-% dim(future_climate)
  c(n_pctl, n_member_t, n_month_t, n_species) %<-% dim(thermal_tolerance)
  
  # Validate dimension consistency
  if(n_member != n_member_t || n_month != n_month_t) {
    stop(paste("Member/month dimensions mismatch between",
               "future_climate and thermal_tolerance"))
  }
  
  # --------------------------------------------------------------------------
  # Step 2: Percentile Selection & Tolerance Preparation
  # --------------------------------------------------------------------------
  # Handle percentile selection
  if(is.null(selected_percentile)) {
    tolerance_array <- thermal_tolerance
    pct_used <- 1:n_pctl
  } else {
    pct_used <- which(thermal_data$percentiles == selected_percentile)
    tolerance_array <- thermal_tolerance[pct_used,,,, drop = FALSE]
    n_pctl <- length(pct_used)
  }
  
  # --------------------------------------------------------------------------
  # Step 3: Data Reshaping for Efficient Computation
  # --------------------------------------------------------------------------
  # Reshape climate data to [member, month, lon*lat, year]
  climate_flat <- array(
    aperm(future_climate, c(1, 5, 2, 3, 4)),
    dim = c(n_member, n_month, n_lon*n_lat, n_year)
  )
  
  
  # Reshape species mask to [lon*lat, species]
  species_flat <- array(
    aperm(species_mask, c(2, 1, 3)),
    dim = c(n_lon*n_lat, n_species)
  )
  
  # Reshape tolerance array to [percentile, member, month, species]
  tolerance_flat <- array(
    tolerance_array,
    dim = c(n_pctl, n_member, n_month, n_species)
  )
  
  # --------------------------------------------------------------------------
  # Step 4: Parallel Exceedance Calculation
  # --------------------------------------------------------------------------
  library(parallel)
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl))
  
  clusterExport(cl, c("climate_flat", "species_flat", "tolerance_flat"))
  
  # exceedance_list <- parLapply(cl, 1:n_species, function(s) {
    exceedance_list <- lapply(1:n_species, function(s) {
      
    presence_mask <- species_flat[, s] == 1
    n_cells <- sum(presence_mask)
    if(n_cells == 0) return(NULL)
    
    # Initialize result array [pctl, member, month, cell, year]
    res <- array(0, dim = c(n_pctl, n_member, n_month, n_cells, n_year))
    
    for(p in 1:n_pctl) {
      for(m in 1:n_month) {
        thresholds <- tolerance_flat[p, , m, s]  # [member]
        
        # Get climate data [member, cell, year]
        climate_slice <- climate_flat[, m, presence_mask, ]
        
        # Compare and store cell-level results
        res[p, , m, , ] <- climate_slice > thresholds
      }
    }
    res
  })
    
    # --------------------------------------------------------------------------
  
    # Metadata structure
    result <- list(
      dims = list(
        percentile = pct_names,
        member = member_names,
        lat = lat_coords,
        lon = lon_coords,
        month = month_names,
        species = species_names,
        year = year_sequence
      ),
      # Store presence cells and their exceedance data
      data = lapply(1:n_species, function(s) {
        cells <- which(species_mask[,,s] == 1, arr.ind = TRUE)
        if (nrow(cells) == 0) return(NULL)
        
        list(
          cells = cells,  # Matrix of [lat, lon] indices
          # 5D array: [pctl, member, month, cell_idx, year]
          values = exceedance_list[[s]]
        )
      })
    )
    
    
    
    convert_to_7d <- function(exceedance_list, species_mask, pct_names) {
      # Get dimensions from inputs
      n_pctl <- dim(exceedance_list[[1]])[1]
      n_member <- dim(exceedance_list[[1]])[2]
      n_month <- dim(exceedance_list[[1]])[3]
      n_year <- dim(exceedance_list[[1]])[5]
      n_species <- dim(species_mask)[3]
      n_lat <- dim(species_mask)[1]
      n_lon <- dim(species_mask)[2]
      
      # Initialize 7D array with NA
      result_array <- array(NA, 
                              dim = list(
                              percentile = n_pctl,
                              member = n_member,
                              lat = n_lat,
                              lon = n_lon,
                              month = n_month,
                              species = n_species,
                              year = n_year
                            ))
      
      # Get presence cells for each species
      presence_cells <- lapply(1:n_species, function(s) {
        which(species_mask[,,s] == 1, arr.ind = TRUE)
      })
      
      # Populate the 7D array
      for (s in 1:n_species) {
        if (is.null(exceedance_list[[s]])) next
        
        cells <- presence_cells[[s]]
        if (nrow(cells) == 0) next
        
        # For each presence cell in this species' range
        for (cell_idx in 1:nrow(cells)) {
          lat <- cells[cell_idx, 1]
          lon <- cells[cell_idx, 2]
          
          # Extract values for this cell across all dimensions
          # exceedance_list[[s]] dimensions: [pctl, member, month, cell, year]
          result_array[,,lat,lon,,s,] <- exceedance_list[[s]][,,,cell_idx,]
        }
      }
      
      return(result_array)
    }
    
  # --------------------------------------------------------------------------
  # Step 5: Result Consolidation
  # --------------------------------------------------------------------------
  # Convert list to array [percentile, member, month, year, species]
  final_result <-     list(
    # Metadata
    dims = list(
      percentile = dimnames(thermal_tolerance)[[1]][pct_idx],
      member = dimnames(future_climate)[[1]],
      month = month.abb[1:n_month],
      year = dimnames(future_climate)[[4]],
      species = dimnames(species_mask)[[3]]
    ),
    
    # Cell coordinates for each species
    presence_cells = lapply(1:n_species, function(s) {
      which(species_mask[,,s] == 1, arr.ind = TRUE)  # [lat, lon] indices
    }),
    
    # Exceedance data (list of 5D arrays)
    data = exceedance_list
  )
  return(final_array)
}