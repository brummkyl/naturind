#' Ecosystem Services
#'
#' @param ll A character vector containing file paths of raster files (e.g., spatial prioritization outputs) in '.tif' format.
#' @param puid A SpatRaster layer with planning units
#' @param all.pa A SpatRaster layer with all existing protected areas
#' @param st.pa A SpatRaster layer with all existing strict protected areas
#' @param n2k A SpatRaster layer with all existing Natura 2000 protected areas (or any other unique subset of PAs)
#' @param countries An sf vector layer containing individual countries
#' @param bioregion An sf vector layer containing individual biogeographic regions
#' @param pollination A SpatRaster layer quantifying pollination supply in units of m^2
#' @param capdem A SpatRaster layer quantifying species-based ES supply
#' @param project.name A character vector specifying the name of your project
#' @param version A character vector specifying the version (e.g., date) of your workflow
#'
#' @returns A formatted .csv file saved to the "performance_outputs/version" directory
#' @export
#'
#' @examples
#' \dontrun{
#' pi.es(ll = "C:/.../solution_file.tif", puid, all.pa, st.pa, n2k, countries,
#'        bioregion, pollination, capdem, project.name = "Example Project",
#'        version = "Version 1")
#'}

pi.es <- function(ll, puid, all.pa, st.pa, n2k, countries, bioregion, pollination = rs3, capdem = sp.reg.es, project.name, version){

  rs3 <- pollination
  sp.reg.es <- capdem

  # Transform to ensure alignment of CRS, and run assertions
  if(terra::crs(puid) != terra::crs(rs3)){
    rs3 <- rs3 |> terra::project(puid)
  }
  assertthat::assert_that(crs(rs3) == crs(puid))

  if(terra::crs(puid) != terra::crs(bioregion)){
    bioregion <- bioregion |> sf::st_transform(crs = sf::st_crs(puid))
  }
  assertthat::assert_that(sf::st_crs(bioregion) == sf::st_crs(puid))

  if(terra::crs(puid) != terra::crs(countries)){
    countries <- countries |> sf::st_transform(crs = sf::st_crs(puid))
  }
  assertthat::assert_that(sf::st_crs(countries) == sf::st_crs(puid))

  # Resample to align resolution and origin of ES and puid rasters
  if(any(terra::res(rs3) != terra::res(puid))) {
    rs3 <- terra::resample(rs3,puid,method="bilinear")
  }
  assertthat::assert_that(all.equal(res(rs3), res(puid)))

  puid <- puid %>% mask(countries)
  puid_area <- puid

  ee <- ll %>% rast() %>% extend(puid) %>% mask(countries)
  existing_pa_names <- c(basename(ll))

  ee <- c(ee, all.pa, st.pa)
  existing_pa_names <- c(existing_pa_names, "Existing all PA", "Existing strict PA")

  # Make an area grid - ha
  puid <- puid %>% mask(ee[[1]])
  puid_area <- puid %>% mask(ee[[1]])
  puid_area_ha <- puid_area * 100
  puid_area_m <- puid_area * 1e6

  # Check whether all layers already have the same NA footprint as ee[[1]]
  # This is used to adjust baseline existing PA calculations, for consistency with the provided planning solutions
  # One use case, for example, is to allow shifting between Europe and EU extents
  reference_na <- is.na(values(ee[[1]]))

  na_match <- all(vapply(
    seq_len(nlyr(ee)),
    function(i) {
      identical(reference_na, is.na(values(ee[[i]])))
    },
    logical(1)
  ))

  # Only mask the stack if at least one layer differs
  if (!na_match) {
    message("NA footprints differ between ee layers; masking stack to align with the first input layer...")
    ee <- ee %>% mask(ee[[1]])
  }

  #ee <- ee %>% mask(puid)

  message(paste0("Number of planning units in solution files = ", sum(!is.na(values(ee[[1]])))))
  message(paste0("Number of planning units in puid file = ", sum(!is.na(values(puid)))))

  #ee <- ee %>% mask(ee[[1]])

  rs3 <- rs3 %>% extend(puid)
  rs3 <- rs3 %>% crop(ee[[1]]) %>% mask(ee[[1]])

  ee_wrapped <- terra::wrap(ee)

  puid_area_m_wrapped <- terra::wrap(puid_area_m)
  rs3_wrapped <- terra::wrap(rs3)

  num_cores <- 1  # Adjust the number of cores as needed
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  clusterEvalQ(cl, {
    library(terra)
    library(exactextractr)
    library(sf)
    library(dplyr)
    library(tictoc)
    library(foreach)
  })

  #--- Export necessary variables to cluster workers ---
  clusterExport(cl, c("ee_wrapped", "bioregion", "countries",
                      "version", "puid_area_m_wrapped", "rs3_wrapped", "existing_pa_names",
                      "project.name"),
                envir = environment())

  all_results <- list()

  gc()

  tic()
  #rs.europe.results <- data.frame()
  rs.europe.results <- foreach(i = 1:nlyr(ee), .packages = c("terra", "exactextractr", "dplyr", "sf", "tictoc", "assertthat", "foreach")) %dopar% {
    #for(i in 1){
    # loop through the solution rasters
    ee <- terra::unwrap(ee_wrapped)
    pa.sol <- ee[[i]]
    file.name <- existing_pa_names[[i]]

    puid_area_m <- terra::unwrap(puid_area_m_wrapped)
    rs3 <- terra::unwrap(rs3_wrapped)

    rs.europe <- data.frame()
    rs.bioregion <- data.frame()
    rs.countries <- data.frame()

    # multiply the proportional solution by the planning unit area and ES layer
    pa.sol.rs3 <- pa.sol * puid_area_m * rs3
    pa.ov.rs3 <- puid_area_m * rs3

    # Overall sum
    rs3.sol <- global(pa.sol.rs3,"sum",na.rm=TRUE)[[1]]
    rs3.ov <- global(pa.ov.rs3,"sum",na.rm=TRUE)[[1]]

    # Overall normalization
    rs3.ov.prop <- (rs3.sol/rs3.ov)[[1]]

    o <- data.frame(project = paste0(project.name,"__",version),
                    scenario = file.name,
                    region = "Europe",
                    variable = "Pollination conserved by the solution",
                    variable_variant = "Total available in region",
                    unit = "%",
                    value = rs3.ov.prop * 100)

    rs.europe <- rbind(rs.europe, o)

    for (j in 1:nrow(bioregion)) {
      # loop through the bioregions
      region1 <- bioregion[j,] %>% sf::st_make_valid()

      # now we run zonal stats to quantify the selected and total available ES value in each region
      rs3.sol <- exactextractr::exact_extract(pa.sol.rs3, region1, fun="sum")
      rs3.ov <- exactextractr::exact_extract(pa.ov.rs3, region1, fun="sum")

      # and here we normalize by the total available ES within each region
      rs3.out.prop <- rs3.sol / rs3.ov

      o.bioregion <- data.frame(project = paste0(project.name,"__",version),
                                scenario = file.name,
                                region = paste0(as.character(bioregion[[j,"code"]])),
                                variable = "Pollination conserved by the solution",
                                variable_variant = "Total available in region",
                                unit = "%",
                                value = rs3.out.prop * 100)

      rs.bioregion <- rbind(rs.bioregion, o.bioregion)
    }

    gc()

    for (j in 1:nrow(countries)) {
      # loop through the bioregions
      region1 <- countries[j,] %>% sf::st_make_valid()

      # now we run zonal stats to quantify the selected and total available ES value in each region
      rs3.sol <- exactextractr::exact_extract(pa.sol.rs3, region1, fun="sum")
      rs3.ov <- exactextractr::exact_extract(pa.ov.rs3, region1, fun="sum")

      # and here we normalize by the total available ES within each region
      rs3.out.prop <- rs3.sol / rs3.ov

      o.countries <- data.frame(project = paste0(project.name,"__",version),
                                scenario = file.name,
                                region = paste0(as.character(countries[[j,"NAME_0"]])),
                                variable = "Pollination conserved by the solution",
                                variable_variant = "Total available in region",
                                unit = "%",
                                value = rs3.out.prop * 100)

      rs.countries <- rbind(rs.countries, o.countries)
    }

    gc()

    list(europe = rs.europe, bioregion = rs.bioregion, country = rs.countries)
  }
  all_results <- c(all_results, rs.europe.results)

  stopCluster(cl)

  rs.out.europe <- do.call(rbind, lapply(all_results, `[[`, "europe"))
  rs.out.bioregion <- do.call(rbind, lapply(all_results, `[[`, "bioregion"))
  rs.out.countries   <- do.call(rbind, lapply(all_results, `[[`, "country"))

  toc()

  sp.reg.es.aligned <- list()
  target_crs <- terra::crs(puid)  # Get the target CRS once

  for (k in 1:length(sp.reg.es)) {
    sp.es.lyr <- terra::rast(sp.reg.es[[k]])

    # Check and reproject if needed
    if(terra::crs(sp.es.lyr) != target_crs){
      sp.es.lyr <- sp.es.lyr |> terra::project(target_crs)
    }
    assertthat::assert_that(crs(sp.es.lyr) == target_crs)

    sp.es.lyr <- sp.es.lyr %>% extend(puid)
    sp.es.lyr <- sp.es.lyr %>% crop(ee[[1]]) %>% mask(ee[[1]])

    # Wrap for parallel processing
    sp.reg.es.aligned[[k]] <- terra::wrap(sp.es.lyr)
  }

  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  clusterEvalQ(cl, {
    library(terra)
    library(exactextractr)
    library(sf)
    library(dplyr)
    library(tictoc)
    library(foreach)
  })


  #--- Export necessary variables to cluster workers ---
  clusterExport(cl, c("ee_wrapped", "bioregion", "countries",
                      "version", "sp.reg.es", "sp.reg.es.aligned", "existing_pa_names", "project.name"),
                envir = environment())

  all_results2 <- list()

  gc()

  tic()
  sp.rs.europe.results <- foreach(i = 1:nlyr(ee), .packages = c("terra", "exactextractr", "dplyr", "sf", "tictoc", "assertthat", "foreach")) %dopar% {

    ee <- terra::unwrap(ee_wrapped)
    pa.sol <- ee[[i]]
    file.name <- existing_pa_names[[i]]

    sp.rs.europe <- data.frame()
    sp.rs.bioregion <- data.frame()
    sp.rs.countries <- data.frame()

    for (k in 1:length(sp.reg.es.aligned)){
      sp.es.lyr <- terra::unwrap(sp.reg.es.aligned[[k]])
      sp.es.lyr.sol <- pa.sol * sp.es.lyr

      tot.sol <- global(sp.es.lyr.sol,"sum",na.rm=TRUE)[[1]]
      tot.ov <- global(sp.es.lyr,"sum",na.rm=TRUE)[[1]]
      ov.prop <- (tot.sol/tot.ov)[[1]]

      o2 <- data.frame(project = paste0(project.name,"__",version),
                       scenario = file.name,
                       region = "Europe",
                       variable = paste0(tools::file_path_sans_ext(basename(sp.reg.es[[k]]))),
                       variable_variant = "Total available in region",
                       unit = "%",
                       value = ov.prop * 100)

      sp.rs.europe <- rbind(sp.rs.europe, o2)


      for (j in 1:nrow(bioregion)) {
        # loop through the bioregions
        region1 <- bioregion[j,] %>% sf::st_make_valid()

        # now we run zonal stats to quantify the selected and total available ES value in each region
        reg.sol <- exactextractr::exact_extract(sp.es.lyr.sol, region1, fun="sum")
        reg.ov <- exactextractr::exact_extract(sp.es.lyr, region1, fun="sum")

        reg.out.prop <- (reg.sol / reg.ov)[[1]]

        # and here we normalize by the total available ES within each region

        o.bioregion2 <- data.frame(project = paste0(project.name,"__",version),
                                   scenario = file.name,
                                   region = paste0(as.character(bioregion[[j,"code"]])),
                                   variable = paste0(tools::file_path_sans_ext(basename(sp.reg.es[[k]]))),
                                   variable_variant = "Total available in region",
                                   unit = "%",
                                   value = reg.out.prop * 100)

        sp.rs.bioregion <- rbind(sp.rs.bioregion, o.bioregion2)
      }

      gc()

      for (j in 1:nrow(countries)) {
        # loop through the bioregions
        region1 <- countries[j,] %>% sf::st_make_valid()

        # now we run zonal stats to quantify the selected and total available ES value in each region
        reg.sol <- exactextractr::exact_extract(sp.es.lyr.sol, region1, fun="sum")
        reg.ov <- exactextractr::exact_extract(sp.es.lyr, region1, fun="sum")

        reg.out.prop <- (reg.sol / reg.ov)[[1]]

        # and here we normalize by the total available ES within each region

        o.countries2 <- data.frame(project = paste0(project.name,"__",version),
                                   scenario = file.name,
                                   region = paste0(as.character(countries[[j,"NAME_0"]])),
                                   variable = paste0(tools::file_path_sans_ext(basename(sp.reg.es[[k]]))),
                                   variable_variant = "Total available in region",
                                   unit = "%",
                                   value = reg.out.prop * 100)

        sp.rs.countries <- rbind(sp.rs.countries, o.countries2)

      }
      gc()
    }
    list(europe = sp.rs.europe, bioregion = sp.rs.bioregion, country = sp.rs.countries)

  }
  all_results2 <- c(all_results2, sp.rs.europe.results)

  stopCluster(cl)

  sp.rs.out.europe <- do.call(rbind, lapply(all_results2, `[[`, "europe"))
  sp.rs.out.bioregion <- do.call(rbind, lapply(all_results2, `[[`, "bioregion"))
  sp.rs.out.countries   <- do.call(rbind, lapply(all_results2, `[[`, "country"))

  toc()

  # Bind results for europe, bioregions and countries, then export to output directory
  out.rs.final <- rbind(rs.out.europe, sp.rs.out.europe, rs.out.bioregion, sp.rs.out.bioregion, rs.out.countries, sp.rs.out.countries)

  outdir <- paste0("performance_outputs/", version)
  dir.create(outdir, showWarnings = FALSE)

  #if(consensus.sol == 1){
  #  ofname <- paste0(outdir.consensus, "/", project.name,"__",version, "__indicator_ecosystemservices.csv")
  #} else if (consensus.sol == 0){
    ofname <- paste0(outdir, "/", project.name,"__",version, "__indicator_ecosystemservices.csv")
  #}

  write.csv(out.rs.final, ofname, row.names = FALSE)

  message("Ecosystem service supply calculations done. Don't forget to create/update variable definitions")
}
