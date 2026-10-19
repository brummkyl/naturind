#' Conserved Climate Resilient Area
#'
#' @param ll A character vector containing file paths of raster files (e.g., spatial prioritization outputs) in '.tif' format.
#' @param puid A SpatRaster layer with planning units
#' @param all.pa A SpatRaster layer with all existing protected areas
#' @param st.pa A SpatRaster layer with all existing strict protected areas
#' @param n2k A SpatRaster layer with all existing Natura 2000 protected areas (or any other unique subset of PAs)
#' @param countries An sf vector layer containing individual countries
#' @param bioregion An sf vector layer containing individual biogeographic regions
#' @param clim.res A continuous SpatRaster layer with larger values representing areas that are more climate resilient
#' @param project.name A character vector specifying the name of your project
#' @param version A character vector specifying the version (e.g., date) of your workflow
#'
#' @returns A formatted .csv file saved to the "performance_outputs/version" directory
#' @export
#'
#' @examples
#' pi.climate(ll = "C:/.../solution_file.tif", puid, all.pa, st.pa, n2k,
#'             countries, bioregion, clim.res, project.name = "Example Project",
#'             version = "Version 1")

pi.climate <- function(ll, puid, all.pa, st.pa, n2k, countries, bioregion, clim.res, project.name, version){

  ee <- ll %>% rast() %>% extend(puid)
  existing_pa_names <- c(basename(ll))

  # Transform to ensure alignment of CRS

  if(terra::crs(puid_area) != terra::crs(clim.res)){
    clim.res <- clim.res |> terra::project(puid_area)
  }
  assertthat::assert_that(crs(clim.res) == crs(puid_area))

  if(terra::crs(puid) != terra::crs(all.pa)){
    all.pa <- all.pa |> terra::project(puid)
  }
  assertthat::assert_that(crs(all.pa) == crs(puid))

  if(terra::crs(puid) != terra::crs(st.pa)){
    st.pa <- st.pa |> terra::project(puid)
  }
  assertthat::assert_that(crs(st.pa) == crs(puid))

  if(terra::crs(puid) != terra::crs(n2k)){
    n2k <- n2k |> terra::project(puid)
  }
  assertthat::assert_that(crs(n2k) == crs(puid))


  # stack layers
  ee <- c(ee, all.pa, st.pa, n2k)
  existing_pa_names <- c(existing_pa_names, "Existing all PA", "Existing strict PA", "Existing N2K")

  # Make an area grid - ha
  puid <- puid %>% mask(ee[[1]])
  puid_area <- puid %>% mask(ee[[1]])

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

  #ee <- ee %>% mask(ee[[1]])

  clim.res <- clim.res %>% extend(puid)
  clim.res <- clim.res %>% crop(ee[[1]]) %>% mask(ee[[1]])

  # Resample to align resolution and origin of climate resilience and puid_area rasters
  # TODO: Add a check so that we only resample the climate raster if resolutions differ
  clim.res <- terra::resample(clim.res,puid,method="bilinear")
  assertthat::assert_that(all.equal(res(clim.res), res(puid)), msg = "Rasters differ in resolution")
  clim.res <- clim.res %>% terra::crop(puid) %>% terra::mask(puid)

  # extract any given quantile from the climate resilience rankmap, and then convert this to a binary format
  quant <- global(clim.res, \(i) quantile(i, climate.quantile, na.rm=T))
  clim.res.binary <- ifel(clim.res <= quant[[1]], NA, 1)
  #plot(clim.res.binary, col = "darkgreen", main = paste0("Top ", 100-(climate.quantile*100), "% climate resilience"), ext = ext(clim.res.binary)*1.1)


  #######################################################################################

  # EXISTING PROTECTED AREAS

  #######################################################################################

  # bioregions and countries
  climate.europe.existing <- data.frame()
  climate.bioregion.existing <- data.frame()
  climate.countries.existing <- data.frame()
  pb <- progress::progress_bar$new(total = nlyr(ee))

  tic()
  for (i in 1:nlyr(ee)) {
    pb$tick()
    # loop through the solution rasters
    pa.sol <- ee[[i]]
    file.name <- existing_pa_names[[i]]

    # calculate the selected area
    pa.sol.area <- pa.sol #* puid_area
    # calculate the climate resilient area
    climarea <- puid_area * clim.res.binary
    # calculate the selected climate resilient area
    pa.sol.clim <- pa.sol * clim.res.binary #*  puid_area

    # Overall sum
    ov.eur.area <- global(puid_area,"sum",na.rm=TRUE)[[1]]
    ov.sol.area <- global(pa.sol.area,"sum",na.rm=TRUE)[[1]]
    ov.clim.area <- global(climarea,"sum",na.rm=TRUE)[[1]]
    ov.clim.out <- global(pa.sol.clim,"sum",na.rm=TRUE)[[1]]
    # normalization
    # normalize by the total area of each region
    clim.out.norm <- (ov.clim.out / as.numeric(ov.eur.area))
    # normalize by the selected area within each region
    clim.out.norm2 <- (ov.clim.out / as.numeric(ov.sol.area))
    # normalize by the climate resilient area within each region
    clim.out.norm3 <- (ov.clim.out / as.numeric(ov.clim.area))

    o <- data.frame(project = paste0(project.name,"__",version),
                    scenario = file.name,
                    region = "Europe",
                    variable = "Climate resilient area selected in the solution",
                    variable_variant = "Total area of region",
                    unit = "%",
                    value = clim.out.norm)
    o2 <- o %>% mutate(
      variable_variant = "Selected area in region",
      unit = "%",
      value = clim.out.norm2)
    o3 <- o %>% mutate(
      variable_variant = "Climate resilient area in region",
      unit = "%",
      value = clim.out.norm3)

    climate.europe.existing <- rbind(climate.europe.existing, o, o2, o3)

    for (j in 1:nrow(bioregion)) {
      # loop through the bioregions
      region1 <- bioregion[j,] %>% sf::st_make_valid()

      # denominators for normalization
      # calculate the total area of each region, in km2
      region1.area <- exactextractr::exact_extract(puid_area, region1, fun="sum")
      # calculate the selected area within each region
      area.sel <- exactextractr::exact_extract(pa.sol.area, region1, fun="sum")
      # calculate the climate resilient area within each region
      climarea.sel <- exactextractr::exact_extract(climarea, region1, fun="sum")
      # numerator for normalization
      # calculate the selected climate resilient area within each region
      clim.out <- exactextractr::exact_extract(pa.sol.clim, region1, fun="sum")
      # normalization
      # normalize by the total area of each region
      clim.out.norm <- (clim.out / as.numeric(region1.area))
      # normalize by the selected area within each region
      clim.out.norm2 <- (clim.out / as.numeric(area.sel))
      # normalize by the climate resilient area within each region
      clim.out.norm3 <- (clim.out / as.numeric(climarea.sel))

      o <- data.frame(project = paste0(project.name,"__",version),
                      scenario = file.name,
                      region = paste0(as.character(bioregion[[j,"code"]])),
                      variable = "Climate resilient area selected in the solution",
                      variable_variant = "Total area of region",
                      unit = "%",
                      value = clim.out.norm)
      o2 <- o %>% mutate(
        variable_variant = "Selected area in region",
        unit = "%",
        value = clim.out.norm2)
      o3 <- o %>% mutate(
        variable_variant = "Climate resilient area in region",
        unit = "%",
        value = clim.out.norm3)

      climate.bioregion.existing <- rbind(climate.bioregion.existing, o, o2, o3)
    }

    for (j in 1:nrow(countries)) {
      # loop through the countries
      region1 <- countries[j,] %>% sf::st_make_valid()

      # denominators for normalization
      # calculate the total area of each region, in km2
      region1.area <- exactextractr::exact_extract(puid_area, region1, fun="sum")
      # calculate the selected area within each region
      area.sel <- exactextractr::exact_extract(pa.sol.area, region1, fun="sum")
      # calculate the climate resilient area within each region
      climarea.sel <- exactextractr::exact_extract(climarea, region1, fun="sum")
      # numerator for normalization
      # calculate the selected climate resilient area within each region
      clim.out <- exactextractr::exact_extract(pa.sol.clim, region1, fun="sum")
      # normalization
      # normalize by the total area of each region
      clim.out.norm <- (clim.out / as.numeric(region1.area))
      # normalize by the selected area within each region
      clim.out.norm2 <- (clim.out / as.numeric(area.sel))
      # normalize by the climate resilient area within each region
      clim.out.norm3 <- (clim.out / as.numeric(climarea.sel))

      o <- data.frame(project = paste0(project.name,"__",version),
                      scenario = file.name,
                      region = paste0(as.character(countries[[j,"NAME_0"]])),
                      variable = "Climate resilient area selected in the solution",
                      variable_variant = "Total area of region",
                      unit = "%",
                      value = clim.out.norm)
      o2 <- o %>% mutate(
        variable_variant = "Selected area in region",
        unit = "%",
        value = clim.out.norm2)
      o3 <- o %>% mutate(
        variable_variant = "Climate resilient area in region",
        unit = "%",
        value = clim.out.norm3)

      climate.countries.existing <- rbind(climate.countries.existing, o, o2, o3)
    }
  }
  toc()

  # Bind results for bioregions and countries, then export to output directory
  out.climate.final <- rbind(climate.europe.existing, climate.bioregion.existing, climate.countries.existing)
  # Convert from proportions to percentages
  out.climate.final$value <- out.climate.final$value * 100

  if(consensus.sol == 1){
    ofname <- paste0(outdir.consensus, "/", project.name,"__",version, "__indicator_climateresilience.csv")
  } else if (consensus.sol == 0){
    ofname <- paste0(outdir, "/", project.name,"__",version, "__indicator_climateresilience.csv")
  }

  write.csv(out.climate.final, ofname, row.names = FALSE)

  message("Climate resilience calculation done. Don't forget to create/update variable definitions")

}
