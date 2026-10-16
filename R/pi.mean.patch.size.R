#' Mean/Median Patch Size
#'
#' @param ll A character vector containing file paths of raster files (e.g., spatial prioritization outputs) in '.tif' format.
#' @param puid A SpatRaster layer with planning units
#' @param all.pa A SpatRaster layer with all existing protected areas
#' @param st.pa A SpatRaster layer with all existing strict protected areas
#' @param n2k A SpatRaster layer with all existing Natura 2000 protected areas (or any other unique subset of PAs)
#' @param countries An sf vector layer containing individual countries
#' @param bioregion An sf vector layer containing individual biogeographic regions
#' @param threshold Majority rule threshold for binarizing a proportional protected area raster
#' @param project.name A character vector specifying the name of your project
#' @param version A character vector specifying the version (e.g., date) of your workflow
#'
#' @returns A formatted .csv file saved to the "performance_outputs/version" directory
#' @export
#'
#' @examples
#' pi.mean.patch.size(ll = "C:/.../solution_file.tif", puid, all.pa, st.pa, n2k, countries, bioregion,
#'                     threshold = 0.5, project.name = "Example Project",
#'                     version = "Version 1")

pi.mean.patch.size <- function(ll, puid, all.pa, st.pa, n2k, countries, bioregion, threshold = prop_binary_thresh, project.name, version){

  prop_binary_thresh <- threshold

  ee <- ll %>% rast() %>% extend(puid)
  existing_pa_names <- c(basename(ll))

  # Transform to ensure alignment of CRS
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

  # Mask puid layer to the solution layer to ensure spatial alignment
  puid <- puid %>% mask(ee[[1]])

  # Check whether all layers already have the same NA footprint as ee[[1]] - i.e., the first solution layer
  # This is used to adjust baseline existing PA calculations, for consistency with the provided planning solutions
  # One use case, for example, is to allow dynamic shifting between Europe and EU extents
  reference_na <- is.na(values(ee[[1]]))

  na_match <- all(vapply(
    seq_len(nlyr(ee)),
    function(i) {identical(reference_na, is.na(values(ee[[i]])))},
    logical(1)
  ))

  # Only mask the stack if at least one layer differs
  if (!na_match) {
    message("NA footprints differ between ee layers; masking stack to align with the first input layer...")
    ee <- ee %>% mask(ee[[1]])
  }
  #ee <- list(all.pa, st.pa, n2k)
  ee <- rast(ll) %>% mask(countries)

  # calculations

  size.europe.existing <- data.frame()
  size.bioregion.existing <- data.frame()
  size.countries.existing <- data.frame()
  pb <- progress::progress_bar$new(total = nlyr(ee))

  tic()
  for (i in 1:nlyr(ee)) {
    pb$tick()
    # loop through the solution rasters
    pa.sol <- ee[[i]]
    file.name <- existing_pa_names[[i]]

    pa.sol <- ifel(pa.sol>prop_binary_thresh, 1, NA)
    pa.patches <- patches(pa.sol, directions = 8, zeroAsNA = FALSE, allowGaps = TRUE)

    patch.areas <- expanse(pa.patches, unit = "km", byValue = TRUE)
    median.patch.size <- median(patch.areas$area, na.rm=T)
    mean.patch.size <- mean(patch.areas$area, na.rm=T)

    o <- data.frame(project = paste0(project.name,"__",version),
                    scenario = file.name,
                    region = "Europe",
                    variable = "Median patch size",
                    variable_variant = "Median patch size in region",
                    unit = "km2",
                    value = median.patch.size)

    o2 <- o %>% mutate(variable = "Mean patch size",
                       variable_variant = "Mean patch size in region",
                       value = mean.patch.size)

    size.europe.existing <- rbind(size.europe.existing, o, o2)

    for (j in 1:nrow(bioregion)) {
      # loop through the bioregions
      # TODO: See if we need to be clipping and masking region 1 to puid_area
      region1 <- bioregion[j,] %>% sf::st_make_valid()

      bi.pa.sol <- pa.sol %>% mask(region1)

      pa.patches <- patches(bi.pa.sol, directions = 8, zeroAsNA = FALSE, allowGaps = FALSE)

      patch.areas <- expanse(pa.patches, unit = "km", byValue = TRUE)
      median.patch.size <- median(patch.areas$area, na.rm=T)
      mean.patch.size <- mean(patch.areas$area, na.rm=T)

      o <- data.frame(project = paste0(project.name,"__",version),
                      scenario = file.name,
                      region = paste0(as.character(bioregion[[j,"code"]])),
                      variable = "Median patch size",
                      variable_variant = "Median patch size in region",
                      unit = "km2",
                      value = median.patch.size)

      o2 <- o %>% mutate(variable = "Mean patch size",
                         variable_variant = "Mean patch size in region",
                         value = mean.patch.size)

      size.bioregion.existing <- rbind(size.bioregion.existing, o, o2)
    }

    for (j in 1:nrow(countries)) {
      # loop through the countries
      # TODO: See if we need to be clipping and masking region 1 to puid_area
      region1 <- countries[j,] %>% sf::st_make_valid()

      co.pa.sol <- pa.sol %>% mask(region1)

      if (!is.na(global(co.pa.sol,"sum",na.rm=T))){

        pa.patches <- patches(co.pa.sol, directions = 8, zeroAsNA = FALSE, allowGaps = FALSE)

        patch.areas <- expanse(pa.patches, unit = "km", byValue = TRUE)
        median.patch.size <- median(patch.areas$area, na.rm=T)
        mean.patch.size <- mean(patch.areas$area, na.rm=T)

      } else if (is.na(global(co.pa.sol,"sum",na.rm=T))) {
        median.patch.size = NA
        mean.patch.size = NA
      }

      o <- data.frame(project = paste0(project.name,"__",version),
                      scenario = file.name,
                      region = paste0(as.character(countries[[j,"NAME_0"]])),
                      variable = "Median patch size",
                      variable_variant = "Median patch size in region",
                      unit = "km2",
                      value = median.patch.size)

      o2 <- o %>% mutate(variable = "Mean patch size",
                         variable_variant = "Mean patch size in region",
                         value = mean.patch.size)

      size.countries.existing <- rbind(size.countries.existing, o, o2)
    }
  }
  toc()

  # without bioregions
  out.size.final <- rbind(size.europe.existing, size.countries.existing)

  if(consensus.sol == 1){
    ofname <- paste0(outdir.consensus, "/", project.name,"__",version, "__indicator_relativepatchsize.csv")
  } else if (consensus.sol == 0){
    ofname <- paste0(outdir, "/", project.name,"__",version, "__indicator_relativepatchsize.csv")
  }

  write.csv(out.size.final, ofname, row.names = FALSE)

  message("Patch size calculations done. Don't forget to create/update variable definitions")
}
