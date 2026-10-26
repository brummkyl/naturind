#' Mean/Median Patch Size
#'
#' @param ll A character vector containing file paths of raster files (e.g., spatial prioritization outputs) in '.tif' format.
#' @param puid A SpatRaster layer with planning units
#' @param countries An sf vector layer containing individual countries
#' @param bioregion An sf vector layer containing individual biogeographic regions
#' @param threshold A double numeric, ranging between 0 and 1, defining the threshold for binarizing a proportional protected area raster
#' @param project.name A character vector specifying the name of your project
#' @param version A character vector specifying the version (e.g., date) of your workflow
#'
#' @returns A formatted .csv file saved to the "performance_outputs/version" directory
#' @export
#'
#' @examples
#' \dontrun{
#' pi.mean.patch.size(ll = "C:/.../solution_file.tif", puid, countries, bioregion,
#'                     threshold = 0, project.name = "Example Project",
#'                     version = "Version 1")
#'}

pi.mean.patch.size <- function(ll, puid, countries, bioregion, threshold = 0, project.name, version){

  prop_binary_thresh <- threshold

  # Transform to ensure alignment of CRS
  # if(terra::crs(puid) != terra::crs(all.pa)){
  #   all.pa <- all.pa |> terra::project(puid)
  # }
  # assertthat::assert_that(crs(all.pa) == crs(puid))
  #
  # if(terra::crs(puid) != terra::crs(st.pa)){
  #   st.pa <- st.pa |> terra::project(puid)
  # }
  # assertthat::assert_that(crs(st.pa) == crs(puid))
  #
  # if(terra::crs(puid) != terra::crs(n2k)){
  #   n2k <- n2k |> terra::project(puid)
  # }
  # assertthat::assert_that(crs(n2k) == crs(puid))

  puid <- puid %>% mask(countries)
  puid_area <- puid

  ee <- ll %>% rast() %>% extend(puid) %>% mask(countries)
  existing_pa_names <- c(basename(ll))

  # # stack layers
  # ee <- c(ee, all.pa, st.pa, n2k)
  # existing_pa_names <- c(existing_pa_names, "Existing all PA", "Existing strict PA", "Existing N2K")

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
  #ee <- rast(ll) %>% mask(countries)
  #ee <- ee %>% mask(puid)

  co.dissolve <- sf::st_union(countries)

  message(paste0("Number of planning units in solution files = ", sum(!is.na(values(ee[[1]])))))
  message(paste0("Number of planning units in puid file = ", sum(!is.na(values(puid)))))

  # calculations

  size.europe.existing <- data.frame()
  size.bioregion.existing <- data.frame()
  size.countries.existing <- data.frame()
  pb <- progress::progress_bar$new(total = nlyr(ee))

  tic()
  for (i in 1:nlyr(ee)) {
    pb$tick()
    # loop through the solution rasters
    pa.sol_prop <- ee[[i]]
    file.name <- existing_pa_names[[i]]

    # Use the threshold only to decide which cells belong to a patch
    patch_mask <- terra::ifel(pa.sol_prop > prop_binary_thresh, 1, NA)

    pa.patches <- terra::patches(
      patch_mask,
      directions = 8,
      zeroAsNA = TRUE,
      allowGaps = TRUE
    )

    # Weight each cell's area by its proportional coverage
    #cell_area_km2 <- terra::cellSize(pa.sol_prop, unit = "km")
    area_km2 <- pa.sol_prop #* cell_area_km2
    names(area_km2) <- "area_km2"

    # Sum proportional area within each patch
    patch.areas <- terra::zonal(
      area_km2,
      pa.patches,
      fun = "sum",
      na.rm = TRUE
    )

    median.patch.size <- median(patch.areas$area_km2, na.rm = TRUE)
    mean.patch.size <- mean(patch.areas$area_km2, na.rm = TRUE)

    #
    # pa.sol <- ifel(pa.sol>prop_binary_thresh, 1, NA)
    # pa.patches <- patches(pa.sol, directions = 8, zeroAsNA = FALSE, allowGaps = TRUE)
    #
    # patch.areas <- expanse(pa.patches, unit = "km", byValue = TRUE)
    # median.patch.size <- median(patch.areas$area, na.rm=T)
    # mean.patch.size <- mean(patch.areas$area, na.rm=T)

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

      bi.pa.sol_prop <- terra::mask(pa.sol_prop, region1)
      bi.patch_mask <- terra::ifel(
        bi.pa.sol_prop > prop_binary_thresh,
        1,
        NA
      )

      bi.patches <- terra::patches(
        bi.patch_mask,
        directions = 8,
        zeroAsNA = TRUE,
        allowGaps = TRUE
      )

      bi.area <- bi.pa.sol_prop
      names(bi.area) <- "area_km2"

      patch.areas <- terra::zonal(
        bi.area,
        bi.patches,
        fun = "sum",
        na.rm = TRUE
      )

      if (nrow(patch.areas) > 0) {
        median.patch.size <- median(patch.areas$area_km2, na.rm = TRUE)
        mean.patch.size <- mean(patch.areas$area_km2, na.rm = TRUE)
      } else {
        median.patch.size <- NA_real_
        mean.patch.size <- NA_real_
      }

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
      gc()
    }

    for (j in 1:nrow(countries)) {
      # loop through the countries
      # TODO: See if we need to be clipping and masking region 1 to puid_area
      region1 <- countries[j,] %>% sf::st_make_valid()

      co.pa.sol_prop <- terra::mask(pa.sol_prop, region1)

      co.patch_mask <- terra::ifel(
        co.pa.sol_prop > prop_binary_thresh,
        1,
        NA
      )

      co.patches <- terra::patches(
        co.patch_mask,
        directions = 8,
        zeroAsNA = TRUE,
        allowGaps = TRUE
      )

      co.area <- co.pa.sol_prop
      names(co.area) <- "area_km2"

      patch.areas <- terra::zonal(
        co.area,
        co.patches,
        fun = "sum",
        na.rm = TRUE
      )

      if (nrow(patch.areas) > 0) {
        median.patch.size <- median(patch.areas$area_km2, na.rm = TRUE)
        mean.patch.size <- mean(patch.areas$area_km2, na.rm = TRUE)
      } else {
        median.patch.size <- NA_real_
        mean.patch.size <- NA_real_
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
      gc()
    }
    gc()
  }
  toc()

  # combine calculations for continental, bioregion, and country extents
  out.size.final <- rbind(size.europe.existing, size.bioregion.existing, size.countries.existing)

  outdir <- paste0("performance_outputs/", version)
  dir.create(outdir, showWarnings = FALSE)

  #if(consensus.sol == 1){
  #  ofname <- paste0(outdir.consensus, "/", project.name,"__",version, "__indicator_relativepatchsize.csv")
  #} else if (consensus.sol == 0){
    ofname <- paste0(outdir, "/", project.name,"__",version, "__indicator_relativepatchsize.csv")
  #}

  write.csv(out.size.final, ofname, row.names = FALSE)

  message("Patch size calculations done. Don't forget to create/update variable definitions")
}
