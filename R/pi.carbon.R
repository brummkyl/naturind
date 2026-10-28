#' Conserved Carbon Stock and Sequestration Potential
#'
#' @param ll A character vector containing file paths of raster files (e.g., spatial prioritization outputs) in '.tif' format.
#' @param puid A SpatRaster layer, reflecting the amount of land area per planning unit
#' @param countries An sf vector layer containing individual countries, with names specified in a column named "code"
#' @param bioregion An sf vector layer containing individual biogeographic regions, with names specified in a column named "code"
#' @param carbon.stock A SpatRaster layer with carbon stock data in units of tC/ha
#' @param carbon.sequestration A SpatRaster layer with carbon sequestration potential data in units of tC/ha/year
#' @param project.name A character vector specifying the name of your project
#' @param version A character vector specifying the version (e.g., date) of your workflow
#'
#' @returns A formatted .csv file saved to the "performance_outputs/version" directory
#' @export
#'
#' @examples
#' \dontrun{
#' pi.carbon(
#'   ll = "C:/.../solution_file.tif", puid,
#'   countries, bioregion, carbon.stock = carbon,
#'   carbon.sequestration = carbon.seq,
#'   project.name = "Example Project", version = "Version 1"
#' )
#' }
pi.carbon <- function(ll, puid, countries, bioregion, carbon.stock = carbon, carbon.sequestration = carbon.seq, project.name, version) {
  carbon <- carbon.stock
  carbon.seq <- carbon.sequestration

  # Transform to ensure alignment of CRS, and run assertions
  if (terra::crs(puid) != terra::crs(carbon)) {
    carbon <- carbon |> terra::project(puid)
  }
  assertthat::assert_that(crs(carbon) == crs(puid))

  if (terra::crs(puid) != terra::crs(carbon.seq)) {
    carbon.seq <- carbon.seq |> terra::project(puid)
  }
  assertthat::assert_that(crs(carbon.seq) == crs(puid))

  puid <- puid %>% mask(countries)
  puid_area <- puid

  ee <- ll %>%
    rast() %>%
    extend(puid) %>%
    mask(countries)
  existing_pa_names <- c(basename(ll))

  # Make an area grid - ha
  puid <- puid %>% mask(ee[[1]])
  puid_area <- puid %>% mask(ee[[1]])
  puid_area_ha <- puid_area * 100

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

  # ee <- ee %>% mask(puid)

  carbon <- carbon %>% mask(ee[[1]])
  carbon.seq <- carbon.seq %>% mask(ee[[1]])

  # Resample to align resolution and origin of carbon and puid_area_ha rasters
  # carbon storage
  if (any(terra::res(carbon) != terra::res(puid_area_ha))) {
    carbon <- terra::resample(carbon, puid_area_ha, method = "bilinear")
  }
  assertthat::assert_that(all.equal(res(carbon), res(puid_area_ha)))
  # carbon sequestration
  if (any(terra::res(carbon.seq) != terra::res(puid_area_ha))) {
    carbon.seq <- terra::resample(carbon.seq, puid_area_ha, method = "bilinear")
  }
  assertthat::assert_that(all.equal(res(carbon.seq), res(puid_area_ha)))

  # for this indicator, we only consider areas in which carbon sequestration is positive
  carbon.seq <- ifel(carbon.seq < 0, 0, carbon.seq)

  # union countries to represent the continental extent
  co.dissolve <- sf::st_union(countries)

  message(paste0("Number of planning units in solution files = ", sum(!is.na(values(ee[[1]])))))
  message(paste0("Number of planning units in puid file = ", sum(!is.na(values(puid)))))

  # output containers
  carbon.europe.existing <- data.frame()
  carbon.bioregion.existing <- data.frame()
  carbon.countries.existing <- data.frame()
  pb <- progress::progress_bar$new(total = nlyr(ee))

  tic()
  # loop through the solution rasters
  for (i in 1:nlyr(ee)) {
    pb$tick()
    pa.sol <- ee[[i]]
    file.name <- existing_pa_names[[i]]

    # multiply the proportional solution by the planning unit area by the carbon stock raster
    # output quantifies tonnes of carbon per cell conserved in the solution
    pa.sol.area <- pa.sol * 100 #* puid_area_ha
    pa.sol.carbon <- pa.sol * 100 * carbon
    pa.sol.carbon.seq <- pa.sol * 100 * carbon.seq
    ov.carbon <- puid_area_ha * carbon
    ov.carbon.seq <- puid_area_ha * carbon.seq

    # Overall sum
    ov.area <- exactextractr::exact_extract(pa.sol.area, co.dissolve, fun = "sum")
    ov.carbon.tot <- exactextractr::exact_extract(ov.carbon, co.dissolve, fun = "sum")
    ov.carbon.seq.tot <- exactextractr::exact_extract(ov.carbon.seq, co.dissolve, fun = "sum")
    carbon.ov <- exactextractr::exact_extract(pa.sol.carbon, co.dissolve, fun = "sum")
    carbon.seq.ov <- exactextractr::exact_extract(pa.sol.carbon.seq, co.dissolve, fun = "sum")

    # Overall normalization
    carbon.ov.norm <- (carbon.ov / ov.area)[[1]]
    carbon.seq.ov.norm <- (carbon.seq.ov / ov.area)[[1]]
    carbon.ov.prop <- (carbon.ov / ov.carbon.tot)[[1]]
    carbon.seq.ov.prop <- (carbon.seq.ov / ov.carbon.seq.tot)[[1]]

    o <- data.frame(
      project = paste0(project.name, "__", version),
      scenario = file.name,
      region = "Europe",
      variable = "Carbon stock conserved by the solution",
      variable_variant = "Overall value",
      unit = "tC",
      value = carbon.ov
    )
    o2 <- o %>% mutate(
      variable_variant = "Selected area in region",
      unit = "tC/ha",
      value = carbon.ov.norm
    )
    o3 <- o %>% mutate(
      variable_variant = "Total available in region",
      unit = "%",
      value = carbon.ov.prop * 100
    )
    s <- o %>% mutate(
      variable = "Carbon sequestration conserved by the solution",
      variable_variant = "Overall value",
      unit = "tC/year",
      value = carbon.seq.ov
    )
    s2 <- s %>% mutate(
      variable_variant = "Selected area in region",
      unit = "tC/ha/year",
      value = carbon.seq.ov.norm
    )
    s3 <- s %>% mutate(
      variable_variant = "Total available in region",
      unit = "%",
      value = carbon.seq.ov.prop * 100
    )

    carbon.europe.existing <- rbind(carbon.europe.existing, o, o2, o3, s, s2, s3)

    for (j in 1:nrow(bioregion)) {
      # loop through the bioregions
      region1 <- bioregion[j, ] %>% sf::st_make_valid()

      # calculate the area selected in the solution
      area.sel <- exactextractr::exact_extract(pa.sol.area, region1, fun = "sum")
      reg.carb.tot <- exactextractr::exact_extract(ov.carbon, region1, fun = "sum")
      reg.carb.seq.tot <- exactextractr::exact_extract(ov.carbon.seq, region1, fun = "sum")

      # now we extract the sum in tonnes of carbon conserved at the scale of the bioregion/country
      carb.out <- exactextractr::exact_extract(pa.sol.carbon, region1, fun = "sum")
      carb.seq.out <- exactextractr::exact_extract(pa.sol.carbon.seq, region1, fun = "sum")

      # and here we normalize by the area of the bioregion/country that was selected in the solution
      carb.out.norm <- (carb.out / as.numeric(area.sel))
      carb.sequest.norm <- (carb.seq.out / as.numeric(area.sel))
      carb.out.prop <- (carb.out / as.numeric(reg.carb.tot))
      carb.sequest.prop <- (carb.seq.out / as.numeric(reg.carb.seq.tot))

      o <- data.frame(
        project = paste0(project.name, "__", version),
        scenario = file.name,
        region = paste0(as.character(bioregion[[j, "code"]])),
        variable = "Carbon stock conserved by the solution",
        variable_variant = "Overall value",
        unit = "tC",
        value = carb.out
      )

      o2 <- o %>% mutate(
        variable_variant = "Selected area in region",
        unit = "tC/ha",
        value = carb.out.norm
      )

      o3 <- o %>% mutate(
        variable_variant = "Total available in region",
        unit = "%",
        value = carb.out.prop * 100
      )

      s <- o %>% mutate(
        variable = "Carbon sequestration conserved by the solution",
        variable_variant = "Overall value",
        unit = "tC/year",
        value = carb.seq.out
      )

      s2 <- s %>% mutate(
        variable_variant = "Selected area in region",
        unit = "tC/ha/year",
        value = carb.sequest.norm
      )

      s3 <- s %>% mutate(
        variable_variant = "Total available in region",
        unit = "%",
        value = carb.sequest.prop * 100
      )

      carbon.bioregion.existing <- rbind(carbon.bioregion.existing, o, s, o2, s2, o3, s3)
    }

    for (j in 1:nrow(countries)) {
      # loop through the countries
      region1 <- countries[j, ] %>% sf::st_make_valid()

      # calculate the area selected in the solution
      area.sel <- exactextractr::exact_extract(pa.sol.area, region1, fun = "sum")
      reg.carb.tot <- exactextractr::exact_extract(ov.carbon, region1, fun = "sum")
      reg.carb.seq.tot <- exactextractr::exact_extract(ov.carbon.seq, region1, fun = "sum")

      # now we run zonal stats to quantify the sum in tonnes of carbon conserved at the scale of the bioregion/country
      carb.out <- exactextractr::exact_extract(pa.sol.carbon, region1, fun = "sum")
      carb.seq.out <- exactextractr::exact_extract(pa.sol.carbon.seq, region1, fun = "sum")

      # and here we normalize by the area of the bioregion/country that was selected in the solution
      carb.out.norm <- (carb.out / as.numeric(area.sel))
      carb.sequest.norm <- (carb.seq.out / as.numeric(area.sel))
      carb.out.prop <- (carb.out / as.numeric(reg.carb.tot))
      carb.sequest.prop <- (carb.seq.out / as.numeric(reg.carb.seq.tot))

      o <- data.frame(
        project = paste0(project.name, "__", version),
        scenario = file.name,
        region = paste0(as.character(countries[[j, "code"]])),
        variable = "Carbon stock conserved by the solution",
        variable_variant = "Overall value",
        unit = "tC",
        value = carb.out
      )

      o2 <- o %>% mutate(
        variable_variant = "Selected area in region",
        unit = "tC/ha",
        value = carb.out.norm
      )

      o3 <- o %>% mutate(
        variable_variant = "Total available in region",
        unit = "%",
        value = carb.out.prop * 100
      )

      s <- o %>% mutate(
        variable = "Carbon sequestration conserved by the solution",
        variable_variant = "Overall value",
        unit = "tC/year",
        value = carb.seq.out
      )

      s2 <- s %>% mutate(
        variable_variant = "Selected area in region",
        unit = "tC/ha/year",
        value = carb.sequest.norm
      )

      s3 <- s %>% mutate(
        variable_variant = "Total available in region",
        unit = "%",
        value = carb.sequest.prop * 100
      )

      carbon.countries.existing <- rbind(carbon.countries.existing, o, s, o2, s2, o3, s3)
    }
  }
  toc()

  # Bind results for europe, bioregions and countries, then export to output directory
  out.carb.final <- rbind(carbon.europe.existing, carbon.bioregion.existing, carbon.countries.existing)

  outdir <- paste0("performance_outputs/", version)
  dir.create(outdir, showWarnings = FALSE)

  ofname <- paste0(outdir, "/", project.name, "__", version, "__indicator_conservedcarbon.csv")
  write.csv(out.carb.final, ofname, row.names = FALSE)

  message("Carbon stock and sequestration calculations done. Don't forget to create/update variable definitions")
}
