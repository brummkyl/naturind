#' Irreplaceability
#'
#' @param ll A character vector containing file paths of raster files (e.g., spatial prioritization outputs) in '.tif' format.
#' @param puid A SpatRaster layer with planning units
#' @param all.pa A SpatRaster layer with all existing protected areas
#' @param st.pa A SpatRaster layer with all existing strict protected areas
#' @param n2k A SpatRaster layer with all existing Natura 2000 protected areas (or any other unique subset of PAs)
#' @param countries An sf vector layer containing individual countries
#' @param bioregion An sf vector layer containing individual biogeographic regions
#' @param gamma.irr A SpatRaster layer containing gamma irreplaceability scores for planning units
#' @param irr.quantile A double specifying the gamma irreplaceability quantile
#' @param project.name A character vector specifying the name of your project
#' @param version A character vector specifying the version (e.g., date) of your workflow
#'
#' @returns A formatted .csv file saved to the "performance_outputs/version" directory
#' @export
#'
#' @examples
#' \dontrun{
#' pi.irreplaceability(ll = "C:/.../solution_file.tif", puid, all.pa, st.pa, n2k, countries,
#'                      bioregion, beta.irr, irr.quantile = 0.9, project.name = "Example Project",
#'                      version = "Version 1")
#'}

pi.irreplaceability <- function(ll, puid, all.pa, st.pa, n2k, countries, bioregion, gamma.irr, irr.quantile, project.name, version){

  beta.irr <- gamma.irr

  # Transform to ensure alignment of CRS
  if(terra::crs(puid) != terra::crs(beta.irr)){
    beta.irr <- beta.irr |> terra::project(puid)
  }
  assertthat::assert_that(crs(beta.irr) == crs(puid_area))

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

  puid <- puid %>% mask(countries)
  puid_area <- puid

  ee <- ll %>% rast() %>% extend(puid)
  existing_pa_names <- c(basename(ll))

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

  ee <- ee %>% mask(puid)
  #ee <- ee %>% mask(ee[[1]])

  beta.irr <- beta.irr %>% extend(puid) %>% crop(ee[[1]]) %>% mask(ee[[1]])

  # Resample to align resolution and origin of irreplaceability and puid_area rasters
  # Only resample if these two layers have different resolutions
  if(!all.equal(res(beta.irr), res(puid_area))) {
    beta.irr <- terra::resample(beta.irr,puid_area,method="bilinear")
  }
  assertthat::assert_that(all.equal(res(beta.irr), res(puid_area)), msg = "Rasters differ in resolution")

  # extract any given quantile from the beta irreplaceability layer, and then convert this to a binary format
  quant <- global(beta.irr, \(i) quantile(i, irr.quantile, na.rm=T))
  beta.irr.binary <- ifel(beta.irr <= quant[[1]], NA, 1)

  #######################################################################################

  # EXISTING PROTECTED AREAS

  #######################################################################################

  # bioregions and countries
  beta.irr.europe.existing <- data.frame()
  beta.irr.bioregion.existing <- data.frame()
  beta.irr.countries.existing <- data.frame()
  pb <- progress::progress_bar$new(total = nlyr(ee))

  tic()
  for (i in 1:nlyr(ee)) {
    pb$tick()
    # loop through the solution rasters
    pa.sol <- ee[[i]]
    file.name <- existing_pa_names[[i]]

    # calculate the selected area
    pa.sol.area <- pa.sol #* puid_area
    # calculate the amount of irreplaceable area
    betaarea <- puid_area * beta.irr.binary
    # calculate the selected amount of irreplaceable area
    pa.sol.beta <- pa.sol * beta.irr.binary #* puid_area

    # Overall sum
    ov.eur.area <- global(puid_area,"sum",na.rm=TRUE)[[1]]
    ov.sol.area <- global(pa.sol.area,"sum",na.rm=TRUE)[[1]]
    ov.beta.area <- global(betaarea,"sum",na.rm=TRUE)[[1]]
    ov.beta.out <- global(pa.sol.beta,"sum",na.rm=TRUE)[[1]]
    # normalization
    # normalize by the total area of each region
    beta.out.norm <- (ov.beta.out / as.numeric(ov.eur.area))
    # normalize by the selected area within each region
    beta.out.norm2 <- (ov.beta.out / as.numeric(ov.sol.area))
    # normalize by the amount of irreplaceable area within each region
    beta.out.norm3 <- (ov.beta.out / as.numeric(ov.beta.area))

    o <- data.frame(project = paste0(project.name,"__",version),
                    scenario = file.name,
                    region = "Europe",
                    variable = "Irreplaceability (gamma) selected in the solution",
                    variable_variant = "Total area of region",
                    unit = "%",
                    value = beta.out.norm)
    o2 <- o %>% mutate(
      variable_variant = "Selected area in region",
      unit = "%",
      value = beta.out.norm2)
    o3 <- o %>% mutate(
      variable_variant = "Irreplaceable (gamma) area in region",
      unit = "%",
      value = beta.out.norm3)

    beta.irr.europe.existing <- rbind(beta.irr.europe.existing, o, o2, o3)

    for (j in 1:nrow(bioregion)) {
      # loop through the bioregions
      region1 <- bioregion[j,] %>% sf::st_make_valid()

      # denominators for normalization
      # calculate the total area of each region, in km2
      region1.area <- exactextractr::exact_extract(puid_area, region1, fun="sum")
      # calculate the selected area within each region
      area.sel <- exactextractr::exact_extract(pa.sol.area, region1, fun="sum")
      # calculate the amount of irreplaceable area within each region
      betaarea.sel <- exactextractr::exact_extract(betaarea, region1, fun="sum")
      # numerator for normalization
      # calculate the selected amount of irreplaceable area within each region
      beta.out <- exactextractr::exact_extract(pa.sol.beta, region1, fun="sum")
      # normalization
      # normalize by the total area of each region
      beta.out.norm <- (beta.out / as.numeric(region1.area))
      # normalize by the selected area within each region
      beta.out.norm2 <- (beta.out / as.numeric(area.sel))
      # normalize by the amount of irreplaceable area within each region
      beta.out.norm3 <- (beta.out / as.numeric(betaarea.sel))

      o <- data.frame(project = paste0(project.name,"__",version),
                      scenario = file.name,
                      region = paste0(as.character(bioregion[[j,"code"]])),
                      variable = "Irreplaceability (gamma) selected in the solution",
                      variable_variant = "Total area of region",
                      unit = "%",
                      value = beta.out.norm)
      o2 <- o %>% mutate(
        variable_variant = "Selected area in region",
        unit = "%",
        value = beta.out.norm2)
      o3 <- o %>% mutate(
        variable_variant = "Irreplaceable (gamma) area in region",
        unit = "%",
        value = beta.out.norm3)

      beta.irr.bioregion.existing <- rbind(beta.irr.bioregion.existing, o, o2, o3)
    }

    for (j in 1:nrow(countries)) {
      # loop through the countries
      region1 <- countries[j,] %>% sf::st_make_valid()

      # denominators for normalization
      # calculate the total area of each region, in km2
      region1.area <- exactextractr::exact_extract(puid_area, region1, fun="sum")
      # calculate the selected area within each region
      area.sel <- exactextractr::exact_extract(pa.sol.area, region1, fun="sum")
      # calculate the amount of irreplaceable area within each region
      betaarea.sel <- exactextractr::exact_extract(betaarea, region1, fun="sum")
      # numerator for normalization
      # calculate the selected amount of irreplaceable area within each region
      beta.out <- exactextractr::exact_extract(pa.sol.beta, region1, fun="sum")
      # normalization
      # normalize by the total area of each region
      beta.out.norm <- (beta.out / as.numeric(region1.area))
      # normalize by the selected area within each region
      beta.out.norm2 <- (beta.out / as.numeric(area.sel))
      # normalize by the amount of irreplaceable area within each region
      beta.out.norm3 <- (beta.out / as.numeric(betaarea.sel))

      o <- data.frame(project = paste0(project.name,"__",version),
                      scenario = file.name,
                      region = paste0(as.character(countries[[j,"NAME_0"]])),
                      variable = "Irreplaceability (gamma) selected in the solution",
                      variable_variant = "Total area of region",
                      unit = "%",
                      value = beta.out.norm)
      o2 <- o %>% mutate(
        variable_variant = "Selected area in region",
        unit = "%",
        value = beta.out.norm2)
      o3 <- o %>% mutate(
        variable_variant = "Irreplaceable (gamma) area in region",
        unit = "%",
        value = beta.out.norm3)

      beta.irr.countries.existing <- rbind(beta.irr.countries.existing, o, o2, o3)
    }
  }
  toc()

  # Bind results for bioregions and countries, then export to output directory
  out.beta.irr.final <- rbind(beta.irr.europe.existing, beta.irr.bioregion.existing, beta.irr.countries.existing)
  # Convert from proportions to percentages
  out.beta.irr.final$value <- out.beta.irr.final$value * 100

  outdir <- paste0("performance_outputs/", version)
  dir.create(outdir, showWarnings = FALSE)

  #if(consensus.sol == 1){
  #  ofname <- paste0(outdir.consensus, "/", project.name,"__",version, "__indicator_irreplaceability.csv")
  #} else if (consensus.sol == 0){
    ofname <- paste0(outdir, "/", project.name,"__",version, "__indicator_irreplaceability.csv")
  #}

  write.csv(out.beta.irr.final, ofname, row.names = FALSE)

  message("Irreplaceability calculation done. Don't forget to create/update variable definitions")
}
