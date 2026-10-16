pi.mean.area <- function(ll, puid, all.pa, st.pa, n2k, countries, bioregion, project.name, version){

  ee <- ll %>% rast() %>% extend(puid)
  existing_pa_names <- c(basename(ll))

  puid_area <- puid

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

  ee <- c(ee, all.pa, st.pa, n2k)
  existing_pa_names <- c(existing_pa_names, "Existing all PA", "Existing strict PA", "Existing N2K")

  puid <- puid %>% mask(ee[[1]])
  puid_area <- puid_area %>% mask(ee[[1]])

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

  #######################################################################################

  # EXISTING PROTECTED AREAS

  #######################################################################################

  # bioregions and countries
  area.europe.existing <- data.frame()
  area.bioregion.existing <- data.frame()
  area.countries.existing <- data.frame()
  pb <- progress::progress_bar$new(total = nlyr(ee))

  tic()
  for (i in 1:nlyr(ee)) {
    pb$tick()
    # loop through the solution rasters
    pa.sol <- ee[[i]]
    file.name <- existing_pa_names[[i]]

    # calculate the selected area; for existing pa, do not multiply by land area
    # this step is only needed for proportional solutions
    pa.sol.area <- pa.sol

    # Overall sum
    ov.eur.area <- global(puid_area,"sum",na.rm=TRUE)[[1]]
    ov.sol.area <- global(pa.sol.area,"sum",na.rm=TRUE)[[1]]

    # normalization
    # normalize by the total area of each region
    area.out.norm <- (ov.sol.area / as.numeric(ov.eur.area))

    o <- data.frame(project = paste0(project.name,"__",version),
                    scenario = file.name,
                    region = "Europe",
                    variable = "Relative amount of area covered",
                    variable_variant = "Total area of region",
                    unit = "%",
                    value = area.out.norm)

    area.europe.existing <- rbind(area.europe.existing, o)

    for (j in 1:nrow(bioregion)) {
      # loop through the bioregions
      # TODO: See if we need to be clipping and masking region 1 to puid_area
      region1 <- bioregion[j,] %>% sf::st_make_valid()

      # denominators for normalization
      # calculate the total area of each region, in km2
      region1.area <- exactextractr::exact_extract(puid_area, region1, fun="sum")
      # calculate the selected area within each region
      area.sel <- exactextractr::exact_extract(pa.sol.area, region1, fun="sum")

      # normalization
      # normalize by the total area of each region
      area.out.norm <- (area.sel / as.numeric(region1.area))
      # normalize by the selected area within each region
      area.out.norm2 <- (area.sel / as.numeric(ov.eur.area))

      o <- data.frame(project = paste0(project.name,"__",version),
                      scenario = file.name,
                      region = paste0(as.character(bioregion[[j,"code"]])),
                      variable = "Relative amount of area covered",
                      variable_variant = "Total area of region",
                      unit = "%",
                      value = area.out.norm)

      o2 <- o %>% mutate(
        variable_variant = "Total area of Europe",
        unit = "%",
        value = area.out.norm2)

      area.bioregion.existing <- rbind(area.bioregion.existing, o, o2)
    }

    for (j in 1:nrow(countries)) {
      # loop through the countries
      # TODO: See if we need to be clipping and masking region 1 to puid_area
      region1 <- countries[j,] %>% sf::st_make_valid()

      # denominators for normalization
      # calculate the total area of each region, in km2
      region1.area <- exactextractr::exact_extract(puid_area, region1, fun="sum")
      # calculate the selected area within each region
      area.sel <- exactextractr::exact_extract(pa.sol.area, region1, fun="sum")

      # normalization
      # normalize by the total area of each region
      area.out.norm <- (area.sel / as.numeric(region1.area))
      # normalize by the selected area within each region
      area.out.norm2 <- (area.sel / as.numeric(ov.eur.area))

      o <- data.frame(project = paste0(project.name,"__",version),
                      scenario = file.name,
                      region = paste0(as.character(countries[[j,"NAME_0"]])),
                      variable = "Relative amount of area covered",
                      variable_variant = "Total area of region",
                      unit = "%",
                      value = area.out.norm)

      o2 <- o %>% mutate(
        variable_variant = "Total area of Europe",
        unit = "%",
        value = area.out.norm2)

      area.countries.existing <- rbind(area.countries.existing, o, o2)
    }
  }
  toc()

  # Bind results for bioregions and countries, then export to output directory
  out.area.final <- rbind(area.europe.existing, area.bioregion.existing, area.countries.existing)
  # Convert from proportions to percentages
  out.area.final$value <- out.area.final$value * 100

  if(consensus.sol == 1){
    ofname <- paste0(outdir.consensus, "/", project.name,"__",version, "__indicator_relativeareacoverage.csv")
  } else if (consensus.sol == 0){
    ofname <- paste0(outdir, "/", project.name,"__",version, "__indicator_relativeareacoverage.csv")
  }

  write.csv(out.area.final, ofname, row.names = FALSE)
  message("Area calculation done. Don't forget to create/update variable definitions")
}
