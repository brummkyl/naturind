pi.carbon <- function(ll, puid, all.pa, st.pa, n2k, countries, bioregion, carbon.stock = carbon, carbon.sequestration = carbon.seq, project.name, version){

  carbon <- carbon.stock
  carbon.seq <- carbon.sequestration

  ee <- ll %>% rast() %>% extend(puid)
  existing_pa_names <- c(basename(ll))

  # Transform to ensure alignment of CRS, and run assertions
  if(terra::crs(puid) != terra::crs(carbon)){
    carbon <- carbon |> terra::project(puid)
  }
  assertthat::assert_that(crs(carbon) == crs(puid))

  if(terra::crs(puid) != terra::crs(carbon.seq)){
    carbon.seq <- carbon.seq |> terra::project(puid)
  }
  assertthat::assert_that(crs(carbon.seq) == crs(puid))

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

  #ee <- ee %>% mask(ee[[1]])

  carbon <- carbon %>% mask(ee[[1]])
  carbon.seq <- carbon.seq %>% mask(ee[[1]])

  # Resample to align resolution and origin of carbon and puid_area_ha rasters
  # carbon storage
  if(any(terra::res(carbon) != terra::res(puid_area_ha))) {
    carbon <- terra::resample(carbon,puid_area_ha,method="bilinear")
  }
  assertthat::assert_that(all.equal(res(carbon), res(puid_area_ha)))
  # carbon sequestration
  if(any(terra::res(carbon.seq) != terra::res(puid_area_ha))) {
    carbon.seq <- terra::resample(carbon.seq,puid_area_ha,method="bilinear")
  }
  assertthat::assert_that(all.equal(res(carbon.seq), res(puid_area_ha)))

  # for this indicator, we only consider areas in which carbon sequestration is positive
  carbon.seq <- ifel(carbon.seq < 0, 0, carbon.seq)

  #######################################################################################

  # EXISTING PROTECTED AREAS

  #######################################################################################

  # bioregions and countries
  carbon.europe.existing <- data.frame()
  carbon.bioregion.existing <- data.frame()
  carbon.countries.existing <- data.frame()
  pb <- progress::progress_bar$new(total = nlyr(ee))

  tic()
  for (i in 1:nlyr(ee)) {
    pb$tick()
    # loop through the solution rasters
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
    ov.area <- global(pa.sol.area,"sum",na.rm=TRUE)[[1]]
    ov.carbon.tot <- global(ov.carbon,"sum",na.rm=TRUE)[[1]]
    ov.carbon.seq.tot <- global(ov.carbon.seq,"sum",na.rm=TRUE)[[1]]
    carbon.ov <- global(pa.sol.carbon,"sum",na.rm=TRUE)[[1]]
    carbon.seq.ov <- global(pa.sol.carbon.seq,"sum",na.rm=TRUE)[[1]]
    # Overall normalization
    carbon.ov.norm <- (carbon.ov/ov.area)[[1]]
    carbon.seq.ov.norm <- (carbon.seq.ov/ov.area)[[1]]
    carbon.ov.prop <- (carbon.ov/ov.carbon.tot)[[1]]
    carbon.seq.ov.prop <- (carbon.seq.ov/ov.carbon.seq.tot)[[1]]

    o <- data.frame(project = paste0(project.name,"__",version),
                    scenario = file.name,
                    region = "Europe",
                    variable = "Carbon stock conserved by the solution",
                    variable_variant = "Overall value",
                    unit = "tC",
                    value = carbon.ov)
    o2 <- o %>% mutate(
      variable_variant = "Selected area in region",
      unit = "tC/ha",
      value = carbon.ov.norm)
    o3 <- o %>% mutate(
      variable_variant = "Total available in region",
      unit = "%",
      value = carbon.ov.prop * 100)
    s <- o %>% mutate(
      variable = "Carbon sequestration conserved by the solution",
      variable_variant = "Overall value",
      unit = "tC/year",
      value = carbon.seq.ov)
    s2 <- s %>% mutate(
      variable_variant = "Selected area in region",
      unit = "tC/ha/year",
      value = carbon.seq.ov.norm)
    s3 <- s %>% mutate(
      variable_variant = "Total available in region",
      unit = "%",
      value = carbon.seq.ov.prop * 100)

    carbon.europe.existing <- rbind(carbon.europe.existing, o, o2, o3, s, s2, s3)

    for (j in 1:nrow(bioregion)) {
      # loop through the bioregions
      region1 <- bioregion[j,] %>% sf::st_make_valid()

      # calculate the area selected in the solution
      area.sel <- exactextractr::exact_extract(pa.sol.area, region1, fun="sum")
      reg.carb.tot <- exactextractr::exact_extract(ov.carbon, region1, fun="sum")
      reg.carb.seq.tot <- exactextractr::exact_extract(ov.carbon.seq, region1, fun="sum")
      # now we run zonal stats to quantify the sum in tonnes of carbon conserved at the scale of the bioregion/country
      carb.out <- exactextractr::exact_extract(pa.sol.carbon, region1, fun="sum")
      carb.seq.out <- exactextractr::exact_extract(pa.sol.carbon.seq, region1, fun="sum")
      # and here we normalize by the area of the bioregion/country that was selected in the solution
      carb.out.norm <- (carb.out / as.numeric(area.sel))
      carb.sequest.norm <- (carb.seq.out / as.numeric(area.sel))
      carb.out.prop <- (carb.out / as.numeric(reg.carb.tot))
      carb.sequest.prop <- (carb.seq.out / as.numeric(reg.carb.seq.tot))

      o <- data.frame(project = paste0(project.name,"__",version),
                      scenario = file.name,
                      region = paste0(as.character(bioregion[[j,"code"]])),
                      variable = "Carbon stock conserved by the solution",
                      variable_variant = "Overall value",
                      unit = "tC",
                      value = carb.out)
      o2 <- o %>% mutate(
        variable_variant = "Selected area in region",
        unit = "tC/ha",
        value = carb.out.norm)
      o3 <- o %>% mutate(
        variable_variant = "Total available in region",
        unit = "%",
        value = carb.out.prop * 100)
      s <- o %>% mutate(
        variable = "Carbon sequestration conserved by the solution",
        variable_variant = "Overall value",
        unit = "tC/year",
        value = carb.seq.out)
      s2 <- s %>% mutate(
        variable_variant = "Selected area in region",
        unit = "tC/ha/year",
        value = carb.sequest.norm)
      s3 <- s %>% mutate(
        variable_variant = "Total available in region",
        unit = "%",
        value = carb.sequest.prop * 100)

      carbon.bioregion.existing <- rbind(carbon.bioregion.existing, o, s, o2, s2, o3, s3)
    }

    for (j in 1:nrow(countries)) {
      # loop through the bioregions
      region1 <- countries[j,] %>% sf::st_make_valid()

      # calculate the area selected in the solution
      area.sel <- exactextractr::exact_extract(pa.sol.area, region1, fun="sum")
      reg.carb.tot <- exactextractr::exact_extract(ov.carbon, region1, fun="sum")
      reg.carb.seq.tot <- exactextractr::exact_extract(ov.carbon.seq, region1, fun="sum")
      # now we run zonal stats to quantify the sum in tonnes of carbon conserved at the scale of the bioregion/country
      carb.out <- exactextractr::exact_extract(pa.sol.carbon, region1, fun="sum")
      carb.seq.out <- exactextractr::exact_extract(pa.sol.carbon.seq, region1, fun="sum")
      # and here we normalize by the area of the bioregion/country that was selected in the solution
      carb.out.norm <- (carb.out / as.numeric(area.sel))
      carb.sequest.norm <- (carb.seq.out / as.numeric(area.sel))
      carb.out.prop <- (carb.out / as.numeric(reg.carb.tot))
      carb.sequest.prop <- (carb.seq.out / as.numeric(reg.carb.seq.tot))

      o <- data.frame(project = paste0(project.name,"__",version),
                      scenario = file.name,
                      region = paste0(as.character(countries[[j,"NAME_0"]])),
                      variable = "Carbon stock conserved by the solution",
                      variable_variant = "Overall value",
                      unit = "tC",
                      value = carb.out)
      o2 <- o %>% mutate(
        variable_variant = "Selected area in region",
        unit = "tC/ha",
        value = carb.out.norm)
      o3 <- o %>% mutate(
        variable_variant = "Total available in region",
        unit = "%",
        value = carb.out.prop * 100)
      s <- o %>% mutate(
        variable = "Carbon sequestration conserved by the solution",
        variable_variant = "Overall value",
        unit = "tC/year",
        value = carb.seq.out)
      s2 <- s %>% mutate(
        variable_variant = "Selected area in region",
        unit = "tC/ha/year",
        value = carb.sequest.norm)
      s3 <- s %>% mutate(
        variable_variant = "Total available in region",
        unit = "%",
        value = carb.sequest.prop * 100)

      carbon.countries.existing <- rbind(carbon.countries.existing, o, s, o2, s2, o3, s3)
    }
  }
  toc()

  # Bind results for europe, bioregions and countries, then export to output directory
  out.carb.final <- rbind(carbon.europe.existing, carbon.bioregion.existing, carbon.countries.existing)

  if(consensus.sol == 1){
    ofname <- paste0(outdir.consensus, "/", project.name,"__",version, "__indicator_conservedcarbon.csv")
  } else if (consensus.sol == 0){
    ofname <- paste0(outdir, "/", project.name,"__",version, "__indicator_conservedcarbon.csv")
  }

  write.csv(out.carb.final, ofname, row.names = FALSE)

  message("Carbon stock and sequestration calculations done. Don't forget to create/update variable definitions")
}
