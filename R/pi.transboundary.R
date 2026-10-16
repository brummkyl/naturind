pi.transboundary <- function(ll, puid, all.pa, st.pa, n2k, countries, bioregion, tb.buffers, project.name, version){

  ee <- ll %>% rast() %>% extend(puid)
  existing_pa_names <- c(basename(ll))

  # Transform to ensure alignment of CRS, and run assertions
  if(terra::crs(puid) != st_crs(tb.buffers)){
    tb.buffers <- tb.buffers |> sf::st_transform(crs = st_crs(puid))
  }
  assertthat::assert_that(st_crs(tb.buffers) == st_crs(puid))

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
    function(i) {identical(reference_na, is.na(values(ee[[i]])))},
    logical(1)
  ))

  # Only mask the stack if at least one layer differs
  if (!na_match) {
    message("NA footprints differ between ee layers; masking stack to align with the first input layer...")
    ee <- ee %>% mask(ee[[1]])
  }

  #######################################################################################

  # EXISTING PROTECTED AREAS

  #######################################################################################

  tb.europe.existing <- data.frame()
  pb <- progress::progress_bar$new(total = nlyr(ee))

  tic()
  for (i in 1:nlyr(ee)) {
    pb$tick()
    # loop through the solution rasters
    pa.sol <- ee[[i]]
    file.name <- existing_pa_names[[i]]
    # multiply the proportional solution by the planning unit area by the carbon stock raster
    # output quantifies tonnes of carbon per cell conserved in the solution

    # Area of protected areas
    pa.sol.area <- pa.sol
    # Total area of protected areas contained in the solution
    ov.area <- global(pa.sol.area, "sum", na.rm=T)[[1]]
    # Total area of transboundary buffers
    tb.area <- as.numeric(st_area(tb.buffers) / 1e6)
    # Total area of protected areas contained within transboundary country buffers
    buffer.exp <- exactextractr::exact_extract(pa.sol.area, tb.buffers, fun="sum")[[1]]
    # what proportion of protected areas are inside of transboundary regions?
    buffer.prop <- (buffer.exp / ov.area)[[1]]
    # what proportion of transboundary buffer area is protected
    buffer.prop2 <- as.numeric((buffer.exp / tb.area)[[1]])

    o <- data.frame(project = paste0(project.name,"__",version),
                    scenario = file.name,
                    region = "Europe",
                    variable = "Transboundary area conserved by the solution",
                    variable_variant = "Selected area in region",
                    unit = "%",
                    value = buffer.prop * 100)
    o2 <- o %>% mutate(
      variable_variant = "Transboundary area in region",
      unit = "%",
      value = buffer.prop2 * 100)

    tb.europe.existing <- rbind(tb.europe.existing, o, o2)
  }
  toc()

  # Bind results for europe - existing PAs and europe - full NaturaConnect solutions, then export to output directory
  out.tb.final <- rbind(tb.europe.existing)

  if(consensus.sol == 1){
    ofname <- paste0(outdir.consensus, "/", project.name,"__",version, "__indicator_transboundaryprotection.csv")
  } else if (consensus.sol == 0){
    ofname <- paste0(outdir, "/", project.name,"__",version, "__indicator_transboundaryprotection.csv")
  }

  write.csv(out.tb.final, ofname, row.names = FALSE)

  message("Transboundary protected area calculations done. Don't forget to create/update variable definitions")
}
