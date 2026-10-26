#' Opportunity Costs
#'
#' @param ll A character vector containing file paths of raster files (e.g., expansion priorities) in '.tif' format.
#' @param puid A SpatRaster layer with planning units
#' @param countries An sf vector layer containing individual countries
#' @param bioregion An sf vector layer containing individual biogeographic regions
#' @param for.opp A SpatRaster layer with opportunity costs for forestry
#' @param ag.opp A SpatRaster layer with opportunity costs for agriculture
#' @param project.name A character vector specifying the name of your project
#' @param version A character vector specifying the version (e.g., date) of your workflow
#'
#' @returns A formatted .csv file saved to the "performance_outputs/version" directory
#' @export
#'
#' @examples
#' \dontrun{
#' pi.oppcosts(ll = "C:/.../solution_file.tif", puid, countries, bioregion,
#'              project.name = "Example Project", version = "Version 1")
#'}

pi.oppcosts <- function(ll, puid, countries, bioregion, for.opp, ag.opp, project.name, version){

  #strict <- st.pa

  # Transform to ensure alignment of CRS, and run assertions
  if(terra::crs(puid) != terra::crs(for.opp)){
    for.opp <- for.opp |> terra::project(puid, method = "bilinear")
  }
  assertthat::assert_that(crs(for.opp) == crs(puid))

  if(terra::crs(puid) != terra::crs(ag.opp)){
    ag.opp <- ag.opp |> terra::project(puid, method = "bilinear")
  }
  assertthat::assert_that(crs(ag.opp) == crs(puid))

  if(terra::crs(puid) != terra::crs(bioregion)){
    bioregion <- bioregion |> sf::st_transform(crs = sf::st_crs(puid))
  }
  assertthat::assert_that(sf::st_crs(bioregion) == sf::st_crs(puid))

  if(terra::crs(puid) != terra::crs(countries)){
    countries <- countries |> sf::st_transform(crs = sf::st_crs(puid))
  }
  assertthat::assert_that(sf::st_crs(countries) == sf::st_crs(puid))

  # if(terra::crs(puid) != terra::crs(all.pa)){
  #   all.pa <- all.pa |> terra::project(puid)
  # }
  # assertthat::assert_that(crs(puid) == crs(all.pa))
  #
  # # N2K
  # if(terra::crs(puid) != terra::crs(n2k)){
  #   n2k <- n2k |> terra::project(puid)
  # }
  # assertthat::assert_that(crs(puid) == crs(n2k))
  #
  # # Strict PAs
  # if(terra::crs(puid) != terra::crs(strict)){
  #   strict <- strict |> terra::project(puid)
  # }
  # assertthat::assert_that(crs(puid) == crs(strict))

  # Resample to align resolution and origin of opportunity cost and puid layers

  # forest opportunity costs
  for.opp_res <- terra::resample(for.opp,puid,method="bilinear")
  assertthat::assert_that(all.equal(res(for.opp_res), res(puid)))

  # ag opportunity costs
  ag.opp_res <- terra::resample(ag.opp,puid,method="bilinear")
  assertthat::assert_that(all.equal(res(ag.opp_res), res(puid)))

  co.dissolve <- sf::st_union(countries)

  # Check CRS
  # All PAs
  # if(terra::crs(puid) != terra::crs(all_base)){
  #   all_base <- all_base |> terra::project(puid)
  # }
  # assertthat::assert_that(crs(puid) == crs(all_base))
  # # N2K
  # if(terra::crs(puid) != terra::crs(n2k_base)){
  #   n2k_base <- n2k_base |> terra::project(puid)
  # }
  # assertthat::assert_that(crs(puid) == crs(n2k_base))
  # # Strict all PAs
  # if(terra::crs(puid) != terra::crs(strict_all_base)){
  #   strict_all_base <- strict_all_base |> terra::project(puid)
  # }
  # assertthat::assert_that(crs(puid) == crs(strict_all_base))
  # # Strict N2K PAs
  # if(terra::crs(puid) != terra::crs(strict_n2k_base)){
  #   strict_n2k_base <- strict_n2k_base |> terra::project(puid)
  # }
  # assertthat::assert_that(crs(puid) == crs(strict_n2k_base))

  # calculations

  oppcosts.europe <- data.frame()
  oppcosts.bioregion <- data.frame()
  oppcosts.countries <- data.frame()
  pb <- progress::progress_bar$new(total = length(ll))

  tic()
  for (i in 1:length(ll)) {
    pb$tick()
    # loop through the solution rasters
    pa.sol <- paste0(ll[[i]]) %>% terra::rast()

    puid <- puid %>% mask(pa.sol)
    puid_area <- puid

    # all.pa <- all.pa %>% mask(pa.sol)
    # st.pa <- st.pa %>% mask(pa.sol)
    # n2k <- n2k %>% mask(pa.sol)

    # all_base <- all_base %>% mask(pa.sol)
    # strict_all_base <- strict_all_base %>% mask(pa.sol)
    # strict_n2k_base <- strict_n2k_base %>% mask(pa.sol)
    # n2k_base <- n2k_base %>% mask(pa.sol)

    for.opp_res <- for.opp_res %>% mask(pa.sol)
    ag.opp_res <- ag.opp_res %>% mask(pa.sol)

    # # SET PARAMETERS
    # if(consensus.sol == 1){
    #   if(grepl("conservation",basename(ll[[i]]))){
    #     pa.sol <- (pa.sol * puid_area) + all.pa
    #     pa.sol <- ifel(pa.sol > 1, 1, pa.sol)
    #     pa <- all_base
    #     var.variant = "All"
    #   } else if (grepl("strict",basename(ll[[i]]))){
    #     pa.sol <- (pa.sol * puid_area) + strict
    #     pa.sol <- ifel(pa.sol > 1, 1, pa.sol)
    #     pa <- strict_all_base
    #     var.variant = "All"
    #   }
    # }  else if(consensus.sol == 0){
    #   if (grepl("all", basename(ll[[i]]), ignore.case = TRUE) && !grepl("strict", basename(ll[[i]]), ignore.case = TRUE)) {
    #     pa = all_base
    #     var.variant = "All"
    #   } else if (grepl("N2K", basename(ll[[i]]), ignore.case = TRUE) && !grepl("strict", basename(ll[[i]]), ignore.case = TRUE)) {
    #     pa = n2k_base
    #     var.variant = "N2K"
    #   } else if (grepl("all", basename(ll[[i]]), ignore.case = TRUE) && grepl("strict", basename(ll[[i]]), ignore.case = TRUE)) {
    #     pa = strict_all_base
    #     var.variant = "Strict"
    #   } else if (grepl("N2K", basename(ll[[i]]), ignore.case = TRUE) && grepl("strict", basename(ll[[i]]), ignore.case = TRUE)) {
    #     pa = strict_n2k_base
    #     var.variant = "Strict"
    #   }
    # }

    # reduce to expansion priorities only
    pa.sol.exp <- pa.sol #- pa
    # multiply the proportional solution by the planning unit area by the opportunity cost raster
    # output quantifies opportunity cost (in euros) for areas selected in the solution

    pa.sol.exp.area <- pa.sol.exp * 100 #* puid_area_ha
    pa.exp.for.opp <- pa.sol.exp.area * for.opp_res
    pa.exp.ag.opp <- pa.sol.exp.area * ag.opp_res

    # Overall calculations for Europe
    pa.exp.europe <- exactextractr::exact_extract(pa.sol.exp.area, co.dissolve, fun="sum")
    for.oppcosts.europe <- exactextractr::exact_extract(pa.exp.for.opp, co.dissolve, fun="sum")
    ag.oppcosts.europe <- exactextractr::exact_extract(pa.exp.ag.opp, co.dissolve, fun="sum")

    # pa.exp.europe <- global(pa.sol.exp.area,"sum",na.rm=TRUE)
    # for.oppcosts.europe <- global(pa.exp.for.opp,"sum",na.rm=TRUE)[[1]]
    # ag.oppcosts.europe <- global(pa.exp.ag.opp,"sum",na.rm=TRUE)[[1]]

    # Overall - normalization

    for.oppcosts.eur.norm <- (for.oppcosts.europe / as.numeric(pa.exp.europe))[[1]]
    ag.oppcosts.eur.norm <- (ag.oppcosts.europe / as.numeric(pa.exp.europe))[[1]]

    o <- data.frame(project = paste0(project.name,"__",version),
                    scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
                    region = "Europe",
                    variable = "Foregone opportunity costs (Forestry)",
                    variable_variant = "Total value",
                    unit = "EUR",
                    value = max(for.oppcosts.europe,0))
    o2 <- o %>% mutate(
      variable_variant = "Selected expansion area in region",
      unit = "EUR/ha",
      value = max(for.oppcosts.eur.norm,0))
    s <- o %>% mutate(
      variable = "Foregone opportunity costs (Agriculture)",
      variable_variant = "Total value",
      unit = "EUR",
      value = max(ag.oppcosts.europe,0))
    s2 <- s %>% mutate(
      variable_variant = "Selected expansion area in region",
      unit = "EUR/ha",
      value = max(ag.oppcosts.eur.norm,0))

    oppcosts.europe <- rbind(oppcosts.europe, o, s, o2, s2)

    # BIOREGIONS

    for (j in 1:nrow(bioregion)) {
      # loop through the bioregions
      region1 <- bioregion[j,] %>% sf::st_make_valid()

      # calculate the area of expansion priorities in the solution
      area.sel <- exactextractr::exact_extract(pa.sol.exp.area, region1, fun="sum")

      # now we run zonal stats to quantify the sum of opportunity costs,
      # within expansion priorities only, per bioregion/country
      for.oppcosts <- exactextractr::exact_extract(pa.exp.for.opp, region1, fun="sum")
      ag.oppcosts <- exactextractr::exact_extract(pa.exp.ag.opp, region1, fun="sum")

      # and here we normalize by the expansion area of the bioregion/country
      # that was selected in the solution
      for.oppcosts.norm <- (for.oppcosts / as.numeric(area.sel))
      ag.oppcosts.norm <- (ag.oppcosts / as.numeric(area.sel))

      o <- data.frame(project = paste0(project.name,"__",version),
                      scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
                      region = paste0(as.character(bioregion[[j,"code"]])),
                      variable = "Foregone opportunity costs (Forestry)",
                      variable_variant = "Total value",
                      unit = "EUR",
                      value = max(for.oppcosts,0))
      o2 <- o %>% mutate(
        variable_variant = "Selected expansion area in region",
        unit = "EUR/ha",
        value = max(for.oppcosts.norm,0))
      s <- o %>% mutate(
        variable = "Foregone opportunity costs (Agriculture)",
        variable_variant = "Total value",
        unit = "EUR",
        value = max(ag.oppcosts,0))
      s2 <- s %>% mutate(
        variable_variant = "Selected expansion area in region",
        unit = "EUR/ha",
        value = max(ag.oppcosts.norm,0))

      oppcosts.bioregion <- rbind(oppcosts.bioregion, o, s, o2, s2)
    }

    # COUNTRIES

    for (j in 1:nrow(countries)) {
      # loop through the bioregions
      region1 <- countries[j,] %>% sf::st_make_valid()

      # calculate the area of expansion priorities in the solution
      area.sel <- exactextractr::exact_extract(pa.sol.exp.area, region1, fun="sum")

      # now we run zonal stats to quantify the sum of opportunity costs,
      # within expansion priorities only, per bioregion/country
      for.oppcosts <- exactextractr::exact_extract(pa.exp.for.opp, region1, fun="sum")
      ag.oppcosts <- exactextractr::exact_extract(pa.exp.ag.opp, region1, fun="sum")

      # and here we normalize by the expansion area of the bioregion/country
      # that was selected in the solution
      for.oppcosts.norm <- (for.oppcosts / as.numeric(area.sel))
      ag.oppcosts.norm <- (ag.oppcosts / as.numeric(area.sel))

      o <- data.frame(project = paste0(project.name,"__",version),
                      scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
                      region = paste0(as.character(countries[[j,"NAME_0"]])),
                      variable = "Foregone opportunity costs (Forestry)",
                      variable_variant = "Total value",
                      unit = "EUR",
                      value = max(for.oppcosts,0))
      o2 <- o %>% mutate(
        variable_variant = "Selected expansion area in region",
        unit = "EUR/ha",
        value = max(for.oppcosts.norm,0))
      s <- o %>% mutate(
        variable = "Foregone opportunity costs (Agriculture)",
        variable_variant = "Total value",
        unit = "EUR",
        value = max(ag.oppcosts,0))
      s2 <- s %>% mutate(
        variable_variant = "Selected expansion area in region",
        unit = "EUR/ha",
        value = max(ag.oppcosts.norm,0))

      oppcosts.countries <- rbind(oppcosts.countries, o, s, o2, s2)
    }
  }
  toc()

  # Bind results for bioregions and countries, then export to output directory
  out.oppcosts.final <- rbind(oppcosts.europe, oppcosts.bioregion, oppcosts.countries)

  outdir <- paste0("performance_outputs/", version)
  dir.create(outdir, showWarnings = FALSE)

  #if(consensus.sol == 1){
  #  ofname <- paste0(outdir.consensus, "/", project.name,"__",version, "__indicator_foregoneopportunitycosts.csv")
  #} else if (consensus.sol == 0){
    ofname <- paste0(outdir, "/", project.name,"__",version, "__indicator_foregoneopportunitycosts.csv")
  #}

  write.csv(out.oppcosts.final, ofname, row.names = FALSE)

  message("Foregone opportunity cost calculations done. Don't forget to create/update variable definitions")
}
