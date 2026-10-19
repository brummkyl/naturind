#' Protected Area Enlargement
#'
#' @param ll A character vector containing file paths of raster files (e.g., spatial prioritization outputs) in '.tif' format.
#' @param puid A SpatRaster layer with planning units
#' @param all.pa A SpatRaster layer with all existing protected areas
#' @param st.pa A SpatRaster layer with all existing strict protected areas
#' @param n2k A SpatRaster layer with all existing Natura 2000 protected areas (or any other unique subset of PAs)
#' @param variants A dataframe containing, at minimum, a column named "expand_on", which describes whether expansion priorities expand on
#' @param countries An sf vector layer containing individual countries
#' @param bioregion An sf vector layer containing individual biogeographic regions
#' @param pae_threshold A numeric specifying the distance threshold used to characterize expansion areas as enlargements or new establishments
#' @param project.name A character vector specifying the name of your project
#' @param version A character vector specifying the version (e.g., date) of your workflow
#'
#' @returns A formatted .csv file saved to the "performance_outputs/version" directory
#' @export
#'
#' @examples
#' \dontrun{
#' pi.enlargement(ll = "C:/.../solution_file.tif", puid, all.pa, st.pa, n2k, variants,
#'                 countries, bioregion, pae_threshold, project.name = "Example Project",
#'                 version = "Version 1")
#'}

pi.enlargement <- function(ll, puid, all.pa, st.pa, n2k, variants,
                           countries, bioregion, pae_threshold, project.name, version){


  # Script for calculating Queen's distance transition layers (i.e., how far away is each PU from the nearest PA?)

  # Create a raster of distances expressing location of PUs relative to existing PAs
  # NOTE: Distances only need to be recalculated when using a different or recently updated PA layer

  # Calculation of the raster quantifying distances from each PU to its nearest PA - all PAs
  {
    all.pa.bi <- ifel(all.pa > 0, 1, 0)
    # check that protected areas exist
    if (terra::global(all.pa.bi > 0, "sum", na.rm = TRUE)[1, 1] == 0) {
      stop("No PA cells found in `all.pa.bi`.")
    }

    # build an 8-neighbor TransitionLayer with a uniform movement cost (= cellsize)
    # this is the "queen" neighborhood approach (compared to directions = 4 as "rook" neighborhood)
    base_r <- raster::raster(puid)
    tr <- gdistance::transition(base_r, function(x) 1, directions = 8)
    tr <- gdistance::geoCorrection(tr, type = "c", multpl = FALSE)

    # extract coordinates of all protected area cells
    pa_cells <- which(terra::values(all.pa.bi) > 0)
    pa_xy <- raster::xyFromCell(base_r, pa_cells)

    # accumulated cost (Queen accumulated cost distance) to nearest PA
    # in map units (meters)
    manh_rl <- gdistance::accCost(tr, pa_xy)
    manh_sr <- terra::rast(manh_rl)
    manh_sr[!is.finite(manh_sr)] <- NA
    plot(manh_sr)

    # align CRS with that of planning units
    terra::crs(manh_sr) <- terra::crs(puid)
    assertthat::assert_that(terra::crs(puid) == terra::crs(manh_sr))
    #ofname <- paste0(path_input, "QueenDistances_allPA_",grain,".tif")
    #writeRaster(manh_sr, ofname, overwrite=T)
    mh.all <- manh_sr
    mh.strict.all <- manh_sr

  }

  # Calculation of the raster quantifying distances from each PU to its nearest PA - N2K only
  {
    n2k.bi <- ifel(n2k > 0, 1, 0)
    # check that protected areas exist
    if (terra::global(n2k.bi > 0, "sum", na.rm = TRUE)[1, 1] == 0) {
      stop("No PA cells found in `n2k.bi`.")
    }

    # build an 8-neighbor TransitionLayer with a uniform movement cost (= cellsize)
    # this is the "queen" neighborhood approach (compared to directions = 4 as "rook" neighborhood)
    base_r <- raster::raster(puid)
    tr <- gdistance::transition(base_r, function(x) 1, directions = 8)
    tr <- gdistance::geoCorrection(tr, type = "c", multpl = FALSE)

    # extract coordinates of all protected area cells
    pa_cells <- which(terra::values(n2k.bi) > 0)
    pa_xy <- raster::xyFromCell(base_r, pa_cells)

    # accumulated cost (Queen accumulated cost distance) to nearest PA
    # in map units (meters)
    manh_rl <- gdistance::accCost(tr, pa_xy)
    manh_sr <- terra::rast(manh_rl)
    manh_sr[!is.finite(manh_sr)] <- NA
    plot(manh_sr)

    # align CRS with that of planning units
    terra::crs(manh_sr) <- terra::crs(puid)
    assertthat::assert_that(terra::crs(puid) == terra::crs(manh_sr))
    #ofname <- paste0(path_input, "QueenDistances_N2K_",grain,".tif")
    #writeRaster(manh_sr, ofname, overwrite=T)
    mh.n2k <- manh_sr
    mh.strict.n2k <- manh_sr

  }

  # Calculation of the raster quantifying distances from each PU to its nearest PA - Strict-all PAs only
  {
    strict_all.bi <- ifel(st.pa > 0, 1, 0)
    # check that protected areas exist
    if (terra::global(strict_all.bi > 0, "sum", na.rm = TRUE)[1, 1] == 0) {
      stop("No PA cells found in `strict_all.bi`.")
    }

    # build an 8-neighbor TransitionLayer with a uniform movement cost (= cellsize)
    # this is the "queen" neighborhood approach (compared to directions = 4 as "rook" neighborhood)
    base_r <- raster::raster(puid)
    tr <- gdistance::transition(base_r, function(x) 1, directions = 8)
    tr <- gdistance::geoCorrection(tr, type = "c", multpl = FALSE)

    # extract coordinates of all protected area cells
    pa_cells <- which(terra::values(strict_all.bi) > 0)
    pa_xy <- raster::xyFromCell(base_r, pa_cells)

    # accumulated cost (Queen accumulated cost distance) to nearest PA
    # in map units (meters)
    manh_rl <- gdistance::accCost(tr, pa_xy)
    manh_sr <- terra::rast(manh_rl)
    manh_sr[!is.finite(manh_sr)] <- NA
    plot(manh_sr)

    # align CRS with that of planning units
    terra::crs(manh_sr) <- terra::crs(puid)
    assertthat::assert_that(terra::crs(puid) == terra::crs(manh_sr))
    #ofname <- paste0(path_input, "QueenDistances_Strict_all_",grain,".tif")
    #writeRaster(manh_sr, ofname, overwrite=T)
    #mh.strict.all <- manh_sr
    #mh.strict.n2k <- manh_sr
  }

  # # Calculation of the raster quantifying distances from each PU to its nearest PA - Strict-n2k PAs only
  # {
  #   strict_n2k.bi <- ifel(st.pa > 0, 1, 0)
  #   # check that protected areas exist
  #   if (terra::global(strict_n2k.bi > 0, "sum", na.rm = TRUE)[1, 1] == 0) {
  #     stop("No PA cells found in `strict_n2k.bi`.")
  #   }
  #
  #   # build an 8-neighbor TransitionLayer with a uniform movement cost (= cellsize)
  #   # this is the "queen" neighborhood approach (compared to directions = 4 as "rook" neighborhood)
  #   base_r <- raster::raster(puid)
  #   tr <- gdistance::transition(base_r, function(x) 1, directions = 8)
  #   tr <- gdistance::geoCorrection(tr, type = "c", multpl = FALSE)
  #
  #   # extract coordinates of all protected area cells
  #   pa_cells <- which(terra::values(strict_n2k.bi) > 0)
  #   pa_xy <- raster::xyFromCell(base_r, pa_cells)
  #
  #   # accumulated cost (Queen accumulated cost distance) to nearest PA
  #   # in map units (meters)
  #   manh_rl <- gdistance::accCost(tr, pa_xy)
  #   manh_sr <- terra::rast(manh_rl)
  #   manh_sr[!is.finite(manh_sr)] <- NA
  #   plot(manh_sr)
  #
  #   # align CRS with that of planning units
  #   crs(manh_sr) <- crs(puid)
  #   assertthat::assert_that(crs(puid) == crs(manh_sr))
  #   ofname <- paste0(path_input, "QueenDistances_Strict_n2k_",grain,".tif")
  #   writeRaster(manh_sr, ofname, overwrite=T)
  # }

  # Check CRS
  # # All PAs
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

  # Check CRS of distance rasters
  if(terra::crs(puid) != terra::crs(mh.all)){
    mh.all <- mh.all |> terra::project(puid)
  }
  assertthat::assert_that(crs(puid) == crs(mh.all))

  if(terra::crs(puid) != terra::crs(mh.n2k)){
    mh.n2k <- mh.n2k |> terra::project(puid)
  }
  assertthat::assert_that(crs(puid) == crs(mh.n2k))

  if(terra::crs(puid) != terra::crs(mh.strict.all)){
    mh.strict.all <- mh.strict.all |> terra::project(puid)
  }
  assertthat::assert_that(crs(puid) == crs(mh.strict.all))

  if(terra::crs(puid) != terra::crs(mh.strict.n2k)){
    mh.strict.n2k <- mh.strict.n2k |> terra::project(puid)
  }
  assertthat::assert_that(crs(puid) == crs(mh.strict.n2k))

  # Bioregion and Country Summaries
  out.europe <- data.frame()
  out.bioregion <- data.frame()
  out.countries <- data.frame()
  pb <- progress::progress_bar$new(total = length(ll))

  tic()
  # 1:length(ll)
  for (i in 1:length(ll)) {
    pb$tick()

    ll_base <- file_path_sans_ext(basename(ll[i]))  # strip extension only

    # find which basematch values are a prefix of ll_base
    is_match <- mapply(function(bm) {
      grepl(paste0("^", bm, "($|_)"), ll_base)
    }, variants$basematch)

    match_idx <- which(is_match)

    # here this is the full solution (expansion priorities are back-calculated below)
    pa.sol <- paste0(ll[[i]]) %>% terra::rast()
    # area selected in solution

    puid <- puid %>% mask(pa.sol)
    puid_area <- puid

    all.pa <- all.pa %>% mask(pa.sol)
    st.pa <- st.pa %>% mask(pa.sol)
    n2k <- n2k %>% mask(pa.sol)

    #all_base <- all_base %>% mask(pa.sol)
    mh.all <- mh.all %>% mask(pa.sol)
    #strict_all_base <- strict_all_base %>% mask(pa.sol)
    mh.strict.all <- mh.strict.all %>% mask(pa.sol)
    #strict_n2k_base <- strict_n2k_base %>% mask(pa.sol)
    mh.strict.n2k <- mh.strict.n2k %>% mask(pa.sol)
    #n2k_base <- n2k_base %>% mask(pa.sol)
    mh.n2k <- mh.n2k %>% mask(pa.sol)

    # SET PARAMETERS

    # if(consensus.sol == 1){
    #   if(variants$expand_on[match_idx] == "all"){ #grepl("conservation",basename(ll[[i]]))){
    #     pa.sol <- (pa.sol * puid_area) + all.pa
    #     pa.sol <- ifel(pa.sol > 1, 1, pa.sol)
    #     pa <- all.pa
    #     mh = mh.all
    #   } else if (grepl("strict",basename(ll[[i]]))){
    #     pa.sol <- (pa.sol * puid_area) + st.pa
    #     pa.sol <- ifel(pa.sol > 1, 1, pa.sol)
    #     pa <- st.pa
    #     mh <- mh.strict.all
    #   }
    #if(consensus.sol == 0){
      if (variants$expand_on[match_idx] == "all" && !grepl("strict", basename(ll[[i]]), ignore.case = TRUE)){ #grepl("all", basename(ll[[i]]), ignore.case = TRUE) && !grepl("strict", basename(ll[[i]]), ignore.case = TRUE)) {
        pa = all.pa
        mh = mh.all
        #var.variant = "All"
      } else if (variants$expand_on[match_idx] == "N2K" && !grepl("strict", basename(ll[[i]]), ignore.case = TRUE)){ #grepl("N2K", basename(ll[[i]]), ignore.case = TRUE) && !grepl("strict", basename(ll[[i]]), ignore.case = TRUE)) {
        pa = n2k
        mh = mh.n2k
        #var.variant = "N2K"
      } else if (variants$expand_on[match_idx] == "all" && grepl("strict", basename(ll[[i]]), ignore.case = TRUE)){ #grepl("all", basename(ll[[i]]), ignore.case = TRUE) && grepl("strict", basename(ll[[i]]), ignore.case = TRUE)) {
        pa = st.pa
        mh = mh.strict.all
        #var.variant = "Strict"
      } else if (variants$expand_on[match_idx] == "N2K" && grepl("strict", basename(ll[[i]]), ignore.case = TRUE)){ #grepl("N2K", basename(ll[[i]]), ignore.case = TRUE) && grepl("strict", basename(ll[[i]]), ignore.case = TRUE)) {
        pa = st.pa
        mh = mh.strict.n2k
        #var.variant = "Strict"
      }
    #}

    pa.sol.area <- pa.sol #* puid_area
    ov.sol.full <- ifel(pa.sol.area + pa > 1, 1, pa.sol.area + pa)
    ov.sol.area <- global(ov.sol.full,"sum",na.rm=TRUE)

    # GLOBAL CALCULATIONS

    # area of existing PAs
    out <- pa #* puid_area
    out.europe.sum <- global(out,"sum",na.rm=TRUE)[[1]]

    # 0 = (!selected AND !protected) | (proportion selected = proportion protected),
    # which in either case, does not represent an expansion or a new PA establishment.
    # values > 0 represent either expansion of existing PA | establishment of new PA
    sol <- pa.sol #- pa
    # reclassify to retain instances of PA expansion | establishment only
    sol <- ifel(sol > 0, sol, 0) #* puid_area
    # total expansion area for solution across all of Europe
    exp.sol.europe <- global(sol,"sum",na.rm=TRUE)[[1]]

    # we can use our Queen distance rasters to filter and differentiate between cases of expansion vs. establishment of a new PA
    # for this we rely on our pae threshold defined in parameter settings
    # here, values of 1 will count as an expansion
    # catching NAs here to ensure correct  of new PA establishments on isolated islands
    pae.expansion <- ifel(!is.na(mh) & mh <= pae_threshold, 1, 0)
    # here, values of 1 will count as a newly established PA
    pae.estab <- ifel(is.na(mh) | mh > pae_threshold, 1, 0)

    # mask to quantify selected pa expansions
    sol.expansion <- sol * pae.expansion
    # mask to quantify selected pa establishments
    sol.estab <- sol * pae.estab

    # normalized by proportion of total solution area, including area of existing PAs
    sol.expansion.europe <- (global(sol.expansion,"sum",na.rm=TRUE)/ov.sol.area[[1]])[[1]]
    sol.estab.europe <- (global(sol.estab,"sum",na.rm=TRUE)/ov.sol.area[[1]])[[1]]
    sol.exist.europe <- (global(out,"sum",na.rm=TRUE)/ov.sol.area[[1]])[[1]]

    # normalized by proportion of expansion area
    sol.expansion.europe.exp <- (global(sol.expansion,"sum",na.rm=TRUE)/exp.sol.europe)[[1]]
    sol.estab.europe.exp <- (global(sol.estab,"sum",na.rm=TRUE)/exp.sol.europe)[[1]]

    # Check that expansion, establishment, and existing add up to 100% of solution area
    total.fraction <- sum(c(sol.expansion.europe,
                            sol.estab.europe,
                            sol.exist.europe), na.rm = TRUE)
    assertthat::assert_that(abs(total.fraction - 1) < 1e-2)

    # and that expansion and establishment, when normalized by total expansion area, add up to 100% of the expansion area
    expansion.fraction <- sum(c(sol.expansion.europe.exp,
                                sol.estab.europe.exp), na.rm = TRUE)
    assertthat::assert_that(abs(expansion.fraction - 1) < 1e-2)

    # format the results
    eur.expansion <- data.frame(project = paste0(project.name,"__",version),
                                scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
                                region = "Europe",
                                variable = "Protected area expansion index",
                                variable_variant = "Selected area in region",
                                unit = "%",
                                value = sol.expansion.europe)

    eur.exp2 <- eur.expansion %>% mutate(variable_variant = "Expansion area in region",
                                         unit = "%",
                                         value = sol.expansion.europe.exp)

    eur.establishment <- data.frame(project = paste0(project.name,"__",version),
                                    scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
                                    region = "Europe",
                                    variable = "Protected area establishment index",
                                    variable_variant = "Selected area in region",
                                    unit = "%",
                                    value = sol.estab.europe)

    eur.est2 <- eur.establishment %>% mutate(variable_variant = "Expansion area in region",
                                             unit = "%",
                                             value = sol.estab.europe.exp)

    out.europe <- rbind(out.europe, eur.expansion, eur.establishment, eur.exp2, eur.est2)

    # BIOREGION SUMMARIES

    for (j in 1:nrow(bioregion)) {
      region1 <- bioregion[j,] %>% sf::st_make_valid()

      # calculate solution area for normalization
      sol.area.full <- ifel(pa.sol.area + pa > 1, 1, pa.sol.area + pa)

      sol.area <- exactextractr::exact_extract(sol.area.full, region1, fun="sum")
      exp.area <- exactextractr::exact_extract(sol, region1, fun="sum")

      if(sol.area > 0) {
        # expansions
        sol.expansion.region <- exactextractr::exact_extract(sol.expansion, region1, fun="sum") / sol.area
        sol.expansion.region.exp <- exactextractr::exact_extract(sol.expansion, region1, fun="sum") / exp.area

        # newly established PAs
        sol.estab.region <- exactextractr::exact_extract(sol.estab, region1, fun="sum") / sol.area
        sol.estab.region.exp <- exactextractr::exact_extract(sol.estab, region1, fun="sum") / exp.area

        # percent of current PAs
        out.region <- exactextractr::exact_extract(out, region1, fun="sum") / sol.area
        # check that the sum of PA expansion, PA establishments, and existing PAs is effectively 1
        # here we allow some room in case of rounding errors
        total.fraction <- sum(c(sol.expansion.region,
                                sol.estab.region,
                                out.region), na.rm = TRUE)
        assertthat::assert_that(abs(total.fraction - 1) < 1e-2)

        # and that expansion and establishment, when normalized by total expansion area, add up to 100% of the expansion area
        if(exp.area > 0){
          expansion.fraction <- sum(c(sol.expansion.region.exp,
                                      sol.estab.region.exp), na.rm = TRUE)
          assertthat::assert_that(abs(expansion.fraction - 1) < 1e-2)
        }
        else {
          sol.expansion.region <- NA
          sol.estab.region <- NA
          sol.expansion.region.exp <- NA
          sol.estab.region.exp <- NA
        }

      } else {
        sol.expansion.region <- NA
        sol.estab.region <- NA
        sol.expansion.region.exp <- NA
        sol.estab.region.exp <- NA
      }

      bio.expansion <- data.frame(project = paste0(project.name,"__",version),
                                  scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
                                  region = paste0(as.character(bioregion[[j,"code"]])),
                                  variable = "Protected area expansion index",
                                  variable_variant = "Selected area in region",
                                  unit = "%",
                                  value = sol.expansion.region)

      bio.exp2 <- bio.expansion %>% mutate(variable_variant = "Expansion area in region",
                                           unit = "%",
                                           value = sol.expansion.region.exp)

      bio.establishment <- data.frame(project = paste0(project.name,"__",version),
                                      scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
                                      region = paste0(as.character(bioregion[[j,"code"]])),
                                      variable = "Protected area establishment index",
                                      variable_variant = "Selected area in region",
                                      unit = "%",
                                      value = sol.estab.region)

      bio.est2 <- bio.establishment %>% mutate(variable_variant = "Expansion area in region",
                                               unit = "%",
                                               value = sol.estab.region.exp)

      out.bioregion <- rbind(out.bioregion, bio.expansion, bio.establishment, bio.exp2, bio.est2)
    }

    # COUNTRY SUMMARIES

    for (j in 1:nrow(countries)) {
      region1 <- countries[j,] %>% sf::st_make_valid()

      # calculate solution area for normalization
      sol.area.full <- ifel(pa.sol.area + pa > 1, 1, pa.sol.area + pa)

      sol.area <- exactextractr::exact_extract(sol.area.full, region1, fun="sum")
      exp.area <- exactextractr::exact_extract(sol, region1, fun="sum")

      if(sol.area > 0) {
        # expansions
        sol.expansion.region <- exactextractr::exact_extract(sol.expansion, region1, fun="sum") / sol.area
        sol.expansion.region.exp <- exactextractr::exact_extract(sol.expansion, region1, fun="sum") / exp.area

        # newly established PAs
        sol.estab.region <- exactextractr::exact_extract(sol.estab, region1, fun="sum") / sol.area
        sol.estab.region.exp <- exactextractr::exact_extract(sol.estab, region1, fun="sum") / exp.area

        # percent of current PAs
        out.region <- exactextractr::exact_extract(out, region1, fun="sum") / sol.area
        # check that the sum of PA expansion, PA establishments, and existing PAs is effectively 1
        # here we allow some room in case of rounding errors
        total.fraction <- sum(c(sol.expansion.region,
                                sol.estab.region,
                                out.region), na.rm = TRUE)
        assertthat::assert_that(abs(total.fraction - 1) < 1e-2)

        # and that expansion and establishment, when normalized by total expansion area, add up to 100% of the expansion area
        if(exp.area > 0){
          expansion.fraction <- sum(c(sol.expansion.region.exp,
                                      sol.estab.region.exp), na.rm = TRUE)
          assertthat::assert_that(abs(expansion.fraction - 1) < 1e-2)
        } else {
          sol.expansion.region <- NA
          sol.estab.region <- NA
          sol.expansion.region.exp <- NA
          sol.estab.region.exp <- NA
        }

      } else {
        sol.expansion.region <- NA
        sol.estab.region <- NA
        sol.expansion.region.exp <- NA
        sol.estab.region.exp <- NA
      }

      c.expansion <- data.frame(project = paste0(project.name,"__",version),
                                scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
                                region = paste0(as.character(countries[[j,"NAME_0"]])),
                                variable = "Protected area expansion index",
                                variable_variant = "Selected area in region",
                                unit = "%",
                                value = sol.expansion.region)

      c.exp2 <- c.expansion %>% mutate(variable_variant = "Expansion area in region",
                                       unit = "%",
                                       value = sol.expansion.region.exp)

      c.establishment <- data.frame(project = paste0(project.name,"__",version),
                                    scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
                                    region = paste0(as.character(countries[[j,"NAME_0"]])),
                                    variable = "Protected area establishment index",
                                    variable_variant = "Selected area in region",
                                    unit = "%",
                                    value = sol.estab.region)

      c.est2 <- c.establishment %>% mutate(variable_variant = "Expansion area in region",
                                           unit = "%",
                                           value = sol.estab.region.exp)

      out.countries <- rbind(out.countries, c.expansion, c.establishment, c.exp2, c.est2)
    }
  }
  toc()

  # Bind results for bioregions and countries, then export to output directory
  out.enlargement.final <- rbind(out.europe, out.bioregion, out.countries)
  # Convert from proportions to percentages
  out.enlargement.final$value <- out.enlargement.final$value * 100

  outdir <- paste0("performance_outputs/", version)
  dir.create(outdir, showWarnings = FALSE)

  #if(consensus.sol == 1){
  #  ofname <- paste0(outdir.consensus, "/", project.name,"__",version, "__indicator_paenlargement.csv")
  #} else if (consensus.sol == 0){
    ofname <- paste0(outdir, "/", project.name,"__",version, "__indicator_paenlargement.csv")
  #}

  write.csv(out.enlargement.final, ofname, row.names = FALSE)

  message("PA enlargement and establishment calculations done. Don't forget to create/update variable definitions")
}
