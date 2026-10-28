#' Protected Area Enlargement
#'
#' @param ll A character vector containing file paths of raster files (e.g., spatial prioritization outputs) in '.tif' format. This layer should contain the full solution (expansion priorities + existing protected areas), not just the expansion priorities.
#' @param puid A SpatRaster layer, reflecting the amount of land area per planning unit
#' @param ref.pa A SpatRaster layer with existing protected areas to be used as the baseline for calculating expansion areas
#' @param st.pa A SpatRaster layer with all existing strict protected areas
#' @param countries An sf vector layer containing individual countries, with names specified in a column named "code"
#' @param bioregion An sf vector layer containing individual biogeographic regions, with names specified in a column named "code"
#' @param ref.pa.name A character string specifying which protected area raster should be used as the baseline reference layer when calculating Queen's distances. For example, are you expanding on "all.pa" or "n2k"?
#' @param strict A logical, specifying whether your ll files represent strict expansion priorities
#' @param pae_threshold An integer specifying the distance threshold (in meters) used to characterize expansion areas as enlargements or new establishments
#' @param project.name A character vector specifying the name of your project
#' @param version A character vector specifying the version (e.g., date) of your workflow
#'
#' @returns A formatted .csv file saved to the "performance_outputs/version" directory
#' @export
#'
#' @examples
#' \dontrun{
#' pi.enlargement(
#'   ll = "C:/.../solution_file.tif", puid, ref.pa, st.pa = NULL,
#'   countries, bioregion, ref.pa.name = "all.pa", strict = FALSE,
#'   pae_threshold, project.name = "Example Project", version = "Version 1"
#' )
#' }
pi.enlargement <- function(ll, puid, ref.pa, st.pa = NULL,
                           countries, bioregion, ref.pa.name = "all.pa", strict = FALSE,
                           pae_threshold, project.name, version) {
  if (isTRUE(strict) && is.null(st.pa)) {
    stop("'st.pa' must be supplied when 'strict = TRUE'", call. = FALSE)
  }

  puid <- puid %>% mask(countries)

  # calculating Queen's distance layer (i.e., how far away is each PU from the nearest PA?)
  message(paste0("Calculating Queen's distances based on layer: ", ref.pa.name))

  all.pa.bi <- ifel(ref.pa > 0, 1, 0)
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

  # align CRS with that of planning units
  terra::crs(manh_sr) <- terra::crs(puid)
  assertthat::assert_that(terra::crs(puid) == terra::crs(manh_sr))
  mh.all <- manh_sr
  mh.strict.all <- manh_sr

  # Check CRS of distance rasters
  if (terra::crs(puid) != terra::crs(mh.all)) {
    mh.all <- mh.all |> terra::project(puid)
  }
  assertthat::assert_that(crs(puid) == crs(mh.all))

  if (terra::crs(puid) != terra::crs(mh.strict.all)) {
    mh.strict.all <- mh.strict.all |> terra::project(puid)
  }
  assertthat::assert_that(crs(puid) == crs(mh.strict.all))

  # union countries to represent continental extent
  co.dissolve <- sf::st_union(countries)

  # output containers
  out.europe <- data.frame()
  out.bioregion <- data.frame()
  out.countries <- data.frame()
  pb <- progress::progress_bar$new(total = length(ll))

  tic()
  for (i in 1:length(ll)) {
    pb$tick()

    # here this is the full solution (expansion priorities are back-calculated below)
    pa.sol <- paste0(ll[[i]]) %>%
      terra::rast() %>%
      mask(countries)

    # area selected in solution
    puid <- puid %>% mask(pa.sol)
    puid_area <- puid

    ref.pa <- ref.pa %>% mask(pa.sol)

    if (!is.null(st.pa)) {
      st.pa <- st.pa %>% mask(pa.sol)
    }

    if (strict == FALSE) {
      pa <- ref.pa
      mh.all <- mh.all %>% mask(pa.sol)
      mh <- mh.all
    } else if (strict == TRUE) {
      pa <- st.pa
      mh.strict.all <- mh.strict.all %>% mask(pa.sol)
      mh <- mh.strict.all
    }

    # pa.sol is the end solution; pa is the reference/baseline PA layer
    baseline <- pa
    out <- baseline
    final <- pa.sol

    # keep gains as selected area; preserve fractional cell values
    sol <- terra::ifel(final > baseline, final - baseline, 0)

    # keep losses separately
    loss <- terra::ifel(baseline > final, baseline - final, 0)

    # final solution already includes the baseline
    ov.sol.full <- final
    ov.sol.area <- exactextractr::exact_extract(ov.sol.full, co.dissolve, fun = "sum")

    # all selected gains, including both expansion and establishment.
    exp.sol.europe <- exactextractr::exact_extract(sol, co.dissolve, fun = "sum")

    # connected patches among gain cells
    patch_id <- terra::patches(
      terra::ifel(sol > 0, 1, NA),
      directions = 8,
      zeroAsNA = TRUE
    )

    # a patch qualifies if any of its cells are within the pae_threshold of the
    # reference PA raster defined above
    near_pa <- !is.na(mh) & mh <= pae_threshold
    patch_values <- terra::values(patch_id, mat = FALSE)
    near_values <- terra::values(near_pa, mat = FALSE)

    expansion_patch_ids <- unique(
      patch_values[near_values & !is.na(patch_values)]
    )

    # apply the patch classification while retaining fractional values in `sol`
    if (length(expansion_patch_ids) == 0) {
      pae.expansion <- terra::ifel(!is.na(patch_id), 0, 0)
    } else {
      pae.expansion <- terra::app(
        patch_id,
        fun = function(x) as.integer(!is.na(x) & x %in% expansion_patch_ids)
      )
    }
    pae.estab <- terra::ifel(pae.expansion == 1, 0, 1)

    # retain the original proportional gain values
    sol.expansion <- sol * pae.expansion
    sol.estab <- sol * pae.estab

    # normalized by proportion of total solution area, including area of existing PAs
    sol.expansion.europe <- exactextractr::exact_extract(sol.expansion, co.dissolve, fun = "sum")[[1]] / ov.sol.area[[1]]
    sol.estab.europe <- exactextractr::exact_extract(sol.estab, co.dissolve, fun = "sum")[[1]] / ov.sol.area[[1]]
    sol.exist.europe <- exactextractr::exact_extract(out, co.dissolve, fun = "sum")[[1]] / ov.sol.area[[1]]
    loss.europe <- exactextractr::exact_extract(loss, co.dissolve, fun = "sum")[[1]] / ov.sol.area[[1]]

    # normalized by proportion of expansion area
    sol.expansion.europe.exp <- exactextractr::exact_extract(sol.expansion, co.dissolve, fun = "sum")[[1]] / exp.sol.europe[[1]]
    sol.estab.europe.exp <- exactextractr::exact_extract(sol.estab, co.dissolve, fun = "sum")[[1]] / exp.sol.europe[[1]]

    # check that expansion, establishment, and existing add up to 100% of solution area
    total.fraction <- sum(c(
      sol.expansion.europe,
      sol.estab.europe,
      sol.exist.europe
    ), na.rm = TRUE) - loss.europe
    assertthat::assert_that(abs(total.fraction - 1) < 1e-2)

    # and that expansion and establishment, when normalized by total expansion area, add up to 100% of the expansion area
    expansion.fraction <- sum(c(
      sol.expansion.europe.exp,
      sol.estab.europe.exp
    ), na.rm = TRUE)
    assertthat::assert_that(abs(expansion.fraction - 1) < 1e-2)

    # format the results
    eur.expansion <- data.frame(
      project = paste0(project.name, "__", version),
      scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
      region = "Europe",
      variable = "Protected area expansion index",
      variable_variant = "Selected area in region",
      unit = "%",
      value = sol.expansion.europe
    )

    eur.exp2 <- eur.expansion %>% mutate(
      variable_variant = "Expansion area in region",
      unit = "%",
      value = sol.expansion.europe.exp
    )

    eur.establishment <- data.frame(
      project = paste0(project.name, "__", version),
      scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
      region = "Europe",
      variable = "Protected area establishment index",
      variable_variant = "Selected area in region",
      unit = "%",
      value = sol.estab.europe
    )

    eur.est2 <- eur.establishment %>% mutate(
      variable_variant = "Expansion area in region",
      unit = "%",
      value = sol.estab.europe.exp
    )

    out.europe <- rbind(out.europe, eur.expansion, eur.establishment, eur.exp2, eur.est2)

    # BIOREGION SUMMARIES

    for (j in 1:nrow(bioregion)) {
      region1 <- bioregion[j, ] %>% sf::st_make_valid()

      # calculate solution area for normalization
      sol.area.full <- ov.sol.full

      sol.area <- exactextractr::exact_extract(sol.area.full, region1, fun = "sum")
      exp.area <- exactextractr::exact_extract(sol, region1, fun = "sum")

      if (sol.area > 0) {
        # expansions
        sol.expansion.region <- exactextractr::exact_extract(sol.expansion, region1, fun = "sum") / sol.area
        sol.expansion.region.exp <- exactextractr::exact_extract(sol.expansion, region1, fun = "sum") / exp.area

        # newly established PAs
        sol.estab.region <- exactextractr::exact_extract(sol.estab, region1, fun = "sum") / sol.area
        sol.estab.region.exp <- exactextractr::exact_extract(sol.estab, region1, fun = "sum") / exp.area

        # percent of current PAs
        out.region <- exactextractr::exact_extract(out, region1, fun = "sum") / sol.area
        loss.region <- exactextractr::exact_extract(loss, region1, fun = "sum") / sol.area

        # check that the sum of PA expansion, PA establishments, and existing PAs is effectively 1
        # here we allow some room in case of rounding errors
        total.fraction <- sum(c(
          sol.expansion.region,
          sol.estab.region,
          out.region
        ), na.rm = TRUE) - loss.region
        assertthat::assert_that(abs(total.fraction - 1) < 1e-2)

        # and that expansion and establishment, when normalized by total expansion area, add up to 100% of the expansion area
        if (exp.area > 0) {
          expansion.fraction <- sum(c(
            sol.expansion.region.exp,
            sol.estab.region.exp
          ), na.rm = TRUE)
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

      bio.expansion <- data.frame(
        project = paste0(project.name, "__", version),
        scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
        region = paste0(as.character(bioregion[[j, "code"]])),
        variable = "Protected area expansion index",
        variable_variant = "Selected area in region",
        unit = "%",
        value = sol.expansion.region
      )

      bio.exp2 <- bio.expansion %>% mutate(
        variable_variant = "Expansion area in region",
        unit = "%",
        value = sol.expansion.region.exp
      )

      bio.establishment <- data.frame(
        project = paste0(project.name, "__", version),
        scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
        region = paste0(as.character(bioregion[[j, "code"]])),
        variable = "Protected area establishment index",
        variable_variant = "Selected area in region",
        unit = "%",
        value = sol.estab.region
      )

      bio.est2 <- bio.establishment %>% mutate(
        variable_variant = "Expansion area in region",
        unit = "%",
        value = sol.estab.region.exp
      )

      out.bioregion <- rbind(out.bioregion, bio.expansion, bio.establishment, bio.exp2, bio.est2)
    }

    # COUNTRY SUMMARIES

    for (j in 1:nrow(countries)) {
      region1 <- countries[j, ] %>% sf::st_make_valid()

      # calculate solution area for normalization
      sol.area.full <- ov.sol.full

      sol.area <- exactextractr::exact_extract(sol.area.full, region1, fun = "sum")
      exp.area <- exactextractr::exact_extract(sol, region1, fun = "sum")

      if (sol.area > 0) {
        # expansions
        sol.expansion.region <- exactextractr::exact_extract(sol.expansion, region1, fun = "sum") / sol.area
        sol.expansion.region.exp <- exactextractr::exact_extract(sol.expansion, region1, fun = "sum") / exp.area

        # newly established PAs
        sol.estab.region <- exactextractr::exact_extract(sol.estab, region1, fun = "sum") / sol.area
        sol.estab.region.exp <- exactextractr::exact_extract(sol.estab, region1, fun = "sum") / exp.area

        # percent of current PAs
        out.region <- exactextractr::exact_extract(out, region1, fun = "sum") / sol.area
        loss.region <- exactextractr::exact_extract(loss, region1, fun = "sum") / sol.area

        # check that the sum of PA expansion, PA establishments, and existing PAs is effectively 1
        # here we allow some room in case of rounding errors
        total.fraction <- sum(c(
          sol.expansion.region,
          sol.estab.region,
          out.region
        ), na.rm = TRUE) - loss.region
        assertthat::assert_that(abs(total.fraction - 1) < 1e-2)

        # and that expansion and establishment, when normalized by total expansion area, add up to 100% of the expansion area
        if (exp.area > 0) {
          expansion.fraction <- sum(c(
            sol.expansion.region.exp,
            sol.estab.region.exp
          ), na.rm = TRUE)
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

      c.expansion <- data.frame(
        project = paste0(project.name, "__", version),
        scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
        region = paste0(as.character(countries[[j, "code"]])),
        variable = "Protected area expansion index",
        variable_variant = "Selected area in region",
        unit = "%",
        value = sol.expansion.region
      )

      c.exp2 <- c.expansion %>% mutate(
        variable_variant = "Expansion area in region",
        unit = "%",
        value = sol.expansion.region.exp
      )

      c.establishment <- data.frame(
        project = paste0(project.name, "__", version),
        scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
        region = paste0(as.character(countries[[j, "code"]])),
        variable = "Protected area establishment index",
        variable_variant = "Selected area in region",
        unit = "%",
        value = sol.estab.region
      )

      c.est2 <- c.establishment %>% mutate(
        variable_variant = "Expansion area in region",
        unit = "%",
        value = sol.estab.region.exp
      )

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

  ofname <- paste0(outdir, "/", project.name, "__", version, "__indicator_paenlargement.csv")
  write.csv(out.enlargement.final, ofname, row.names = FALSE)

  message("PA enlargement and establishment calculations done. Don't forget to create/update variable definitions")
}
