#' Mean Target Achievement
#'
#' @param ll A character vector containing file paths of raster files (e.g., spatial prioritization outputs) in '.tif' format.
#' @param puid A SpatRaster layer, reflecting the amount of land area per planning unit
#' @param puid.num A SpatRaster layer with unique identifiers (i.e., cell positions) for each planning unit
#' @param countries An sf vector layer containing individual countries
#' @param feat.info A dataframe containing attributes for each feature. Useful for sub-setting e.g., species, habitats. Requires columns: 'name' = names of features; 'feat_type' = e.g., species, habitats, or any other character name specified in 'selection'; 'total_range_size' = double numeric quantifying the total distribution or range area of each feature; 'target area' = double numeric quantifying the targeted area to be conserved for each feature
#' @param rij An rij matrix with features as rows and unique planning unit identifiers as columns
#' @param project.name A character vector specifying the name of your project
#' @param version A character vector specifying the version (e.g., date) of your workflow
#' @param countries.sub A dataframe containing two columns: 'names' = country names; 'codes' = code names for countries
#' @param countries.sub2 A dataframe containing two columns: 'id' as in puid.num's unique identifiers; 'names' as in countries.sub
#' @param bioregion.sub A dataframe containing two columns: 'names' = region names; 'codes' = code names for regions
#' @param bioregion.sub2 A dataframe containing two columns: 'id' as in puid.num's unique identifiers; 'names' as in bioregion.sub
#' @param selection A character vector specifying which feature subsets you would like to calculate MTA for (e.g., c("species", "habitats"); "species"; "habitats")
#' @param separate A logical indicating whether you would like to calculate MTA separately for each group specified in the 'selection' argument. e.g., If selection = c("species", "habitats") and separate = TRUE, MTA will be calculated first for species, and then separately for habitats. In this example, if separate = FALSE, then a single estimate of MTA will be calculated for the combination of species and habitats
#' @param regional.metrics A logical of whether metrics should be calculated for countries and biogeographic regions
#' @param intermediate A logical of whether intermediate values should be returned for each individual feature
#'
#' @returns A formatted .csv file saved to the "performance_outputs/version" directory
#' @export
#'
#' @examples
#' \dontrun{
#' pi.mta(
#'   ll = "C:/.../solution_file.tif", puid, puid.num, countries,
#'   feat.info, rij, project.name = "Example Project",
#'   version = "Version 1", countries.sub, countries.sub2,
#'   bioregion.sub, bioregion.sub2, selection = c("species", "habitats"), separate = FALSE,
#'   regional.metrics = TRUE, intermediate = TRUE
#' )
#' }
pi.mta <- function(ll, puid, puid.num, countries, feat.info = p, rij = y, project.name, version,
                   countries.sub, countries.sub2, bioregion.sub, bioregion.sub2, selection = c("species", "habitats"), separate = FALSE,
                   regional.metrics = TRUE, intermediate = TRUE) {
  outdir <- paste0("performance_outputs/", version)
  dir.create(outdir, showWarnings = FALSE)

  # Deep-copy every terra object you might mutate, before doing anything else
  puid <- terra::deepcopy(puid)
  puid.num <- terra::deepcopy(puid.num)

  puid <- puid %>% mask(countries)
  puid.num <- puid.num %>% mask(countries)

  country.id <- countries.sub2$id[!is.na(countries.sub2$MS)]
  bioregion.id <- bioregion.sub2$id[!is.na(bioregion.sub2$BGR)]
  assertthat::assert_that(length(country.id) == length(bioregion.id))

  if (sum(!is.na(terra::values(puid))) != length(country.id)) {
    message("Planning units appear to be missing regional associations in the '.sub2' data frames. Planning unit and solution rasters (ll) will be adjusted accordingly.")
    drop_cells <- setdiff(seq_len(terra::ncell(puid)), country.id)
    puid[drop_cells] <- NA
    puid.num[drop_cells] <- NA
  }

  p <- feat.info
  y <- rij

  ee <- ll %>%
    terra::rast() %>%
    terra::extend(puid) %>%
    mask(countries)
  existing_pa_names <- c(basename(ll))

  # Check whether all layers already have the same NA footprint as ee[[1]]
  # This is used to adjust baseline existing PA calculations, for consistency with the provided planning solutions
  # One use case, for example, is to allow shifting between Europe and EU extents
  reference_na <- is.na(terra::values(ee[[1]]))

  na_match <- all(vapply(
    seq_len(terra::nlyr(ee)),
    function(i) {
      identical(reference_na, is.na(terra::values(ee[[i]])))
    },
    logical(1)
  ))

  # Only mask the stack if at least one layer differs
  if (!na_match) {
    message("NA footprints differ between ee layers; masking stack to align with the first input layer...")
    ee <- ee %>% terra::mask(ee[[1]])
  }

  ee <- ee %>% mask(puid)

  puid <- puid %>% terra::mask(ee[[1]])
  puid_area <- puid

  puid.num <- puid.num %>% terra::mask(puid)

  keep <- which(!is.na(terra::values(puid)))
  col_id <- which(as.numeric(colnames(y)) %in% keep)
  y <- y[, col_id, drop = FALSE]

  tryCatch(
    {
      assertthat::assert_that(ncol(y) == sum(!is.na(terra::values(puid_area))))
    },
    error = function(e) {
      message("First assertion failed, reconciling puid/raster layers with the rij subset: ", conditionMessage(e))

      keep <- as.character(colnames(y))
      keep.cells <- which(!is.na(values(puid.num)) & as.character(terra::values(puid.num)) %in% keep)

      condition <- setdiff(1:terra::ncell(puid.num), keep.cells)

      puid.num[condition] <<- NA
      puid[condition] <<- NA
      puid_area[condition] <<- NA
      assertthat::assert_that(sum(!is.na(terra::values(puid.num))) == sum(!is.na(terra::values(puid_area))))
      assertthat::assert_that(ncol(y) == sum(!is.na(terra::values(puid_area))))

      for (i in 1:terra::nlyr(ee)) {
        ee[[i]][condition] <<- NA
        assertthat::assert_that(sum(!is.na(terra::values(ee[[i]]))) == sum(!is.na(terra::values(puid_area))))
        assertthat::assert_that(sum(!is.na(terra::values(ee[[i]]))) == sum(!is.na(terra::values(ee[[1]]))))
        gc()
      }
    }
  )

  message(paste0("Number of planning units = ", ncol(y)))

  selection <- unique(selection)
  n_passes <- if (separate) length(selection) else 1L

  # One group means to process the selected types together.
  # Multiple groups means to process each selected type in its own pass.
  groups <- if (separate) as.list(selection) else list(selection)

  # output containers
  co.ov.mta.out.existing <- data.frame()
  bi.ov.mta.out.existing <- data.frame()
  ta.list <- list()
  pb <- progress::progress_bar$new(total = terra::nlyr(ee) * n_passes) # Here we multiply by the total length of h specified in the first loop below, depending on whether separate = TRUE or FALSE

  tictoc::tic()

  for (h in groups) {
    h_label <- paste(h, collapse = "_")

    mta.out.existing <- data.frame()
    co.mta.out.existing <- data.frame()
    bi.mta.out.existing <- data.frame()

    for (i in 1:terra::nlyr(ee)) {
      pb$tick()
      pa.sol <- ee[[i]]
      file.name <- existing_pa_names[[i]]

      mat <- terra::values(pa.sol, mat = FALSE, na.rm = TRUE)
      mat_pu <- as.numeric(mat)
      assertthat::assert_that(length(mat_pu) == sum(!is.na(puid[])))

      land_area <- terra::values(puid, mat = FALSE, na.rm = TRUE)

      # what proportion of land area is protected?
      mat_prop <- pmin(mat_pu / as.numeric(land_area), 1)
      mat_prop[is.na(mat_prop)] <- 0

      # filter species from habitats
      x <- p[p$feat_type %in% h, ]

      # Extract the subset features and targets
      features <- data.frame(name = x$name, prop = x$target_area / x$total_range_size)
      targets <- data.frame(name = x$name, value = x$target_area)

      # Remove features for which targets were set to 0
      # This ensures we aren't artificially inflating the mean target achievement
      targets <- targets[targets$value > 0, ]
      features <- features[features$name %in% targets$name, ]

      # Get row indices to subset rij for feature/target alignment
      row_id <- which(rownames(y) %in% targets$name)
      # Subset as dense matrix
      rij2 <- y[row_id, , drop = FALSE]

      # Dimensions check
      assertthat::assert_that(
        all(features$name %in% targets$name),
        nrow(features) == nrow(rij2)
      )

      # Calculate the overall mean target achievement per feature #

      # Account for proportional solutions by multiplying feature/PU by proportion of PU selected
      # Result is interpreted as the amount of each feature represented per PU
      D <- Matrix::Diagonal(x = mat_prop)
      rij_sel <- rij2 %*% D
      rownames(rij_sel) <- rownames(rij2)
      colnames(rij_sel) <- colnames(rij2)

      # Calculate target achievements
      # sum of the amount selected for each feature
      ta <- data.frame(
        feature = names(Matrix::rowSums(rij_sel)),
        pi = as.numeric(Matrix::rowSums(rij_sel))
      )

      # recalculate targets in relation to aligned puid, ll, and rij files
      corrected.targets <- data.frame(name = features$name, value = Matrix::rowSums(rij2) * features$prop)
      drop.features <- corrected.targets$name[corrected.targets$value == 0]
      row.id.drop <- which(!rownames(rij_sel) %in% drop.features)
      rij_sel <- rij_sel[row.id.drop, , drop = FALSE]
      rij2 <- rij2[row.id.drop, , drop = FALSE]
      features <- features[!features$name %in% drop.features, ]

      # Join in ti
      ta <- ta |> dplyr::left_join(
        corrected.targets |> dplyr::select(name, value) |>
          dplyr::rename(ti = value),
        by = c("feature" = "name")
      )

      # mean( min( Conserved amount  for species / Target for species, 1) )
      # pi = total amount of each feature selected in the solution
      # ti = initial (absolute) target
      ta <- ta[ta$ti > 0, ]
      ta$ratio <- ta$pi / ta$ti # could exceed 100%
      ta$final <- pmin(ta$ratio, 1) # capped at 100%
      mta <- mean(ta$final, na.rm = TRUE)

      ta.list[[paste(file.name, h_label, sep = "__")]] <- ta

      m <- data.frame(
        project = paste0(project.name, "__", version),
        scenario = file.name,
        region = "Europe",
        variable = paste0("Mean target achievement (", h_label, ")"),
        variable_variant = "Overall",
        unit = "%",
        value = mta
      )

      mta.out.existing <- rbind(mta.out.existing, m)

      gc()

      if (regional.metrics == TRUE) {
        # Countries
        for (j in 1:nrow(countries.sub)) {
          co.pu <- countries.sub2 %>% dplyr::filter(names == countries.sub$names[j])

          if (nrow(co.pu) == 0) {
            co.mta <- NA
          } else {
            # subset to retain selected PUs within a given region, and retain available
            # area of each feature per region
            co.rij.sel <- rij_sel[, colnames(rij_sel) %in% co.pu$id, drop = FALSE]
            co.rij.all <- rij2[, colnames(rij2) %in% co.pu$id, drop = FALSE]

            # amount selected per feature per region
            co.sel <- rowSums(co.rij.sel)
            # total available per feature per region
            co.available <- rowSums(co.rij.all)

            # regional target adjustment
            co.target <- data.frame(co.target = co.available * features$prop) # same prop as for overall MTA calculation
            co.ta <- as.data.frame(cbind(co.sel, co.available, co.target))
            co.ta <- co.ta[rownames(co.ta) %in% ta$feature, ]
            assertthat::assert_that(all(nrow(co.ta) == nrow(ta)))
            assertthat::assert_that(all(rownames(co.ta) %in% ta$feature))

            # regional contribution to overall target achievement
            co.ov.ta <- as.data.frame(cbind(co.ta, ov.target = ta$ti))

            # remove features not found in the given region
            co.ta <- co.ta[co.ta$co.target > 0, ]
            assertthat::assert_that(all(co.ta$co.available > 0))

            # country-level MTA
            co.ta.ratio <- co.ta$co.sel / co.ta$co.target
            co.ta.final <- pmin(co.ta.ratio, 1)
            co.mta <- mean(co.ta.final, na.rm = TRUE)

            n <- data.frame(
              project = paste0(project.name, "__", version),
              scenario = file.name,
              region = paste0(as.character(countries.sub[[j, "names"]])),
              variable = paste0("Mean target achievement (", h_label, ")"),
              variable_variant = "Region burden share",
              unit = "%",
              value = co.mta
            )

            co.mta.out.existing <- rbind(co.mta.out.existing, n)

            # contribution of each country to meeting the overall target of each feature
            co.ov.out <- data.frame(co.ov.ta$co.sel / co.ov.ta$ov.target)
            co.ov.ta.ratio <- data.frame(rownames(co.ov.ta), co.ov.out)
            colnames(co.ov.ta.ratio) <- c("features", paste0(countries.sub$names[j]))

            if (ncol(co.ov.mta.out.existing) == 0) {
              # initialize with first result
              co.ov.mta.out.existing <- co.ov.ta.ratio
            } else {
              # check row compatibility before cbind
              if (nrow(co.ov.mta.out.existing) == nrow(co.ov.ta.ratio)) {
                co.ov.mta.out.existing <- cbind(co.ov.mta.out.existing, co.ov.ta.ratio[2])
              } else {
                stop("Row mismatch: cannot cbind different lengths")
              }
            }
          }
        }

        # assert that overall MTA can accurately be back-calculated from these regional contributions
        xx <- data.frame(sum = rowSums(co.ov.mta.out.existing[, 2:ncol(co.ov.mta.out.existing)]))
        xx2 <- data.frame(sum = pmin(as.numeric(xx[[1]]), 1, na.rm = FALSE))
        assertthat::assert_that(abs(mta - mean(xx2$sum)) < 1e-6)

        outdir.mta <- paste0(outdir, "/MTA/")
        dir.create(outdir.mta, showWarnings = FALSE)

        ofname <- paste0(outdir.mta, project.name, "__", version, "__", paste0(file.name, "__indicator_MTA_", h_label, "_Country_Contribution_Percentages.csv"))
        utils::write.csv(co.ov.mta.out.existing, ofname, row.names = FALSE)
        co.ov.mta.out.existing <- data.frame()

        gc()

        # Bioregions
        for (k in 1:nrow(bioregion.sub)) {
          bi.pu <- bioregion.sub2 %>% dplyr::filter(names == bioregion.sub$names[k])

          if (nrow(bi.pu) == 0) {
            bi.mta <- NA
          } else {
            # subset to retain selected PUs within a given region
            bi.rij.sel <- rij_sel[, colnames(rij_sel) %in% bi.pu$id, drop = FALSE]
            bi.rij.all <- rij2[, colnames(rij2) %in% bi.pu$id, drop = FALSE]

            # amount selected per feature per region
            bi.sel <- rowSums(bi.rij.sel)
            # total available per feature per region
            bi.available <- rowSums(bi.rij.all)

            # maximum potential contribution of region toward feature target
            bi.target <- data.frame(bi.target = bi.available * features$prop) # same prop as for overall MTA calculations
            bi.ta <- as.data.frame(cbind(bi.sel, bi.available, bi.target))
            bi.ta <- bi.ta[rownames(bi.ta) %in% ta$feature, ]
            assertthat::assert_that(all(nrow(bi.ta) == nrow(ta)))
            assertthat::assert_that(all(rownames(bi.ta) %in% ta$feature))

            bi.ov.ta <- as.data.frame(cbind(bi.ta, ov.target = ta$ti))

            bi.ta <- bi.ta[bi.ta$bi.target > 0, ]
            assertthat::assert_that(all(bi.ta$bi.available > 0))

            # fraction of region contribution achieved
            bi.ta.ratio <- bi.ta$bi.sel / bi.ta$bi.target
            bi.ta.final <- pmin(bi.ta.ratio, 1)
            bi.mta <- mean(bi.ta.final, na.rm = TRUE)

            o <- data.frame(
              project = paste0(project.name, "__", version),
              scenario = file.name, # " (", filtered[z],")"),
              region = paste0(as.character(bioregion.sub[[k, "names"]])),
              variable = paste0("Mean target achievement (", h_label, ")"),
              variable_variant = "Region burden share",
              unit = "%",
              value = bi.mta
            )

            bi.mta.out.existing <- rbind(bi.mta.out.existing, o)

            # contribution of each country to meeting the overall target of each feature
            bi.ov.out <- data.frame(bi.ov.ta$bi.sel / bi.ov.ta$ov.target)
            bi.ov.ta.ratio <- data.frame(rownames(bi.ov.ta), bi.ov.out)
            colnames(bi.ov.ta.ratio) <- c("features", paste0(bioregion.sub$names[k]))

            if (ncol(bi.ov.mta.out.existing) == 0) {
              # initialize with first result
              bi.ov.mta.out.existing <- bi.ov.ta.ratio
            } else {
              # check row compatibility before cbind
              if (nrow(bi.ov.mta.out.existing) == nrow(bi.ov.ta.ratio)) {
                bi.ov.mta.out.existing <- cbind(bi.ov.mta.out.existing, bi.ov.ta.ratio[2])
              } else {
                stop("Row mismatch: cannot cbind different lengths")
              }
            }
          }
        }
        # assert that overall MTA can accurately be back-calculated from these regional contributions
        yy <- data.frame(sum = rowSums(bi.ov.mta.out.existing[, 2:ncol(bi.ov.mta.out.existing)]))
        yy2 <- data.frame(sum = pmin(as.numeric(yy[[1]]), 1, na.rm = FALSE))
        assertthat::assert_that(abs(mta - mean(yy2$sum)) < 1e-6)

        ofname <- paste0(outdir.mta, project.name, "__", version, "__", paste0(file.name, "__indicator_MTA_", h_label, "_Bioregion_Contribution_Percentages.csv"))
        utils::write.csv(bi.ov.mta.out.existing, ofname, row.names = FALSE)
        bi.ov.mta.out.existing <- data.frame()

        gc()
      }
    }

    # Bind results for overall, bioregions, and countries. export to output directory
    out.mta.final <- rbind(mta.out.existing, co.mta.out.existing, bi.mta.out.existing)
    out.mta.final$value <- out.mta.final$value * 100

    ofname <- paste0(outdir, "/", project.name, "__", version, "__indicator_", h_label, "_meantargetachievement.csv")
    utils::write.csv(out.mta.final, ofname, row.names = FALSE)
  }
  tictoc::toc()

  if (intermediate == TRUE) {
    return(ta.list)
  }

  message("Mean target achievement calculation done")
}
