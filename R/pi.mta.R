#' Mean Target Achievement
#'
#' @param ll A character vector containing file paths of raster files (e.g., spatial prioritization outputs) in '.tif' format.
#' @param puid A SpatRaster layer with planning units
#' @param puid.num TBD
#' @param countries TBD
#' @param all.pa A SpatRaster layer with all existing protected areas
#' @param st.pa A SpatRaster layer with all existing strict protected areas
#' @param n2k A SpatRaster layer with all existing Natura 2000 protected areas (or any other unique subset of PAs)
#' @param feat.info A dataframe containing attributes for each feature. Useful for sub-setting e.g., species, habitats.
#' @param rij TBD
#' @param project.name A character vector specifying the name of your project
#' @param version A character vector specifying the version (e.g., date) of your workflow
#' @param countries.sub TBD
#' @param countries.sub2 TBD
#' @param bioregion.sub TBD
#' @param bioregion.sub2 TBD
#' @param outdir Output directory
#'
#' @returns A formatted .csv file saved to the "performance_outputs/version" directory
#' @export
#'
#' @examples
#' \dontrun{
#' pi.mta(ll = "C:/.../solution_file.tif", puid, puid.num, countries, all.pa,
#'         st.pa, n2k, feat.info, rij, project.name = "Example Project",
#'         version = "Version 1", countries.sub, countries.sub2,
#'         bioregion.sub, bioregion.sub2, outdir)
#'}

pi.mta <- function(ll, puid, puid.num, countries, all.pa, st.pa, n2k, feat.info = p, rij = y, project.name, version,
                   countries.sub, countries.sub2, bioregion.sub, bioregion.sub2, outdir){

  outdir <- paste0("performance_outputs/", version)
  dir.create(outdir, showWarnings = FALSE)

  # Deep-copy every terra object you might mutate, before doing anything else
  puid      <- terra::deepcopy(puid)
  puid.num  <- terra::deepcopy(puid.num)
  all.pa    <- terra::deepcopy(all.pa)
  st.pa     <- terra::deepcopy(st.pa)
  n2k       <- terra::deepcopy(n2k)

  puid <- puid %>% mask(countries)
  puid.num <- puid.num %>% mask(countries)

  p <- feat.info
  y <- rij

  ee <- ll %>% terra::rast() %>% terra::extend(puid) #%>% mask(countries)
  existing_pa_names <- c(basename(ll))

  ee <- c(ee, all.pa, st.pa, n2k)
  existing_pa_names <- c(existing_pa_names, "Existing all PA", "Existing strict PA", "Existing N2K")

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

  #ee <- ee %>% mask(ee[[1]])

  puid <- puid %>% terra::mask(ee[[1]])
  puid_area <- puid

  puid.num <- puid.num %>% terra::mask(puid)

  keep <- which(!is.na(terra::values(puid)))
  col_id <- which(as.numeric(colnames(y)) %in% keep)
  y <- y[, col_id, drop=FALSE]

  # TODO: add trycatch, if this fails, reverse and mask based on y columns
  tryCatch({
    assertthat::assert_that(ncol(y) == sum(!is.na(terra::values(puid_area))))
  }, error = function(e) {
    message("First assertion failed, reconciling puid/raster layers with the rij subset: ", conditionMessage(e))

    keep <- as.character(colnames(y))
    keep.cells <- which(!is.na(values(puid.num)) & as.character(terra::values(puid.num)) %in% keep)

    condition <- setdiff(1:terra::ncell(puid.num), keep.cells)

    puid.num[condition] <<- NA
    puid[condition] <<- NA
    puid_area[condition] <<- NA
    assertthat::assert_that(sum(!is.na(terra::values(puid.num))) == sum(!is.na(terra::values(puid_area))))
    assertthat::assert_that(ncol(y) == sum(!is.na(terra::values(puid_area))))

    for(i in 1:terra::nlyr(ee)){
      ee[[i]][condition] <<- NA
      assertthat::assert_that(sum(!is.na(terra::values(ee[[i]]))) == sum(!is.na(terra::values(puid_area))))
      assertthat::assert_that(sum(!is.na(terra::values(ee[[i]]))) == sum(!is.na(terra::values(ee[[1]]))))
      gc()
    }
  })

  #######################################################################################

  # EXISTING PROTECTED AREAS

  #######################################################################################

  # Initialize output containers
  mta.out.existing <- data.frame()
  co.ov.mta.out.existing <- data.frame()
  bi.ov.mta.out.existing <- data.frame()
  co.mta.out.existing <- data.frame()
  bi.mta.out.existing <- data.frame()
  pb <- progress::progress_bar$new(total = terra::nlyr(ee) * 2) # Here we multiply by the total length of h specified in the first loop below

  tictoc::tic()
  for(h in c("species", "habitats")){

    for(i in 1:terra::nlyr(ee)) {
      pb$tick()
      pa.sol <- ee[[i]]
      #pa.sol <- pa.sol %>% mask(puid)
      file.name <- existing_pa_names[[i]]

      mat <- terra::values(pa.sol, mat = FALSE,na.rm=TRUE)
      mat_pu <- as.numeric(mat)
      assertthat::assert_that(length(mat_pu) == sum(!is.na(puid[])))

      land_area <- terra::values(puid, mat = FALSE, na.rm = TRUE)
      # what proportion of land area is currently protected?
      mat_prop <- pmin(mat_pu / as.numeric(land_area), 1)
      mat_prop[is.na(mat_prop)] <- 0

      # filter species from habitats
      x <- p[p$feat_type %in% h,]

      # Extract the subset features and targets
      features <- data.frame(name=x$name,prop=x$target_area/x$total_range_size)
      targets <- data.frame(name = x$name,value = x$target_area)

      # Remove features for which targets were set to 0
      # This ensures we aren't artificially inflating the mean target achievement
      targets <- targets[targets$value > 0,]
      features <- features[features$name %in% targets$name,]

      # Get row indices to subset rij for feature/target alignment
      row_id <- which(rownames(y) %in% targets$name)
      # Subset as dense matrix
      rij2 <- y[row_id, , drop = FALSE]

      # Dimensions check
      assertthat::assert_that(all(features$name %in% targets$name),
                              nrow(features) == nrow(rij2))

      # Calculate the overall mean target achievement per feature #

      # Account for proportional solutions by multiplying feature/PU by proportion of PU selected
      # Result is interpreted as the amount of each feature represented per PU
      D <- Matrix::Diagonal(x = mat_prop)
      rij_sel <- rij2 %*% D
      rownames(rij_sel) <- rownames(rij2)
      colnames(rij_sel) <- colnames(rij2)

      # Calculate target achievements
      # sum of the amount selected for each feature
      ta <- data.frame(feature = names(Matrix::rowSums(rij_sel)),
                       pi = as.numeric(Matrix::rowSums(rij_sel)))

      corrected.targets <- data.frame(name = features$name, value = rowSums(rij2) * features$prop)

      # Join in ti
      ta <- ta |> dplyr::left_join(corrected.targets |> dplyr::select(name, value) |>
                                     dplyr::rename(ti = value),
                                   by = c("feature" = "name") )

      # mean( min( Conserved amount  for species / Target for species, 1) )
      # pi = total amount of each feature selected in the solution
      # ti = initial (absolute) target
      ta <- ta[ta$ti > 0,]
      ta.ratio <- ta$pi/ta$ti
      ta.final <- pmin(ta.ratio, 1)
      mta <- mean(ta.final, na.rm=TRUE)

      m <- data.frame(project = paste0(project.name,"__",version),
                      scenario = file.name,
                      region = "Europe",
                      variable = paste0("Mean target achievement (", h, ")"),
                      variable_variant = "Overall",
                      unit = "%",
                      value = mta)

      mta.out.existing <- rbind(mta.out.existing, m)

      gc()

      # Countries
      for (j in 1:nrow(countries.sub)) {
        co.pu <- countries.sub2 %>% dplyr::filter(names == countries.sub$names[j])

        if(nrow(co.pu) == 0){
          co.mta <- NA
        } else {
          # subset to retain selected PUs within a given region, and retain available
          # area of each feature per region
          co.rij.sel = rij_sel[, colnames(rij_sel) %in% co.pu$id, drop=FALSE]
          co.rij.all = rij2[, colnames(rij2) %in% co.pu$id, drop=FALSE]

          # amount selected per feature per region
          co.sel = rowSums(co.rij.sel)
          # total available per feature per region
          co.available = rowSums(co.rij.all)

          # regional target adjustment
          co.target = data.frame(co.target = co.available * features$prop)  # same prop as for overall MTA calculation
          co.ta = as.data.frame(cbind(co.sel,co.available,co.target))
          co.ta = co.ta[rownames(co.ta) %in% ta$feature,]
          assertthat::assert_that(all(nrow(co.ta) == nrow(ta)))
          assertthat::assert_that(all(rownames(co.ta) %in% ta$feature))

          # regional contribution to overall target achievement
          co.ov.ta = as.data.frame(cbind(co.ta, ov.target = ta$ti))

          # remove features not found in the given region
          co.ta = co.ta[co.ta$co.target > 0,]
          assertthat::assert_that(all(co.ta$co.available > 0))

          # country-level MTA
          co.ta.ratio = co.ta$co.sel / co.ta$co.target
          co.ta.final = pmin(co.ta.ratio, 1)
          co.mta = mean(co.ta.final, na.rm = TRUE)

          n <- data.frame(project = paste0(project.name,"__",version),
                          scenario =  file.name,
                          region = paste0(as.character(countries.sub[[j,"names"]])),
                          variable = paste0("Mean target achievement (", h, ")"),
                          variable_variant = "Region burden share",
                          unit = "%",
                          value = co.mta)

          co.mta.out.existing <- rbind(co.mta.out.existing, n)

          # contribution of each country to meeting the overall target of each feature
          co.ov.out <- data.frame(co.ov.ta$co.sel / co.ov.ta$ov.target)
          co.ov.ta.ratio = data.frame(rownames(co.ov.ta), co.ov.out)
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
      xx = data.frame(sum = rowSums(co.ov.mta.out.existing[,2:ncol(co.ov.mta.out.existing)]))
      xx2 <- data.frame(sum = pmin(as.numeric(xx[[1]]), 1, na.rm = FALSE))
      assertthat::assert_that(abs(mta - mean(xx2$sum)) < 1e-8)

      outdir.mta <- paste0(outdir, "/MTA/")
      dir.create(outdir.mta, showWarnings = FALSE)

      #if(consensus.sol == 1){
      #  ofname <- paste0(outdir.consensus, "/", "MTA/", project.name,"__", version, "__",paste0(file.name, "__indicator_MTA_",h,"_Country_Contribution_Percentages.csv"))
      #} else if (consensus.sol == 0){
        ofname <- paste0(outdir.mta, project.name,"__", version, "__",paste0(file.name, "__indicator_MTA_",h,"_Country_Contribution_Percentages.csv"))
      #}

      utils::write.csv(co.ov.mta.out.existing, ofname, row.names = FALSE)
      co.ov.mta.out.existing <- data.frame()

      gc()

      # Bioregions
      for (k in 1:nrow(bioregion.sub)) {
        bi.pu <- bioregion.sub2 %>% dplyr::filter(names == bioregion.sub$names[k])

        if(nrow(bi.pu) == 0){
          bi.mta <- NA
        } else {
          # subset to retain selected PUs within a given region
          bi.rij.sel = rij_sel[, colnames(rij_sel) %in% bi.pu$id, drop=FALSE]
          bi.rij.all = rij2[, colnames(rij2) %in% bi.pu$id, drop=FALSE]

          # amount selected per feature per region
          bi.sel = rowSums(bi.rij.sel)
          # total available per feature per region
          bi.available = rowSums(bi.rij.all)

          # maximum potential contribution of region toward feature target
          bi.target = data.frame(bi.target = bi.available * features$prop)  # same prop as for overall MTA calculations
          bi.ta = as.data.frame(cbind(bi.sel,bi.available,bi.target))
          bi.ta = bi.ta[rownames(bi.ta) %in% ta$feature,]
          assertthat::assert_that(all(nrow(bi.ta) == nrow(ta)))
          assertthat::assert_that(all(rownames(bi.ta) %in% ta$feature))

          bi.ov.ta = as.data.frame(cbind(bi.ta, ov.target = ta$ti))

          bi.ta = bi.ta[bi.ta$bi.target > 0,]
          assertthat::assert_that(all(bi.ta$bi.available > 0))

          # fraction of region contribution achieved
          bi.ta.ratio = bi.ta$bi.sel / bi.ta$bi.target
          bi.ta.final = pmin(bi.ta.ratio, 1)
          bi.mta = mean(bi.ta.final, na.rm = TRUE)

          o <- data.frame(project = paste0(project.name,"__",version),
                          scenario = file.name,# " (", filtered[z],")"),
                          region = paste0(as.character(bioregion.sub[[k,"names"]])),
                          variable = paste0("Mean target achievement (", h, ")"),
                          variable_variant = "Region burden share",
                          unit = "%",
                          value = bi.mta)

          bi.mta.out.existing <- rbind(bi.mta.out.existing, o)

          # contribution of each country to meeting the overall target of each feature
          bi.ov.out <- data.frame(bi.ov.ta$bi.sel / bi.ov.ta$ov.target)
          bi.ov.ta.ratio = data.frame(rownames(bi.ov.ta), bi.ov.out)
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
      yy = data.frame(sum = rowSums(bi.ov.mta.out.existing[,2:ncol(bi.ov.mta.out.existing)]))
      yy2 <- data.frame(sum = pmin(as.numeric(yy[[1]]), 1, na.rm = FALSE))
      assertthat::assert_that(abs(mta - mean(yy2$sum)) < 1e-8)

      #if(consensus.sol == 1){
      #  ofname <- paste0(outdir.consensus, "/", "MTA/", project.name,"__", version, "__", paste0(file.name, "__indicator_MTA_",h,"_Bioregion_Contribution_Percentages.csv"))
      #} else if (consensus.sol == 0){
        ofname <- paste0(outdir.mta, project.name,"__", version, "__", paste0(file.name, "__indicator_MTA_",h,"_Bioregion_Contribution_Percentages.csv"))
      #}

      utils::write.csv(bi.ov.mta.out.existing, ofname, row.names = FALSE)
      bi.ov.mta.out.existing <- data.frame()

      gc()

    }
  }
  tictoc::toc()

  # Bind results for overall, bioregions, and countries. export to output directory
  out.mta.final <- rbind(mta.out.existing, co.mta.out.existing, bi.mta.out.existing)
  out.mta.final$value <- out.mta.final$value * 100

  #if(consensus.sol == 1){
  #  ofname <- paste0(outdir.consensus, "/", project.name,"__",version, "__indicator_meantargetachievement.csv")
  #} else if (consensus.sol == 0){
    ofname <- paste0(outdir, "/", project.name,"__",version, "__indicator_meantargetachievement.csv")
  #}

  utils::write.csv(out.mta.final, ofname, row.names = FALSE)

  message("Mean target achievement calculation done")

}
