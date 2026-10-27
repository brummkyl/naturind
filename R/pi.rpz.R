#' Mean Representation
#'
#' @param ll A character vector containing file paths of raster files (e.g., spatial prioritization outputs) in '.tif' format.
#' @param puid A SpatRaster layer with planning units
#' @param puid.num A SpatRaster layer with unique identifiers (i.e., cell positions) for each planning unit
#' @param countries An sf vector layer containing individual countries
#' @param feat.info A dataframe containing attributes for each feature. Useful for sub-setting e.g., species, habitats. Requires columns: 'name' = names of features; 'feat_type' = e.g., species, habitats, or any other character name specified in 'selection'; 'total_range_size' = double numeric quantifying the total distribution or range area of each feature
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
#' pi.rpz(ll = "C:/.../solution_file.tif", puid, puid.num, countries,
#'         feat.info, rij, project.name = "Example Project",
#'         version = "Version 1", countries.sub, countries.sub2,
#'         bioregion.sub, bioregion.sub2, selection = c("species", "habitats"), separate = FALSE,
#'         regional.metrics = TRUE, intermediate = TRUE)
#'}

pi.rpz <- function(ll, puid, puid.num, countries, feat.info = p, rij = y, project.name, version,
                                   countries.sub, countries.sub2, bioregion.sub, bioregion.sub2,
                                   selection = c("species", "habitats"), separate = FALSE,
                                   regional.metrics = TRUE, intermediate = TRUE){

  outdir <- paste0("performance_outputs/", version)
  dir.create(outdir, showWarnings = FALSE)

  # Deep-copy every terra object you might mutate, before doing anything else
  puid      <- terra::deepcopy(puid)
  puid.num  <- terra::deepcopy(puid.num)

  country.id <- countries.sub2$id[!is.na(countries.sub2$MS)]
  bioregion.id <- bioregion.sub2$id[!is.na(bioregion.sub2$BGR)]
  assertthat::assert_that(length(country.id) == length(bioregion.id))

  if(sum(!is.na(terra::values(puid))) != length(country.id)){
    message("Planning units appear to be missing regional associations in the '.sub2' data frames. Planning unit and solution rasters (ll) will be adjusted accordingly.")
    drop_cells <- setdiff(seq_len(terra::ncell(puid)), country.id)
    puid[drop_cells] <- NA
    puid.num[drop_cells] <- NA
  }

  p <- feat.info
  y <- rij

  ee <- ll %>% rast() %>% extend(puid) %>% mask(countries)
  existing_pa_names <- c(basename(ll))

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

  puid <- puid %>% mask(ee[[1]])
  puid_area <- puid

  puid.num <- puid.num %>% mask(puid)

  keep <- which(!is.na(values(puid)))
  col_id <- which(as.numeric(colnames(y)) %in% keep)
  #col_id <- which(colnames(y) %in% as.character(keep))
  y <- y[, col_id, drop=FALSE]

  # TODO: add trycatch, if this fails, reverse and mask based on y columns
  tryCatch({
    assertthat::assert_that(ncol(y) == sum(!is.na(values(puid_area))))
  }, error = function(e) {
    message("First assertion failed, reconciling puid/raster layers with the rij subset: ", conditionMessage(e))

    keep <- as.character(colnames(y))
    keep.cells <- which(!is.na(values(puid.num)) & as.character(values(puid.num)) %in% keep)

    condition <- setdiff(1:ncell(puid.num), keep.cells)

    puid.num[condition] <<- NA
    puid[condition] <<- NA
    puid_area[condition] <<- NA
    assertthat::assert_that(sum(!is.na(values(puid.num))) == sum(!is.na(values(puid_area))))
    assertthat::assert_that(ncol(y) == sum(!is.na(values(puid_area))))

    for(i in 1:nlyr(ee)){
      ee[[i]][condition] <<- NA
      assertthat::assert_that(sum(!is.na(values(ee[[i]]))) == sum(!is.na(values(puid_area))))
      assertthat::assert_that(sum(!is.na(values(ee[[i]]))) == sum(!is.na(values(ee[[1]]))))
      gc()
    }
  })

  #######################################################################################

  # EXISTING PROTECTED AREAS

  #######################################################################################

  selection <- unique(selection)
  n_passes <- if (separate) length(selection) else 1L

  # One group means process the selected types together.
  # Multiple groups means process each selected type in its own pass.
  groups <- if (separate) as.list(selection) else list(selection)

  # Load the .RData container for individual problems/solutions
  # Loop through target environments here
  evp.out.existing <- data.frame()
  co.ov.evp.out.existing <- data.frame()
  bi.ov.evp.out.existing <- data.frame()
  co.evp.out.existing <- data.frame()
  bi.evp.out.existing <- data.frame()
  ta.list <- list()
  pb <- progress::progress_bar$new(total = nlyr(ee) * 2)

  tic()
  for(h in groups){

    h_label <- paste(h, collapse = "_")

    for(i in 1:nlyr(ee)) {
      pb$tick()

      pa.sol <- ee[[i]]
      #pa.sol <- pa.sol %>% crop(puid) %>% mask(puid)
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
      features <- data.frame(name=x$name,ov.avail=x$total_range_size)

      # Remove features that are not EDGE
      #features <- features[features$name %in% edge.out$edge,]

      # Ensure features$name contains only values that exist in rownames(y)
      edge_features <- intersect(features$name, rownames(y))

      # Get row indices
      #row_id <- which(rownames(y) %in% features$name)
      row_id <- which(rownames(y) %in% edge_features)

      # Subset safely
      rij2 <- y[row_id, , drop = FALSE]

      # Dimensions check
      assertthat::assert_that(all(features$name %in% rownames(rij2)),
                              nrow(features) == nrow(rij2))

      # Calculate the area and total proportion of range conserved for each EDGE feature #

      # Account for proportional solutions by multiplying feature/PU by proportion of PU selected
      # Result is interpreted as the amount of each feature represented per PU
      D <- Diagonal(x = mat_prop)
      rij_sel <- rij2 %*% D
      rownames(rij_sel) <- rownames(rij2)
      colnames(rij_sel) <- colnames(rij2)

      # Calculate proportion conserved
      # sum of the amount selected for each feature
      ta <- data.frame(feature = names(Matrix::rowSums(rij_sel)),
                       pi = as.numeric(Matrix::rowSums(rij_sel)),
                       ov.avail = as.numeric(Matrix::rowSums(rij2)))

      drop.features <- ta$feature[ta$ov.avail == 0]
      row.id.drop <- which(!rownames(rij_sel) %in% drop.features)
      rij_sel <- rij_sel[row.id.drop, , drop = FALSE]
      rij2 <- rij2[row.id.drop, , drop = FALSE]
      features<-features[!features$name %in% drop.features,]

      ta <- ta[ta$ov.avail > 0,]
      ta$ratio <- ta$pi/ta$ov.avail
      ta$final <- pmin(ta$ratio, 1)
      evp <- mean(ta$final, na.rm=TRUE)

      ta.list[[paste(file.name, h, sep = "__")]] <- ta

      m <- data.frame(project = paste0(project.name,"__",version),
                      scenario = file.name,
                      region = "Europe",
                      variable = "Mean representation",
                      variable_variant = "Overall",
                      unit = "%",
                      value = evp)

      evp.out.existing <- rbind(evp.out.existing, m)

      if(regional.metrics == TRUE){

        # Countries
        for (j in 1:nrow(countries.sub)) {
          co.pu <- countries.sub2 %>% filter(names == countries.sub$names[j])

          if(nrow(co.pu) == 0){
            co.evp <- NA
          } else {
            # subset to retain selected PUs within a given country, then retain all available
            co.rij.sel = rij_sel[, colnames(rij_sel) %in% co.pu$id, drop=FALSE]
            co.rij.all = rij2[, colnames(rij2) %in% co.pu$id, drop=FALSE]

            # amount selected per feature in country
            co.sel = rowSums(co.rij.sel)
            #co.sel = as.data.frame(co.sel)

            # total available per feature in country
            co.available = rowSums(co.rij.all)

            # maximum potential contribution of country toward feature target
            co.ta = as.data.frame(cbind(co.sel,co.available))

            assertthat::assert_that(all(nrow(co.ta) == nrow(ta)))
            assertthat::assert_that(all(rownames(co.ta) %in% ta$feature))
            co.ov.ta = as.data.frame(cbind(co.ta, ov.avail = ta$ov.avail))

            co.ta = co.ta[co.ta$co.available > 0,]

            if(nrow(co.ta) == 0) {
              co.evp <- NA
            } else {
              assertthat::assert_that(all(co.ta$co.available > 0))
            }

            # country-level MTA
            co.ta.ratio = co.ta$co.sel / co.ta$co.available
            co.ta.final = pmin(co.ta.ratio, 1)
            co.evp = mean(co.ta.final, na.rm = TRUE)

            n <- data.frame(project = paste0(project.name,"__",version),
                            scenario = file.name,
                            region = paste0(as.character(countries.sub[[j,"names"]])),
                            variable = "Mean representation",
                            variable_variant = "Region burden share",
                            unit = "%",
                            value = co.evp)

            co.evp.out.existing <- rbind(co.evp.out.existing, n)

            # contribution of each country to meeting the overall target of each feature
            co.ov.out.existing <- data.frame(co.ov.ta$co.sel / co.ov.ta$ov.avail)
            co.ov.ta.ratio = data.frame(rownames(co.ov.ta), co.ov.out.existing)
            colnames(co.ov.ta.ratio) <- c("features", paste0(countries.sub$names[j]))

            if (ncol(co.ov.evp.out.existing) == 0) {
              # initialize with first result
              co.ov.evp.out.existing <- co.ov.ta.ratio
            } else {
              # check row compatibility before cbind
              if (nrow(co.ov.evp.out.existing) == nrow(co.ov.ta.ratio)) {
                co.ov.evp.out.existing <- cbind(co.ov.evp.out.existing, co.ov.ta.ratio[2])
              } else {
                stop("Row mismatch: cannot cbind different lengths")
              }
            }

          }
        }

        xx = data.frame(sum = rowSums(co.ov.evp.out.existing[,2:ncol(co.ov.evp.out.existing)]))
        xx2 <- data.frame(sum = pmin(as.numeric(xx[[1]]), 1, na.rm = FALSE))
        assertthat::assert_that(abs(evp - mean(xx2$sum)) < 1e-6)


        outdir.rpz <- paste0(outdir, "/RPZ/")
        dir.create(outdir.rpz, showWarnings = FALSE)

        #ofname <- paste0(outdir, "/", "EVP/", "NaturaConnect","__", version, "__",paste0("Existing ", file.name, " (", filtered[z],")"), "__indicator_EVP_Country_Contribution_Percentages.csv")
        ofname <- paste0(outdir.rpz, project.name,"__", version, "__",paste0(file.name, "__indicator_RPZ_",h_label,"Country_Contribution_Percentages.csv"))

        write.csv(co.ov.evp.out.existing, ofname, row.names = FALSE)
        co.ov.evp.out.existing <- data.frame()

        gc()

        # Bioregions
        for (k in 1:nrow(bioregion.sub)) {
          bi.pu <- bioregion.sub2 %>% filter(names == bioregion.sub$names[k])

          if(nrow(bi.pu) == 0){
            bi.evp <- NA
          } else {
            # subset to retain selected PUs within a given country, then retain all available
            bi.rij.sel = rij_sel[, colnames(rij_sel) %in% bi.pu$id, drop=FALSE]
            bi.rij.all = rij2[, colnames(rij2) %in% bi.pu$id, drop=FALSE]

            # amount selected per feature in country
            bi.sel = rowSums(bi.rij.sel)

            # total available per feature in country
            bi.available = rowSums(bi.rij.all)

            # maximum potential contribution of country toward feature target
            bi.ta = as.data.frame(cbind(bi.sel,bi.available))

            assertthat::assert_that(all(nrow(bi.ta) == nrow(ta)))
            assertthat::assert_that(all(rownames(bi.ta) %in% ta$feature))
            bi.ov.ta = as.data.frame(cbind(bi.ta, ov.avail = ta$ov.avail))

            bi.ta = bi.ta[bi.ta$bi.available > 0,]

            if(nrow(bi.ta) == 0) {
              bi.evp <- NA
            } else {
              assertthat::assert_that(all(bi.ta$bi.available > 0))
            }

            # fraction of country contribution achieved
            bi.ta.ratio = bi.ta$bi.sel / bi.ta$bi.available
            bi.ta.final = pmin(bi.ta.ratio, 1)
            bi.evp = mean(bi.ta.final, na.rm = TRUE)

            o <- data.frame(project = paste0(project.name,"__",version),
                            scenario = file.name,
                            region = paste0(as.character(bioregion.sub[[k,"names"]])),
                            variable = "Mean representation",
                            variable_variant = "Region burden share",
                            unit = "%",
                            value = bi.evp)

            bi.evp.out.existing <- rbind(bi.evp.out.existing, o)

            # contribution of each country to meeting the overall target of each feature
            bi.ov.out.existing <- data.frame(bi.ov.ta$bi.sel / bi.ov.ta$ov.avail)
            bi.ov.ta.ratio = data.frame(rownames(bi.ov.ta), bi.ov.out.existing)
            colnames(bi.ov.ta.ratio) <- c("features", paste0(bioregion.sub$names[k]))

            if (ncol(bi.ov.evp.out.existing) == 0) {
              # initialize with first result
              bi.ov.evp.out.existing <- bi.ov.ta.ratio
            } else {
              # check row compatibility before cbind
              if (nrow(bi.ov.evp.out.existing) == nrow(bi.ov.ta.ratio)) {
                bi.ov.evp.out.existing <- cbind(bi.ov.evp.out.existing, bi.ov.ta.ratio[2])
              } else {
                stop("Row mismatch: cannot cbind different lengths")
              }
            }


          }
        }

        yy = data.frame(sum = rowSums(bi.ov.evp.out.existing[,2:ncol(bi.ov.evp.out.existing)]))
        yy2 <- data.frame(sum = pmin(as.numeric(yy[[1]]), 1, na.rm = FALSE))
        assertthat::assert_that(abs(evp - mean(yy2$sum)) < 1e-6)

        #ofname <- paste0(outdir, "/", "EVP/", "NaturaConnect","__", version, "__",paste0("Existing ", file.name, " (", filtered[z],")"), "__indicator_EVP_Bioregion_Contribution_Percentages.csv")
        ofname <- paste0(outdir.rpz, project.name,"__", version, "__",paste0(file.name, "__indicator_RPZ_",h_label,"Bioregion_Contribution_Percentages.csv"))

        write.csv(bi.ov.evp.out.existing, ofname, row.names = FALSE)
        bi.ov.evp.out.existing <- data.frame()

        gc()
      }
    }
  }
  toc()

  # Bind results for bioregions and countries, then export to output directory
  out.evp.final <- rbind(evp.out.existing, co.evp.out.existing, bi.evp.out.existing)
  out.evp.final$value <- out.evp.final$value * 100

  #ofname <- paste0(outdir, "/", "NaturaConnect","__",version, "__indicator_evolutionarypotential.csv")
    ofname <- paste0(outdir, "/", project.name,"__",version, "__", h_label, "indicator_meanrepresentation.csv")
    write.csv(out.evp.final, ofname, row.names = FALSE)

  if(intermediate == TRUE){
    return(ta.list)
  }

  message("Mean representation calculation done. Don't forget to create/update variable definitions")
}
