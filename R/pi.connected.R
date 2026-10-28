#' Absolute (ProtCon) and Relative (RelConn) Protected and Connected Land Area
#'
#' @param ll A character vector containing file paths of raster files (e.g., spatial prioritization outputs) in '.tif' format.
#' @param puid A SpatRaster layer, reflecting the amount of land area per planning unit
#' @param countries An sf vector layer containing individual countries, with names specified in a column named "code"
#' @param bioregion An sf vector layer containing individual biogeographic regions, with names specified in a column named "code"
#' @param conn.layer A connectivity layer, which will internally be converted to a resistance raster.
#' @param target.res An integer defining the resolution (in meters) at which the connectivity metrics should be calculated. This value should be greater than (or equal to) the resolution of your ll, puid, and conn.layer rasters. If this value is greater than the resolution of your input rasters, it should be completely divisible. For example, if the resolution of your input layers is 1000m, then target.res may be set to "2000", "3000", "4000", "5000", etc...
#' @param bin.threshold A double numeric, ranging between 0 and 1, defining the threshold for binarizing a proportional protected area raster
#' @param euc.threshold An integer specifying the dispersal distance (in meters) used to calculate connectivity metrics.
#' @param n.cores Number of cores to use in parallel processing. The actual number of cores used will be taken as the minimum of 'n.cores' and length(ll)
#' @param project.name A character vector specifying the name of your project
#' @param version A character vector specifying the version (e.g., date) of your workflow
#'
#' @returns A formatted .csv file saved to the "performance_outputs/version" directory
#' @export
#'
#' @examples
#' \dontrun{
#' pi.connected(ll, puid, countries, bioregion, conn.layer,
#'   target.res = 5000,
#'   bin.threshold = 0.5, euc.threshold = 10000, n.cores = 2,
#'   project.name = "Example Project", version = "Version 1"
#' )
#' }
pi.connected <- function(ll, puid, countries, bioregion, conn.layer, target.res = 5000, bin.threshold = 0.5, euc.threshold = 10000, n.cores = 2, project.name, version) {
  puid <- puid %>% mask(countries)
  puid_area <- puid
  target_res <- target.res
  prop_binary_thresh <- bin.threshold
  euclid_threshold <- euc.threshold

  cat(paste0(
    "Starting run at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " with ", length(ll), " raster solutions.\n"
  ))
  flush.console()

  # align projected CRS
  if (terra::crs(puid) != terra::crs(conn.layer)) {
    conn.layer <- conn.layer |> terra::project(puid)
  }
  assertthat::assert_that(crs(puid) == crs(conn.layer))

  # rescale from 0 to 1, and then resample to
  conn.layer <- terra::clamp((conn.layer - global(conn.layer, "min", na.rm = TRUE)[1, 1]) /
    (global(conn.layer, "max", na.rm = TRUE)[1, 1] - global(conn.layer, "min", na.rm = TRUE)[1, 1]), 0, 1)
  conn.layer <- terra::resample(conn.layer, puid, method = "bilinear")
  assertthat::assert_that(all.equal(res(puid), res(conn.layer)))

  # invert to convert from connectivity to resistance
  resist.layer.temp <- 1 - conn.layer
  resist.layer.temp <- resist.layer.temp %>% mask(puid)

  fact <- target_res / res(resist.layer.temp)[1]

  if (fact > 1) {
    message(paste0("Aggregating resistance layer by factor ", fact, " (to ~", target_res / 1000, " km resolution)"))
    resist.layer.temp <- terra::aggregate(resist.layer.temp, fact = fact, fun = "mean", na.rm = TRUE)
  }

  # calculate the mean resistance of the raster
  mean_resistance <- terra::global(resist.layer.temp, fun = "mean", na.rm = TRUE)[1, 1]
  assertthat::assert_that(mean_resistance > 0, msg = "Average resistance should be larger than 0!")

  # convert from the euclidean distance threshold (10km) to a resistance-weighted cost-distance threshold
  # euclidean distance is defined in the parameter settings
  cost_threshold <- euclid_threshold * mean_resistance
  resist.layer.temp2 <- ifel(is.na(resist.layer.temp), 1, resist.layer.temp)
  assertthat::assert_that(!any(is.na(values(resist.layer.temp2))), msg = "Resistance raster contains NA values!")

  # convert to rasterlayer format to accomodate the MK_protconn function
  resist.layer <- raster(resist.layer.temp2)
  assertthat::assert_that(class(resist.layer)[1] == "RasterLayer", msg = "Resistance layer must be of class 'RasterLayer'")

  # loop through solution files
  # convert each raster solution into a binary polygon vector using a selection frequency threshold as defined in parameter settings

  ncores <- min(n.cores, length(ll))
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  cat("Running on", ncores, "cores.\n")

  # clean up cluster and files on function exit
  on.exit(
    {
      stopCluster(cl)
      terra::tmpFiles(remove = TRUE)
    },
    add = TRUE
  )

  clusterEvalQ(cl, {
    library(terra)
    library(sf)
    library(dplyr)
    library(Makurhini)
    library(tictoc)
    library(foreach)
    library(raster)
  })

  # set terra options on each worker (tempdir on local scratch + limit memfrac)
  clusterEvalQ(cl, {
    # create a per-worker temp dir (uses TMPDIR or /scratch as available)
    td <- Sys.getenv("TMPDIR")
    if (td == "") td <- "/scratch"
    mytmp <- file.path(td, paste0("terra_tmp_worker_", Sys.getpid()))
    dir.create(mytmp, showWarnings = FALSE, recursive = TRUE)
    terraOptions(tmpdir = mytmp, memfrac = 0.4, todisk = TRUE) # memfrac ~ fraction of available memory for terra
    # reduce verbose output
    options("terra.mem" = TRUE)
  })

  # save in-memory resistance raster to a temp file
  resist_path <- tempfile(fileext = ".tif")
  terra::writeRaster(resist.layer.temp2, resist_path, overwrite = TRUE)
  puid_path <- tempfile(fileext = ".tif")
  terra::writeRaster(puid_area, puid_path, overwrite = TRUE)

  # export variables to cluster workers
  clusterExport(cl, c(
    "ll", "bioregion", "countries",
    "cost_threshold", "resist_path",
    "prop_binary_thresh", "project.name", "version",
    "puid_path", "target_res"
  ),
  envir = environment()
  )

  # output containers
  out.europe <- data.frame()
  out.bioregion <- data.frame()
  out.country <- data.frame()

  tic()

  log_message <- function(msg) {
    cat(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg, "\n"))
    flush.console()
  }

  # starting message
  log_message(paste0("Starting run with ", length(ll), " raster solutions."))

  # set batch size (same as ncores)
  batch_size <- ncores
  n_solutions <- length(ll)
  all_results <- list()

  # process batches
  for (batch_start in seq(1, n_solutions, by = batch_size)) {
    batch_end <- min(batch_start + batch_size - 1, n_solutions)
    batch_indices <- batch_start:batch_end

    log_message(paste0(
      "Starting batch ", (batch_start - 1) %/% batch_size + 1,
      " (solutions ", batch_start, "-", batch_end, ")"
    ))

    results_batch <- foreach(
      i = batch_indices,
      .packages = c("terra", "sf", "raster", "dplyr", "tictoc", "assertthat", "Makurhini", "tools")
    ) %dopar% {
      puid_area <- terra::rast(puid_path)

      pa.sol <- paste0(ll[[i]]) %>%
        terra::rast() %>%
        mask(puid_area)

      # Coarsen from 1km to some lower resolution (e.g. 2km or 5km)
      fact <- target_res / res(pa.sol)[1]

      # Aggregate to coarser grid
      # Use 'fun = "mean"' to preserve proportion of selected cells, or 'fun = "max"' for binary presence
      pa.sol <- terra::aggregate(pa.sol, fact = fact, fun = "mean", na.rm = TRUE)

      pa.sol <- ifel(pa.sol >= prop_binary_thresh, 1, NA)
      assertthat::assert_that(global(pa.sol, fun = "sum", na.rm = TRUE) > 0)

      pa_patches <- patches(pa.sol, directions = 8, allowGaps = FALSE, zeroAsNA = TRUE)
      pa_patches <- as.polygons(pa_patches) %>% st_as_sf()

      resist.layer.p <- terra::rast(resist_path)

      out.europe.local <- data.frame()
      out.bioregion.local <- data.frame()
      out.country.local <- data.frame()

      # Europe

      europe <- st_union(countries) %>%
        st_make_valid() %>%
        st_as_sf()
      pa.sol.europe <- pa.sol %>%
        terra::crop(europe) %>%
        terra::mask(europe)

      if (is.na(global(pa.sol.europe, fun = "sum", na.rm = T)[[1]])) { # PA selected in solutions should be at least 1km2 for ProtConn function to recognize them as 'nodes'
        rel.conn <- "No protected areas found in region"
        prot.conn <- "No protected areas found in region"
      } else {
        conn.results <- MK_ProtConn(
          nodes = pa_patches,
          region = europe,
          area_unit = "km2",
          distance = list(type = "least-cost", resistance = resist.layer.p),
          distance_thresholds = cost_threshold,
          probability = 0.5,
          plot = FALSE, delta = FALSE, intern = FALSE, protconn_bound = FALSE
        )

        # Try to extract the "Viewer Panel" - some small countries may output different formats if they only contain a single PA
        if ("Protected Connected (Viewer Panel)" %in% names(conn.results)) {
          conn.df <- conn.results[["Protected Connected (Viewer Panel)"]]
        } else {
          conn.df <- data.frame(conn.results)
        }

        names(conn.df) <- c("Index", "Value", "ProtConn indicator", "Percentage")
        val_protconn <- filter(conn.df, `ProtConn indicator` == "ProtConn")$Percentage
        val_prot <- filter(conn.df, `ProtConn indicator` == "Prot")$Percentage
        val_relconn <- filter(conn.df, `ProtConn indicator` == "RelConn")$Percentage

        if (!is.na(val_relconn)) {
          rel.conn <- val_relconn
          val_prot
        } else {
          rel.conn <- "Contains a single PA patch (OR) No PA patches are connected"
        }

        prot.conn <- val_protconn
      }

      n <- data.frame(
        project = paste0(project.name, "__", version),
        scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
        region = "Europe",
        variable = "Protected area connectivity",
        variable_variant = "RelConn",
        unit = "%",
        value = paste0(rel.conn)
      )

      n2 <- data.frame(
        project = paste0(project.name, "__", version),
        scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
        region = "Europe",
        variable = "Protected area connectivity",
        variable_variant = "ProtConn",
        unit = "%",
        value = paste0(prot.conn)
      )

      out.europe.local <- rbind(out.europe.local, n, n2)
      rm(conn.results, conn.df, pa.sol.europe)
      gc()

      # Bioregion
      for (j in 1:nrow(bioregion)) {
        region1 <- bioregion[j, ] %>% st_make_valid()
        assertthat::assert_that(st_is_valid(region1))

        pa.sol.region <- pa.sol %>%
          terra::crop(region1) %>%
          terra::mask(region1)

        if (is.na(global(pa.sol.region, fun = "sum", na.rm = T)[[1]])) { # PA selected in solutions should be at least 1km2 for ProtConn function to recognize them as 'nodes'
          rel.conn <- "No protected areas found in region"
          prot.conn <- "No protected areas found in region"
        } else {
          conn.results <- MK_ProtConn(
            nodes = pa_patches,
            region = region1,
            area_unit = "km2",
            distance = list(type = "least-cost", resistance = resist.layer.p),
            distance_thresholds = cost_threshold,
            probability = 0.5,
            transboundary = 10000,
            transboundary_type = "region",
            plot = FALSE, delta = FALSE, intern = FALSE, protconn_bound = FALSE
          )

          # Try to extract the "Viewer Panel" - some small countries may output different formats if they only contain a single PA
          if ("Protected Connected (Viewer Panel)" %in% names(conn.results)) {
            conn.df <- conn.results[["Protected Connected (Viewer Panel)"]]
          } else {
            conn.df <- data.frame(conn.results)
          }

          names(conn.df) <- c("Index", "Value", "ProtConn indicator", "Percentage")
          val_protconn <- filter(conn.df, `ProtConn indicator` == "ProtConn")$Percentage
          val_prot <- filter(conn.df, `ProtConn indicator` == "Prot")$Percentage
          val_relconn <- filter(conn.df, `ProtConn indicator` == "RelConn")$Percentage

          if (!is.na(val_relconn)) {
            rel.conn <- val_relconn
            val_prot
          } else {
            rel.conn <- "Contains a single PA patch (OR) No PA patches are connected"
          }

          prot.conn <- val_protconn
        }

        o <- data.frame(
          project = paste0(project.name, "__", version),
          scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
          region = paste0(as.character(bioregion[[j, "code"]])),
          variable = "Protected area connectivity",
          variable_variant = "RelConn",
          unit = "%",
          value = paste0(rel.conn)
        )

        o2 <- data.frame(
          project = paste0(project.name, "__", version),
          scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
          region = paste0(as.character(bioregion[[j, "code"]])),
          variable = "Protected area connectivity",
          variable_variant = "ProtConn",
          unit = "%",
          value = paste0(prot.conn)
        )

        out.bioregion.local <- rbind(out.bioregion.local, o, o2)
        rm(conn.results, conn.df, pa.sol.region)
        gc()
      }
      # Country

      for (k in 1:nrow(countries)) {
        region2 <- countries[k, ] %>% st_make_valid()
        assertthat::assert_that(st_is_valid(region2))

        pa.sol.region <- pa.sol %>%
          terra::crop(region2) %>%
          terra::mask(region2)

        if (is.na(global(pa.sol.region, fun = "sum", na.rm = T)[[1]])) { # PA selected in solutions should be at least 1km2 for ProtConn function to recognize them as 'nodes'
          rel.conn <- "No protected areas found in region"
          prot.conn <- "No protected areas found in region"
        } else {
          conn.results <- MK_ProtConn(
            nodes = pa_patches,
            region = region2,
            area_unit = "km2",
            distance = list(type = "least-cost", resistance = resist.layer.p),
            distance_thresholds = cost_threshold,
            probability = 0.5,
            transboundary = 10000,
            transboundary_type = "region",
            plot = FALSE, delta = FALSE, intern = FALSE, protconn_bound = FALSE
          )

          # Try to extract the "Viewer Panel" - some small countries may output different formats if they only contain a single PA
          if ("Protected Connected (Viewer Panel)" %in% names(conn.results)) {
            conn.df <- conn.results[["Protected Connected (Viewer Panel)"]]
          } else {
            conn.df <- data.frame(conn.results)
          }

          names(conn.df) <- c("Index", "Value", "ProtConn indicator", "Percentage")
          val_protconn <- filter(conn.df, `ProtConn indicator` == "ProtConn")$Percentage
          val_prot <- filter(conn.df, `ProtConn indicator` == "Prot")$Percentage
          val_relconn <- filter(conn.df, `ProtConn indicator` == "RelConn")$Percentage

          if (!is.na(val_relconn)) {
            rel.conn <- val_relconn
          } else {
            rel.conn <- "Contains a single PA patch (OR) No PA patches are connected"
          }

          prot.conn <- val_protconn
        }

        p <- data.frame(
          project = paste0(project.name, "__", version),
          scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
          region = paste0(as.character(countries[[k, "code"]])),
          variable = "Protected area connectivity",
          variable_variant = "RelConn",
          unit = "%",
          value = paste0(rel.conn)
        )

        p2 <- data.frame(
          project = paste0(project.name, "__", version),
          scenario = paste0(tools::file_path_sans_ext(basename(ll[[i]]))),
          region = paste0(as.character(countries[[k, "code"]])),
          variable = "Protected area connectivity",
          variable_variant = "ProtConn",
          unit = "%",
          value = paste0(prot.conn)
        )

        out.country.local <- rbind(out.country.local, p, p2)
        rm(conn.results, conn.df, pa.sol.region)
        gc()
      }
      rm(pa.sol, pa_patches, resist.layer.p)
      gc()

      list(europe = out.europe.local, bioregion = out.bioregion.local, country = out.country.local)
    }

    all_results <- c(all_results, results_batch)

    log_message(paste0("Finished batch ", (batch_start - 1) %/% batch_size + 1))
  }

  out.europe <- do.call(rbind, lapply(all_results, `[[`, "europe"))
  out.bioregion <- do.call(rbind, lapply(all_results, `[[`, "bioregion"))
  out.country <- do.call(rbind, lapply(all_results, `[[`, "country"))

  # Bind results for europe, bioregions and countries, then export
  out.conn.final <- rbind(out.europe, out.bioregion, out.country)

  outdir <- paste0("performance_outputs/", version)
  dir.create(outdir, showWarnings = FALSE)

  ofname <- paste0(outdir, "/", project.name, "__", version, "__indicator_protectedconnected_", target_res / 1000, "km.csv")
  write.csv(out.conn.final, ofname, row.names = FALSE)

  toc()

  message("Protected area connectivity calculations done. Don't forget to create/update variable definitions")
}
