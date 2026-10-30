#' Mean/Median Patch Size (Vector)
#'
#' @param ll.v A named list of sf vector layers representing contiguous, non-overlapping protected area polygons. Requires a column named 'area_km2', reflecting the area of each unique patch. For help with preprocessing, see the 'create.flat.layer' function.
#' @param countries An sf vector layer containing individual countries, with names specified in a column named "code"
#' @param bioregion An sf vector layer containing individual biogeographic regions, with names specified in a column named "code"
#' @param min.area An integer, specifying the minimum patch area (in square km) to retain for the calculations. This parameter is useful for removing erroneous patch fragments that may have been generated during geoprocessing. Defaults to 5e-6, which retains patches with an area larger than 5 square meters.
#' @param project.name A character vector specifying the name of your project
#' @param version A character vector specifying the version (e.g., date) of your workflow
#'
#' @returns A formatted .csv file saved to the "performance_outputs/version" directory
#' @export
#'
#' @examples
#' \dontrun{
#' pi.mean.patch.size.vector(ll.v, countries, bioregion, min.area = 5e-6,
#'   project.name = "Example Project", version = "Version 1"
#' )
#' }
pi.mean.patch.size.vector <- function(ll.v, countries, bioregion, min.area = 5e-6, project.name, version) {
  # union countries to represent continental extent
  countries <- countries %>% sf::st_make_valid()
  bioregion <- bioregion %>% sf::st_make_valid()

  co.dissolve <- sf::st_union(countries)

  # output containers
  size.europe.existing <- data.frame()
  size.bioregion.existing <- data.frame()
  size.countries.existing <- data.frame()

  tic()
  for (i in names(ll.v)) {
    message("Processing file: ", i)

    pa.sol.crop <- ll.v[[i]]
    file.name <- names(ll.v[i])

    # geometry-only copy for the intersections
    pa.geom <- sf::st_sf(geometry = sf::st_geometry(pa.sol.crop))

    areas <- pa.sol.crop$area_km2[pa.sol.crop$area_km2 > min.area]

    mean.patch.size <- as.numeric(mean(areas, na.rm = T))
    median.patch.size <- as.numeric(median(areas, na.rm = T))

    o <- data.frame(
      project = paste0(project.name, "__", version),
      scenario = file.name,
      region = "Europe",
      variable = "Mean patch size",
      variable_variant = "Mean patch size in region",
      unit = "km2",
      value = mean.patch.size
    )

    o2 <- o %>% mutate(
      variable = "Median patch size",
      variable_variant = "Median patch size in region",
      value = median.patch.size
    )

    size.europe.existing <- rbind(size.europe.existing, o, o2)

    # fresh progress bar per region
    pb_bioregion <- progress::progress_bar$new(
      total = nrow(bioregion),
      format = "    [:bar] :percent | :bioregion",
      clear = FALSE
    )

    for (j in 1:nrow(bioregion)) {
      # loop through the bioregions
      pb_bioregion$tick(tokens = list(bioregion = bioregion[[j, "code"]]))

      region1 <- bioregion[j, ]

      pre.filt <- sf::st_filter(pa.geom, region1)

      if (nrow(pre.filt) > 0) {
        bi.pa.sol.crop <- pre.filt %>%
          sf::st_intersection(region1) %>%
          sf::st_collection_extract("POLYGON") %>% # drop lines/points from touching borders
          sf::st_cast("POLYGON") # one row per polygon

        areas <- as.numeric(sf::st_area(bi.pa.sol.crop)) / 1e6
        areas <- areas[areas > min.area]
      } else {
        areas <- numeric(0)
      }

      if (length(areas) > 0) {
        mean.patch.size <- mean(areas)
        median.patch.size <- median(areas)
      } else {
        mean.patch.size <- NA
        median.patch.size <- NA
      }

      o <- data.frame(
        project = paste0(project.name, "__", version),
        scenario = file.name,
        region = paste0(as.character(bioregion[[j, "code"]])),
        variable = "Mean patch size",
        variable_variant = "Mean patch size in region",
        unit = "km2",
        value = mean.patch.size
      )

      o2 <- o %>% mutate(
        variable = "Median patch size",
        variable_variant = "Median patch size in region",
        value = median.patch.size
      )

      size.bioregion.existing <- rbind(size.bioregion.existing, o, o2)
      gc()
    }

    # fresh progress bar per region
    pb_countries <- progress::progress_bar$new(
      total = nrow(countries),
      format = "    [:bar] :percent | :country",
      clear = FALSE
    )

    for (j in 1:nrow(countries)) {
      # loop through the countries
      pb_countries$tick(tokens = list(country = countries[[j, "code"]]))

      region1 <- countries[j, ]

      pre.filt <- sf::st_filter(pa.geom, region1)

      if (nrow(pre.filt) > 0) {
        co.pa.sol.crop <- pre.filt %>%
          sf::st_intersection(region1) %>%
          sf::st_collection_extract("POLYGON") %>% # drop lines/points from touching borders
          sf::st_cast("POLYGON") # one row per polygon

        areas <- as.numeric(sf::st_area(co.pa.sol.crop)) / 1e6
        areas <- areas[areas > min.area]
      } else {
        areas <- numeric(0)
      }

      if (length(areas) > 0) {
        mean.patch.size <- mean(areas)
        median.patch.size <- median(areas)
      } else {
        mean.patch.size <- NA
        median.patch.size <- NA
      }

      o <- data.frame(
        project = paste0(project.name, "__", version),
        scenario = file.name,
        region = paste0(as.character(countries[[j, "code"]])),
        variable = "Mean patch size",
        variable_variant = "Mean patch size in region",
        unit = "km2",
        value = mean.patch.size
      )

      o2 <- o %>% mutate(
        variable = "Median patch size",
        variable_variant = "Median patch size in region",
        value = median.patch.size
      )

      size.countries.existing <- rbind(size.countries.existing, o, o2)
      gc()
    }
    gc()
  }
  toc()

  # Bind results for bioregions and countries, then export to output directory
  out.size.final <- rbind(size.europe.existing, size.bioregion.existing, size.countries.existing)

  outdir <- paste0("performance_outputs/", version)
  dir.create(outdir, showWarnings = FALSE)

  ofname <- paste0(outdir, "/", project.name, "__", version, "__indicator_relativepatchsize_vector.csv")
  write.csv(out.size.final, ofname, row.names = FALSE)

  message("Patch size (vector) calculations done. Don't forget to create/update variable definitions")
}
