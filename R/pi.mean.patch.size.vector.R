# ll.v <- list(hist, curr)
# names(ll.v) <- c("2014", "2024")

#' Mean/Median Patch Size (Vector)
#'
#' @param ll.v A named list of SpatVector layers representing protected areas
#' @param countries An sf vector layer containing individual countries, with names specified in a column named "code"
#' @param bioregion An sf vector layer containing individual biogeographic regions, with names specified in a column named "code"
#' @param project.name A character vector specifying the name of your project
#' @param version A character vector specifying the version (e.g., date) of your workflow
#'
#' @returns A formatted .csv file saved to the "performance_outputs/version" directory
#' @export
#'
#' @examples
#' \dontrun{
#' pi.mean.patch.size.vector(ll.v, countries, bioregion,
#'   project.name = "Example Project", version = "Version 1"
#' )
#' }
pi.mean.patch.size.vector <- function(ll.v, countries, bioregion, project.name, version) {
  # union countries to represent continental extent
  co.dissolve <- sf::st_union(countries)

  # output containers
  size.europe.existing <- data.frame()
  size.bioregion.existing <- data.frame()
  size.countries.existing <- data.frame()

  tic()
  for (i in names(ll.v)) {
    message("Processing file: ", i)

    pa.sol <- ll.v[[i]]
    file.name <- names(ll.v[i])

    pa.sol.crop <- pa.sol %>%
      sf::st_make_valid() %>%
      sf::st_union() %>%
      sf::st_cast("POLYGON") %>%
      sf::st_sf(geometry = .) %>%
      sf::st_intersection(co.dissolve)

    pa.sol.crop$area_km2 <- sf::st_area(pa.sol.crop) / 1e6

    mean.patch.size <- as.numeric(mean(pa.sol.crop$area_km2, na.rm = T))
    median.patch.size <- as.numeric(median(pa.sol.crop$area_km2, na.rm = T))

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

      region1 <- bioregion[j, ] %>% sf::st_make_valid()

      bi.pa.sol.crop <- pa.sol.crop %>%
        sf::st_intersection(region1)

      if (nrow(bi.pa.sol.crop) > 0) {
        bi.pa.sol.crop$area_km2 <- sf::st_area(bi.pa.sol.crop) / 1e6
        mean.patch.size <- as.numeric(mean(bi.pa.sol.crop$area_km2, na.rm = T))
        median.patch.size <- as.numeric(median(bi.pa.sol.crop$area_km2, na.rm = T))
      } else if (nrow(bi.pa.sol.crop) == 0) {
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

      region1 <- countries[j, ] %>% sf::st_make_valid()

      co.pa.sol.crop <- pa.sol.crop %>%
        sf::st_intersection(region1)

      if (nrow(co.pa.sol.crop) > 0) {
        co.pa.sol.crop$area_km2 <- sf::st_area(co.pa.sol.crop) / 1e6
        mean.patch.size <- as.numeric(mean(co.pa.sol.crop$area_km2, na.rm = T))
        median.patch.size <- as.numeric(median(co.pa.sol.crop$area_km2, na.rm = T))
      } else if (nrow(co.pa.sol.crop) == 0) {
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
