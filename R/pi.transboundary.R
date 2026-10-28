#' Transboundary Protection
#'
#' @param ll A character vector containing file paths of raster files (e.g., spatial prioritization outputs) in '.tif' format.
#' @param puid A SpatRaster layer, reflecting the amount of land area per planning unit
#' @param countries An sf vector layer containing individual countries
#' @param tb.buffers An sf vector layer containing a transboundary region. See the 'create.tb.buffer' function.
#' @param project.name A character vector specifying the name of your project
#' @param version A character vector specifying the version (e.g., date) of your workflow
#'
#' @returns A formatted .csv file saved to the "performance_outputs/version" directory
#' @export
#'
#' @examples
#' \dontrun{
#' pi.transboundary(
#'   ll = "C:/.../solution_file.tif", puid, countries,
#'   tb.buffers, project.name = "Example Project",
#'   version = "Version 1"
#' )
#' }
pi.transboundary <- function(ll, puid, countries, tb.buffers, project.name, version) {
  ee <- ll %>%
    rast() %>%
    extend(puid) %>%
    mask(countries)
  existing_pa_names <- c(basename(ll))

  # Transform to ensure alignment of CRS, and run assertions
  if (terra::crs(puid) != st_crs(tb.buffers)) {
    tb.buffers <- tb.buffers |> sf::st_transform(crs = st_crs(puid))
  }
  assertthat::assert_that(st_crs(tb.buffers) == st_crs(puid))

  puid <- puid %>% mask(countries)

  # Make an area grid - ha
  puid <- puid %>% mask(ee[[1]])
  puid_area <- puid %>% mask(ee[[1]])

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

  co.dissolve <- sf::st_union(countries)

  message(paste0("Number of planning units in solution files = ", sum(!is.na(values(ee[[1]])))))
  message(paste0("Number of planning units in puid file = ", sum(!is.na(values(puid)))))

  tb.europe.existing <- data.frame()
  pb <- progress::progress_bar$new(total = nlyr(ee))

  tic()
  # loop through the solution rasters
  for (i in 1:nlyr(ee)) {
    pb$tick()
    pa.sol <- ee[[i]]
    file.name <- existing_pa_names[[i]]

    # Area of protected areas
    pa.sol.area <- pa.sol
    # Total area of protected areas contained in the solution
    ov.area <- exactextractr::exact_extract(pa.sol.area, co.dissolve, fun = "sum")

    # Total area of transboundary buffers
    tb.area <- as.numeric(st_area(tb.buffers) / 1e6)
    # Total area of protected areas contained within transboundary country buffers
    buffer.exp <- exactextractr::exact_extract(pa.sol.area, tb.buffers, fun = "sum")[[1]]
    # what proportion of protected areas are inside of transboundary regions?
    buffer.prop <- (buffer.exp / ov.area)[[1]]
    # what proportion of transboundary buffer area is protected
    buffer.prop2 <- as.numeric((buffer.exp / tb.area)[[1]])

    o <- data.frame(
      project = paste0(project.name, "__", version),
      scenario = file.name,
      region = "Europe",
      variable = "Transboundary area conserved by the solution",
      variable_variant = "Selected area in region",
      unit = "%",
      value = buffer.prop * 100
    )

    o2 <- o %>% mutate(
      variable_variant = "Transboundary area in region",
      unit = "%",
      value = buffer.prop2 * 100
    )

    tb.europe.existing <- rbind(tb.europe.existing, o, o2)
  }
  toc()

  # Bind results then export
  out.tb.final <- rbind(tb.europe.existing)

  outdir <- paste0("performance_outputs/", version)
  dir.create(outdir, showWarnings = FALSE)

  ofname <- paste0(outdir, "/", project.name, "__", version, "__indicator_transboundaryprotection.csv")
  write.csv(out.tb.final, ofname, row.names = FALSE)

  message("Transboundary protected area calculations done. Don't forget to create/update variable definitions")
}
