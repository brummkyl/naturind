#' Transboundary buffer file creation
#'
#' @param regions An sf vector layer containing individual regions (e.g., countries)
#' @param buffer.width An integer defining the width of buffers (in meters) in either direction, i.e., the total buffer width will be 2x the defined 'buffer.width'
#' @param progress A logical of whether progress messages should be printed
#' @param plot A logical of whether the resulting layer should be plotted
#' @param save A logical of whether the resulting layer should be saved
#' @param outdir If 'save = TRUE', an absolute file path of where to save the resulting layer. Note: 'outdir' must end with a forward slash, e.g., ("C:/.../folder/")
#'
#' @returns A .gpkg layer (sf object) containing a single feature, representing the transboundary area between regions
#' @export
#'
#' @examples
#' \dontrun{
#' tb.buffers <- create.tb.buffers(regions,
#'   buffer.width = 10000, progress = TRUE,
#'   plot = TRUE, save = FALSE, outdir = NULL
#' )
#' }
create.tb.buffers <- function(regions, buffer.width = 10000, progress = TRUE, plot = TRUE, save = FALSE, outdir = NULL) {
  countries_sf <- regions
  countries_sf <- suppressWarnings(st_make_valid(countries_sf))
  shared_border_buffers <- list()

  # width = buffer distance (in meters), extending from each side of a region border
  # the joint width of transboundary buffers corresponds to width x 2
  width <- buffer.width

  for (i in 1:nrow(countries_sf)) {
    if (isTRUE(progress)) {
      message("Processing buffer ", i, " of ", nrow(countries_sf))
    }

    country_i <- countries_sf[i, ]

    for (j in (i + 1):nrow(countries_sf)) {
      country_j <- countries_sf[j, ]

      if (st_touches(country_i, country_j, sparse = FALSE)[1, 1]) {
        buffer_i <- suppressWarnings(st_buffer(country_i, dist = width))
        buffer_j <- suppressWarnings(st_buffer(country_j, dist = width))

        shared_zone <- suppressWarnings(st_intersection(buffer_i, buffer_j))

        if (nrow(shared_zone) > 0) {
          shared_border_buffers[[length(shared_border_buffers) + 1]] <- shared_zone
        }
      }
    }
  }

  # Combine all shared border buffers
  if (length(shared_border_buffers) > 0) {
    if (isTRUE(progress)) {
      message("Combining shared border buffers...")
    }

    tb.buffers <- do.call(rbind, shared_border_buffers)
    clipped <- suppressWarnings(sf::st_intersection(tb.buffers, countries_sf))
    one_geom <- suppressWarnings(sf::st_union(sf::st_geometry(clipped)))

    tb.buffers <- suppressWarnings(sf::st_sf(geometry = sf::st_sfc(one_geom, crs = sf::st_crs(countries_sf))))

    stopifnot(nrow(tb.buffers) == 1L)

    if (isTRUE(plot)) {
      plot(st_geometry(countries_sf), col = "lightgray", border = "lightgray")
      plot(st_geometry(tb.buffers), col = "red", add = TRUE)
    }

    message("Done! Shared border buffers stored as 'tb.buffers'")

    if (isTRUE(save) && !is.null(outdir)) {
      st_write(tb.buffers, paste0(outdir, width / 1000, "km_TBbuffers.gpkg"), append = FALSE)
    }
  } else {
    message("No shared borders found.")
  }
  return(tb.buffers)
}
