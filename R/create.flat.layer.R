#' Create a Flat Layer of Protected Areas Using 2 Input Vectors
#'
#' @param pa1 A SpatVector layer containing protected areas from source 1
#' @param pa2 A SpatVector layer containing protected areas from source 2
#' @param puid A SpatRaster layer, reflecting the amount of land area per planning unit
#' @param countries An sf vector layer containing individual countries, with names specified in a column named "code"
#' @param n.tiles An integer of the desired grid dimensions to use in processing. Default is set to 10, which generates a 10 x 10 grid (e.g., 100 tiles)
#' @param save A logical of whether the resulting layer should be saved
#' @param outdir If 'save = TRUE', an absolute file path of where to save the resulting layer. Note: 'outdir' must end with a forward slash, e.g., ("C:/.../folder/")
#' @param file.name A character vector of what to name the resulting layer, excluding '.gpkg'
#'
#' @returns A flat .gpkg layer (sf object) containing a single row for each protected area patch. The area of each patch is provided in the column named 'area_km2'
#' @export
#'
#' @examples
#' #'
#' \dontrun{
#' create.flat.layer(pa1, pa2, puid, countries,
#'   n.tiles = 10, save = FALSE,
#'   outdir = NULL, file.name = NULL
#' )
#' }
create.flat.layer <- function(pa1, pa2, puid, countries, n.tiles = 10, save = TRUE, outdir, file.name) {
  nc.pu <- puid
  n2k <- pa1
  natda <- pa2

  n2k <- n2k |> terra::project(crs(nc.pu))
  assertthat::assert_that(crs(nc.pu) == crs(n2k))

  natda <- natda |> terra::project(crs(nc.pu))
  assertthat::assert_that(crs(nc.pu) == crs(natda))

  # Convert SpatVector to sf and extract only the geometries
  n2k_sf <- st_as_sf(n2k)[, "geometry"] # Keep only the geometry column
  natda_sf <- st_as_sf(natda)[, "geometry"]

  # Combine the geometries using rbind
  combined_sf <- rbind(n2k_sf, natda_sf)
  nr <- combined_sf

  # Create spatial grid
  grid <- st_make_grid(nr, n = c(n.tiles, n.tiles))
  length(grid)

  # Turn grid into an sf object
  grid_sf <- st_sf(id = 1:length(grid), geometry = grid)

  # Initialize results list
  chunk_results <- list()

  # Loop through tiles
  for (i in seq_len(nrow(grid_sf))) {
    tile <- grid_sf[i, ]
    message("Processing tile ", i, " of ", nrow(grid_sf))

    # Clip both layers to this tile (bounding box-level intersection)
    base_tile <- suppressWarnings(st_intersection(nr, tile))
    # Union mask for this tile to simplify geometry before difference
    diss_tile <- st_union(base_tile)
    # store tile output in list
    chunk_results[[i]] <- diss_tile
  }

  result <- st_sf(geometry = do.call(c, chunk_results))

  # co.dissolve must exist and be in the same CRS as result
  countries <- sf::st_make_valid(countries) %>% sf::st_transform(sf::st_crs(result))
  co.dissolve <- sf::st_union(countries)

  pa.sol.crop <- result %>%
    sf::st_make_valid() %>%
    sf::st_intersection(co.dissolve) %>%
    sf::st_collection_extract("POLYGON") %>%
    sf::st_union() %>%
    sf::st_cast("POLYGON") %>%
    sf::st_sf(geometry = .)

  pa.sol.crop$area_km2 <- as.numeric(sf::st_area(pa.sol.crop)) / 1e6

  # one-off checksfile.name
  table(sf::st_geometry_type(pa.sol.crop)) # should be POLYGON only
  nrow(pa.sol.crop)
  sum(pa.sol.crop$area_km2)


  if (isTRUE(save)) {
    if (is.null(outdir) || is.null(file.name)) {
      stop("When save = TRUE, provide both outdir and file.name.")
    }

    if (!dir.exists(outdir)) {
      dir.create(outdir, recursive = TRUE)
    }

    sf::st_write(pa.sol.crop, paste0(outdir, file.name, ".gpkg"))
  }

  return(pa.sol.crop)
}
