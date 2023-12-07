#' Occurrence records for an example species
#'
#' @description A dataset containing geographic coordinates of an example
#' species.
#'
#' @format A data frame with 46 rows and 3 columns.
#' \describe{
#'   \item{species}{character, species scientific name.}
#'   \item{longitude}{numeric, longitude values.}
#'   \item{latitude}{numeric, latitude values.}
#' }
#' @source \url{https://www.gbif.org/}
#'
#' @examples
#' data("records", package = "grinnell")
#' head(records)
"records"



#' Example of variables to represent current climate condition in a region
#'
#' A dataset containing raster layers of climatic variables for a region where
#' simulations can be performed in examples.
#'
#' @format A SpatRaster with 108 rows, 84 columns, 9072 cells, and 6 layers:
#' \describe{
#'   \item{Temperature}{temperature, in Celsius degrees times 10.}
#'   \item{Precipitation}{precipitation, in milimeters.}
#' }
#'
#' @source \url{https://www.worldclim.org/data/index.html}
#'
#' @examples
#' variables <- terra::rast(system.file("extdata/variables.tif",
#'                                        package = "grinnell"))
#'
#' terra::plot(variables[[1]])
#' @name variables
NULL



#' Example of variables to represent current climate condition in a region
#'
#' A dataset containing raster layers of climatic variables representing Last
#' Glacial Maximum conditions for a region where simulations can be performed
#' in examples.
#'
#' @format A SpatRaster with 108 rows, 84 columns, 9072 cells, and 6 layers:
#' \describe{
#'   \item{Temperature}{temperature, in Celsius degrees times 10.}
#'   \item{Precipitation}{precipitation, in milimeters.}
#' }
#'
#' @source \url{https://www.worldclim.org/data/index.html}
#'
#' @examples
#' variables_lgm <- terra::rast(system.file("extdata/variables_lgm.tif",
#'                                            package = "grinnell"))
#'
#' terra::plot(variables_lgm[[1]])
#' @name variables_lgm
NULL




#' Example of layer representing dispersal barriers for a species
#'
#' A raster layer representing dispersal barriers for a species in a region
#' where simulations can be performed in examples.
#'
#' @format A SpatRaster with 108 rows, 84 columns, 9072 cells:
#' \describe{
#'   \item{barrier}{barriers are represented with values of 1}
#' }
#'
#' @examples
#' barrier <- terra::rast(system.file("extdata/barrier.tif",
#'                                       package = "grinnell"))
#'
#' terra::plot(barrier)
#' @name barrier
NULL




#' Example of layer representing environmental suitability for a species
#'
#' A raster layer representing distinct levels of suitability for a species
#' in a region where simulations can be performed in examples.
#'
#' @format A SpatRaster with 108 rows, 84 columns, 9072 cells:
#' \describe{
#'   \item{suitability}{values from low = 0 to high = 1}
#' }
#'
#' @examples
#' suitability <- terra::rast(system.file("extdata/suitability.tif",
#'                                           package = "grinnell"))
#'
#' terra::plot(suitability)
#' @name suitability
NULL




#' Example of layer representing future environmental suitability for a species
#'
#' A raster layer representing distinct levels of future suitability for a
#' species in a region where simulations can be performed in examples.
#'
#' @format A SpatRaster with 108 rows, 84 columns, 9072 cells:
#' \describe{
#'   \item{suitability}{values from low = 0 to high = 1}
#' }
#'
#' @examples
#' suitability_fut <- terra::rast(system.file("extdata/suitability_fut.tif",
#'                                               package = "grinnell"))
#'
#' terra::plot(suitability_fut)
#' @name suitability_fut
NULL
