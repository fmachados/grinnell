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
#' @format A RasterStack with 108 rows, 84 columns, 9072 cells, and 6 layers:
#' \describe{
#'   \item{Temperature}{temperature, in Celsius degrees times 10.}
#'   \item{Precipitation}{precipitation, in milimeters.}
#' }
#'
#' @source \url{https://www.worldclim.org/data/index.html}
#'
#' @examples
#' variables <- raster::stack(system.file("extdata/variables.tif",
#'                                        package = "grinnell"))
#'
#' raster::plot(variables[[1]])
#' @name variables
NULL



#' Example of variables to represent current climate condition in a region
#'
#' A dataset containing raster layers of climatic variables representing Last
#' Glacial Maximum conditions for a region where simulations can be performed
#' in examples.
#'
#' @format A RasterStack with 108 rows, 84 columns, 9072 cells, and 6 layers:
#' \describe{
#'   \item{Temperature}{temperature, in Celsius degrees times 10.}
#'   \item{Precipitation}{precipitation, in milimeters.}
#' }
#'
#' @source \url{https://www.worldclim.org/data/index.html}
#'
#' @examples
#' variables_lgm <- raster::stack(system.file("extdata/variables_lgm.tif",
#'                                            package = "grinnell"))
#'
#' raster::plot(variables_lgm[[1]])
#' @name variables_lgm
NULL




#' Example of layer representing dispersal barriers for a species
#'
#' A raster layer representing dispersal barriers for a species in a region
#' where simulations can be performed in examples.
#'
#' @format A RasterLayer with 108 rows, 84 columns, 9072 cells:
#' \describe{
#'   \item{barrier}{barriers are represented with values of 1}
#' }
#'
#' @examples
#' barrier <- raster::stack(system.file("extdata/barrier.tif",
#'                                      package = "grinnell"))
#'
#' raster::plot(barrier[[1]])
#' @name barrier
NULL
