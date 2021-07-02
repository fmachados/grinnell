#' Simulation of species dispersal to future
#'
#' @description Simulation of species dispersal to identify areas
#' that could be accessed and colonized based on environmental
#' suitability and user-defined dispersal parameters.
#'
#' @param suit_layer (character) name of (local) raster layer containing values
#' of suitability for the species of interest in the study region over which the
#' simulation will run. The name of the layer name should include parent
#' directory if needed.
#' @param data (optional) data.frame containing geographic coordinates of
#' occurrences of the species of interest. Columns must be: species, longitude,
#' latitude, in that order.
#' @param suit_forward (optional) name of (local) raster layer(s) containing
#' values of suitability for the species of interest in the study region over
#' which the simulation will run. If more than one, layer names must be
#' ordered from first to last scenario. These layers must have the same
#' extent, resolution, number of cells, and projection than \code{suit_layer}.
#' Layer names should include parent directory if needed.
#' @param barriers (optional) RasterLayer representing dispersal barriers for
#' the species. This layer must have the same extent and projection than
#' \code{suit_layers}. The only values allowed in this layer are 1 and NA;
#' 1 = barrier. Default = NULL.
#' @param starting_proportion (numeric) proportion of \code{data} to be used as
#' starting points for the simulation. Default = 0.5. All data is used if a
#' value of 1 is defined.
#' @param proportion_to_disperse (numeric) proportion of colonized cells from
#' which dispersers will start a new dispersal event; default = 1.
#' @param dispersal_kernel (character) dispersal kernel (dispersal function)
#' used to simulate the movement of the species. Options are: "normal",
#' "log_normal". Default = "normal".
#' @param kernel_spread (numeric) standard deviation for the
#' \code{dispersal_kernel}. Default = 1.
#' @param max_dispersers (numeric) maximum number of dispersers that depart from
#' each colonized pixel. Depending on suitability this number will automatically
#' decrease in areas with low suitability values. Default = 4.
#' @param dispersal_events (numeric) number of dispersal events to happen per
#' scenario. Default = 25.
#' @param replicates (numeric) number of times to repeat the simulation
#' per dispersal event, depending on \code{results_by}. Default = 10.
#' @param threshold (numeric) percentage to be considered when excluding
#' accessed or colonized cells with lower values. Default = 5.
#' @param set_seed (numeric) a seed to be used when sampling \code{data}
#' according to \code{starting_proportion}. Default = 1.
#' @param out_format (character) format of raster layers to be written in
#' \code{output_directory}. Options are "ascii", "GTiff", and "EHdr" = bil.
#' Default = "GTiff".
#' @param output_directory (character) name of the output directory where
#' results should be written. If this directory does not exist, it will be
#' created.
#'
#' @export
#' @importFrom grDevices rgb col2rgb
#' @importFrom raster plot
#'
#' @usage
#' forward_simulation(suit_layer, data = NULL, suit_forward = NULL,
#'                    barriers = NULL, starting_proportion = 0.5,
#'                    proportion_to_disperse = 1,
#'                    dispersal_kernel = "normal",
#'                    kernel_spread = 1, max_dispersers = 4,
#'                    dispersal_events = 25, replicates = 10,
#'                    threshold = 5, set_seed = 1,
#'                    out_format = "GTiff", output_directory)
#'
#' @return
#' A list containing:
#' - occurrences found in suitable areas in the scenario where the simulation
#' started.
#' - all scenarios considered for the simulation
#' - a list with the parameters used during the simulation
#' - a RasterLayer with values representing the number of the
#' dispersal event when areas where accessed
#' - a RasterLayer with values representing the number of the
#' dispersal event when areas where colonized
#' - if defined, the RasterLayer used in \code{barriers}, else NULL
#'
#' The complete set of results derived from data preparation and the simulation
#' is written in \code{output_directory}. These results include the ones
#' mentioned above (except barriers), plus:
#' - if needed, a folder containing results from correcting suitability layer(s)
#' with \code{barriers}
#' - other raster layers representing statistics of accessibility:
#' mean and variance
#' - a figure representing accessed and colonized areas per dispersal events,
#' and the occurrences used for simulation
#' - a simple report from the simulation process
#'
#' The number of dispersal events in results is continuous among scenarios. If
#' 10 dispersal events are defined and multiple scenarios exist in
#' \code{suit_layers}, the first dispersal event in the second scenario will be
#' number 11.
#'
#' @examples
#' # data
#' data("records", package = "grinnell")
#' suitability <- system.file("extdata/suitability.tif", package = "grinnell")
#'
#' # simulation current
#' f_s <- forward_simulation(suit_layer = suitability, data = records,
#'                           dispersal_kernel = "normal",
#'                           kernel_spread = 2, max_dispersers = 2,
#'                           dispersal_events = 15, replicates = 3,
#'                           output_directory = file.path(tempdir(), "eg_fsim"))
#'
#' # simulation current and future
#' \donttest{
#' suitf <- system.file("extdata/suitability_fut.tif", package = "grinnell")
#'
#' f_s1 <- forward_simulation(suit_layer = suitability, data = records,
#'                            suit_forward = suitf, dispersal_kernel = "normal",
#'                            kernel_spread = 2, max_dispersers = 2,
#'                            dispersal_events = 15, replicates = 3,
#'                            output_directory = file.path(tempdir(), "eg_fsim1"))
#'
#' # simulation current and future using dispersal barriers
#' barrier <- raster::raster(system.file("extdata/barrier.tif",
#'                                       package = "grinnell"))
#'
#' f_s2 <- forward_simulation(suit_layer = suitability, data = records,
#'                            suit_forward = suitf, barriers = barrier,
#'                            dispersal_kernel = "normal",
#'                            kernel_spread = 2, max_dispersers = 2,
#'                            dispersal_events = 15, replicates = 3,
#'                            output_directory = file.path(tempdir(), "eg_fsim2"))
#' }

forward_simulation <- function(suit_layer, data = NULL, suit_forward = NULL,
                               barriers = NULL, starting_proportion = 0.5,
                               proportion_to_disperse = 1,
                               dispersal_kernel = "normal",
                               kernel_spread = 1, max_dispersers = 4,
                               dispersal_events = 25, replicates = 10,
                               threshold = 5, set_seed = 1,
                               out_format = "GTiff", output_directory) {
  # --------
  # testing for initial requirements
  if (missing(suit_layer)) {
    stop("Argument 'suit_layer' must be defined")
  }
  if (missing(output_directory)) {
    stop("Argument 'output_directory' must be defined")
  }
  if (!is.null(barriers)) {
    suit_lay <- raster::raster(suit_layer)
    if (suit_lay@extent != barriers@extent) {
      stop("'barriers' and 'suit_layer' must have the same extent")
    }
  }

  # --------
  # preparing data
  message("Preparing data to run simulation...")

  ## defining initial points if data = NULL
  if (is.null(data)) {
    if (is.null(barriers)) {
      suit_lay <- raster::raster(suit_layer)
    }
    noz <- which((suit_lay[] > 0))
    noz <- raster::xyFromCell(suit_lay, noz)
    data <- data.frame(species = "Species", longitude = noz[, 1],
                       latitude = noz[, 2])
  }

  ## output directory and raster file type
  dir.create(output_directory)
  out_dir <- normalizePath(output_directory)

  ftype <- rformat_type(out_format)

  ## suitability files and corrections with barriers if needed
  if (!is.null(barriers)) {
    message("\nSuitability layer(s) will be corrected using barriers")
    barr <- is.na(barriers)

    suit_fol <- paste0(out_dir, "/Suitability_barrier_corrected")
    dir.create(suit_fol)

    if (is.null(suit_forward)) {
      suit_lay <- suit_lay * barr

      s_name <- paste0(suit_fol, "/suitability", ftype)
      raster::writeRaster(suit_lay, filename = s_name, format = out_format)

      suit_name <- normalizePath(s_name)

    } else {
      suits <- c(suit_layer, suit_forward)
      len <- length(suits)
      suit_name <- vapply(1:len, FUN.VALUE = character(1), function(x) {
        if (x > 1) {
          suit_lay <- raster::raster(suits[x])
        }
        suit_lay <- suit_lay * barr

        s_name <- paste0(suit_fol, "/suitability", x, ftype)
        raster::writeRaster(suit_lay, filename = s_name, format = out_format)

        normalizePath(s_name)
      })
    }
  } else {
    if (is.null(suit_forward)) {
      suit_name <- normalizePath(suit_layer)
    } else {
      suit_name <- normalizePath(c(suit_layer, suit_forward))
    }
  }

  ## occurrences relevant for simulation
  suit_bar <- raster::extract(raster::raster(suit_name[1]), data[, 2:3])
  data <- data[suit_bar > 0, ]

  ## write relevant records
  oca_nam <- paste0(out_dir, "/occ_simulation.csv")
  write.csv(data, oca_nam, row.names = FALSE)

  # --------
  # running simulation
  message("")
  res <- dispersal_simulationR(data = data, suit_layers = suit_name,
                               starting_proportion = starting_proportion,
                               proportion_to_disperse = proportion_to_disperse,
                               dispersal_kernel = dispersal_kernel,
                               kernel_spread = kernel_spread,
                               max_dispersers = max_dispersers,
                               dispersal_events = dispersal_events,
                               replicates = replicates,
                               threshold = threshold,
                               results_by = "event", set_seed = set_seed,
                               return = "all", write_to_directory = TRUE,
                               raster_format = out_format,
                               output_directory = out_dir)

  # --------
  # preparing, writing, and outputs
  ## plot and save figure of accessed and colonized in geographic space
  save_event_plot(data, res, barriers, size_proportion = 0.55, output_directory)

  message("\nForward simulation finished\nCheck results in:  ", out_dir, "\n")

  # return
  return(list(Simulation_occurrences = data, Simulation_scenarios = suit_name,
              Summary = res$Summary, A_events = res$A_events,
              C_events = res$C_events, Barriers = barriers))
}
