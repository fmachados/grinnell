#' Simulation of species accessible areas (M) R version
#'
#' @description M_simulation1 generates an area that has been potentially
#' accessible to a species based on simulations of dispersal events determined
#' by environmental suitability and user-defined parameters.
#'
#' @param data data.frame with occurrence records for the species of interest to
#' be used to run the simulation; columns must be: species, longitude, latitude.
#' @param current_variables RasterStack of environmental variables representing
#' "current" conditions (interglacial). Recommended projection WGS84 (EPSG:4326).
#' @param starting_proportion (numeric) proportion of \code{data} located in
#' suitable areas to be used as starting points for the simulation.
#' Default = 0.5. All data is used if a value of 1 is defined.
#' @param sampling_rule (character) rule to be used to sample a
#' \code{starting_proportion} of points to run dispersal simulation steps.
#' Options are: "random" and "suitability". Using the option "suitability"
#' prioritizes records with higher suitability values. Default = "random".
#' @param barriers RasterLayer representing dispersal barriers for the species.
#' This layer must have the same extent and projection than
#' \code{current_variables}. The only values allowed in this layer are 1 and NA;
#' 1 = barrier. Default = NULL.
#' @param scale (logical) whether or not to scale variables while performing
#' principal component analyses.
#' @param center (logical) whether or not to center variables while performing
#' principal component analyses.
#' @param project (logical) whether or not to project environmental suitability
#' to past scenarios. The projection is done to the scenario defined by
#' \code{projection_variables} and to any other scenario resulting from
#' interpolations between current and past conditions. If TRUE,
#' arguments \code{projection_variables}, \code{simulation_period},
#' \code{stable_current}, \code{stable_lgm}, \code{transition_to_lgm},
#' \code{lgm_to_current}, and \code{scenario_span} need to be defined.
#' Default = FALSE.
#' @param projection_variables RasterStack of environmental variables
#' representing the "Last Glacial Maximum" scenario. Variable names, projection,
#' and extent of these layers must be the same than those in
#' \code{current_variables}.
#' @param dispersal_kernel (character) dispersal kernel (dispersal function)
#' used to simulate the movement of the species. Options are: "normal",
#' "log_normal". Default = "normal".
#' @param kernel_spread (numeric) standard deviation for the
#' \code{dispersal_kernel}. Default = 2.
#' @param max_dispersers (numeric) maximum number of dispersers that depart from
#' each colonized pixel. Depending on suitability this number will decrease in
#' less suitable areas.
#' @param suitability_threshold value (percentage) to be used as threshold
#' for suitability; default = 5. Below this value environments are considered
#' unsuitable.
#' @param replicates (numeric) number of times to repeat the simulation
#' per scenario. Default = 10.
#' @param dispersal_events (numeric) number of dispersal events to happen per
#' scenario. Default = 25.
#' @param simulation_period (numeric) time in thousands of years for the
#' complete period of simulation.
#' @param transition_to_lgm (numeric) time in thousands of years for the
#' transition period from current-like (interglacial) to glacial (LGM) climatic
#' conditions.
#' @param stable_lgm (numeric) time in thousands of years for the period when
#' glacial (LGM) conditions are assumed to be relatively stable.
#' @param lgm_to_current (numeric) time in thousands of years for the
#' transition period from glacial (LGM) to current-like (interglacial) climatic
#' conditions.
#' @param stable_current (numeric) time in thousands of years for the period when
#' current-like (interglacial) conditions are assumed to be relatively stable.
#' @param scenario_span (numeric) time in thousands of years that have to pass
#' for changing the scenario. Default = 1 (one thousand years).
#' @param access_threshold (numeric) percentage of frequency values to be
#' considered as highly unlikely to have been visited during process of
#' dispersal. Default = 5.
#' @param out_format (character) format of raster layers to be written in
#' \code{output_directory}. Options are "ascii", "GTiff", and "EHdr" = bil.
#' Default = "GTiff".
#' @param set_seed (numeric) a seed to be used when sampling \code{data}
#' according to \code{starting_proportion}. Default = 1.
#' @param write_all_scenarios (logical) whether or not to write results
#' for all scenarios. The default, FALSE, writes only the final results.
#' @param output_directory (character) name of the output directory to be created
#' in which all results will be written.
#'
#' @return
#' A list containing:
#' - occurrences found in suitable areas in the scenario where the simulation
#' started.
#' - all and ordered scenarios considered for the simulation
#' - a list with the parameters used during the simulation
#' - accessible areas as a RasterLayer (value 1 = accessed)
#' - accessible areas as a SpatialPolygons* object (only accessed areas)
#' - a RasterLayer representing mean values of accessibility frequency among
#' replicates
#' - a RasterLayer representing variance among values of accessibility frequency
#' of all replicates
#' - if defined, the RasterLayer used in \code{barriers}, else NULL
#'
#' The complete set of results derived from data preparation and the simulation
#' is written in \code{output_directory}. These results include the ones
#' mentioned above (except barriers), plus:
#' - a folder containing results from the PCA performed
#' - a folder containing results from the preparation of suitability layer(s)
#' - other raster layers representing statistics of accessibility:
#' mean and variance
#' - a plot representing the accessible areas and the occurrences
#' - a simple report from the simulation process
#'
#' @export
#' @importFrom raster writeRaster stack mask crop trim rasterToPolygons image
#' @importFrom sp proj4string CRS plot
#' @importFrom rgdal writeOGR
#' @importFrom ellipse ellipse
#' @importFrom graphics box legend lines par points
#' @importFrom grDevices dev.off png terrain.colors
#'
#' @usage
#' M_simulation1(data, current_variables, starting_proportion = 0.5,
#'               sampling_rule = "random", barriers = NULL, scale = TRUE,
#'               center = TRUE, project = FALSE, projection_variables = NULL,
#'               dispersal_kernel = "normal", kernel_spread = 1,
#'               max_dispersers = 4, suitability_threshold = 5,
#'               replicates = 10, dispersal_events = 25,
#'               access_threshold = 5, simulation_period = 50,
#'               stable_lgm = 7, transition_to_lgm = 100,
#'               lgm_to_current = 7, stable_current = 13,
#'               scenario_span = 1, out_format = "GTiff", set_seed = 1,
#'               write_all_scenarios = FALSE, output_directory)
#'
#' @details
#' A principal component analysis is performed with \code{current_variables}.
#' Then the three first principal components are used to calculate the
#' suitability layer used in dispersal simulations. Values of suitability are
#' derived from an ellipsoid envelope model created based on occurrence records
#' and principal components. The ellipsoid model is used because it is a simple
#' yet reliable representation of a species ecological niche that does not
#' require a background or pseudo-absences.
#'
#' If \code{barriers} are used, suitability values in the areas where barriers
#' exist become zero. This is, populations cannot establish there and dispersal
#' will be truncated unless dispersal abilities defined by arguments
#' \code{dispersal_kernel} and \code{kernel_spread}, allow the species to
#' overpass the barriers.
#'
#' If \code{project} = TRUE, the simulation will run on a set of scenarios
#' representing glacial-interglacial climate conditions. This set of scenarios
#' are constructed based on interpolations between environmental conditions in
#' \code{current_variables} and \code{projection_variables}. The later set of
#' variables must represent Last Glacial Maximum conditions. Interpolations
#' are linear and depend on the distance between the two initial set of
#' conditions and other parameter defined in \code{simulation_period},
#' \code{stable_current}, \code{stable_lgm}, \code{transition_to_lgm},
#' \code{lgm_to_current}, and \code{scenario_span}.
#'
#' @examples
#' # data
#' data("records", package = "grinnell")
#' variables <- raster::stack(system.file("extdata/variables.tif",
#'                                        package = "grinnell"))
#' # example in current scenario
#' m <- M_simulation1(data = records, current_variables = variables,
#'                    max_dispersers = 2, replicates = 3, dispersal_events = 5,
#'                    output_directory = file.path(tempdir(), "eg_Msim1"))
#'
#' # example under changing climatic conditions (starting from the past)
#' \donttest{
#' variables_lgm <- raster::stack(system.file("extdata/variables_lgm.tif",
#'                                            package = "grinnell"))
#' names(variables_lgm) <- names(variables)
#'
#' m_p <- M_simulation1(data = records, current_variables = variables,
#'                      project = TRUE, projection_variables = variables_lgm,
#'                      kernel_spread = 2, max_dispersers = 2,
#'                      replicates = 3, dispersal_events = 25,
#'                      simulation_period = 25, stable_lgm = 7,
#'                      transition_to_lgm = 3, lgm_to_current = 3,
#'                      stable_current = 7, scenario_span = 3,
#'                      output_directory = file.path(tempdir(), "eg_Msim1_p"))
#'
#' # example under changing conditions, considering dispersal barriers
#' barrier <- raster::raster(system.file("extdata/barrier.tif",
#'                                       package = "grinnell"))
#'
#' m_pb <- M_simulation1(data = records, current_variables = variables,
#'                       barriers = barrier, project = TRUE,
#'                       projection_variables = variables_lgm,
#'                       kernel_spread = 2, max_dispersers = 2,
#'                       replicates = 3, dispersal_events = 25,
#'                       simulation_period = 25, stable_lgm = 7,
#'                       transition_to_lgm = 3, lgm_to_current = 3,
#'                       stable_current = 7, scenario_span = 3,
#'                       output_directory = file.path(tempdir(), "eg_Msim1_pb"))
#' }

M_simulation1 <- function(data, current_variables, starting_proportion = 0.5,
                          sampling_rule = "random", barriers = NULL,
                          scale = TRUE, center = TRUE, project = FALSE,
                          projection_variables = NULL,
                          dispersal_kernel = "normal", kernel_spread = 1,
                          max_dispersers = 4, suitability_threshold = 5,
                          replicates = 10, dispersal_events = 25,
                          access_threshold = 5, simulation_period = 50,
                          stable_lgm = 7, transition_to_lgm = 100,
                          lgm_to_current = 7, stable_current = 13,
                          scenario_span = 1, out_format = "GTiff", set_seed = 1,
                          write_all_scenarios = FALSE, output_directory) {

  # --------
  # testing for initial requirements
  if (missing(data)) {
    stop("Argument 'data' must be defined")
  }
  if (missing(current_variables)) {
    stop("Argument 'current_variables' must be defined")
  }
  if (missing(output_directory)) {
    stop("Argument 'output_directory' must be defined")
  }
  if (!sampling_rule %in% c("random", "suitability")) {
    stop("Argument 'sampling_rule' is not valid")
  }

  n <- raster::nlayers(current_variables)
  if (n < 2) {
    stop("At least 2 variables are needed in 'current_variables'")
  }
  if (project == TRUE) {
    if (missing(projection_variables)) {
      stop("If 'project' = TRUE, argument 'projection_variables' must be defined")
    }
    if (!all(current_variables@extent == projection_variables@extent)) {
      stop("'projection_variables' and 'current_variables' must have the same extent")
    }
    if (!all(names(current_variables) == names(projection_variables))) {
      stop("Variable names in 'projection_variables' and 'current_variables' must bee the same")
    }
  }
  if (!is.null(barriers)) {
    if (current_variables@extent != barriers@extent) {
      stop("'barriers' and 'current_variables' must have the same extent")
    }
  }

  # --------
  # preparation of data in R
  message("Preparing data to run simulation...")

  # output directory and slash type of layers
  dir.create(output_directory)
  ftype <- rformat_type(out_format)

  # occurrences
  sp_nam <- as.character(data[1, 1])
  data <- data[, 2:3]

  # --------
  # principal components; interpolations if needed; and suitability layer(s)
  ## folder for suitability layers
  suit_fol <- paste0(output_directory, "/Suitability_results")
  dir.create(suit_fol)

  npcs <- ifelse(n > 3, 3, n) # number of pcs
  opca_fol <- paste0(output_directory, "/PCA_results") # PCA folder

  if (project == FALSE) {
    ## pca
    current_variables <- pca_raster(variables = current_variables, scale = scale,
                                    center = center, n_pcs = npcs,
                                    project = project, write_to_directory = TRUE,
                                    output_directory = opca_fol)[[2]]

    ## suitability
    suit_mod <- ellipsoid_suitability(data, current_variables,
                                      suitability_threshold)
    suit_layer <- suit_mod[[5]]

    ## barrier consideration
    if (!is.null(barriers)) {
      message("\nSuitability layer will be corrected using barriers")
      barr <- is.na(barriers)
      suit_layer <- suit_layer * barr
    }

    ## write suitability layer
    emodfile <- paste0(suit_fol, "/Ellipsoid_metadata")
    write_ellmeta(suit_mod, emodfile)

    s_name <- paste0(suit_fol, "/suitability", ftype)
    raster::writeRaster(suit_layer, filename = s_name, format = out_format)

    suit_name <- normalizePath(s_name)

  } else {
    pcs <- pca_raster(variables = current_variables, scale = scale,
                      center = center, n_pcs = npcs, project = project,
                      projection_variables = projection_variables,
                      return_projection = TRUE, write_to_directory = TRUE,
                      out_format = out_format,
                      output_directory = opca_fol)[c(2, 3)]

    current_variables <- pcs[[1]]
    projection_variables <- pcs[[2]][[1]]

    ## initial suitability
    suit_mod <- ellipsoid_suitability(data, current_variables,
                                      suitability_threshold, project = project,
                                      projection_variables = projection_variables)

    suit_layer <- suit_mod[[5]][[1]]
    suit_lgm <- suit_mod[[5]][[2]]

    ## barrier consideration
    if (!is.null(barriers)) {
      message("\nSuitability layers will be corrected using barriers")
      barr <- is.na(barriers)
      suit_layer <- suit_layer * barr
      suit_lgm <- suit_lgm * barr
    } else {
      barr <- barriers
    }

    ## write suitability layer current and lgm
    emodfile <- paste0(suit_fol, "/Ellipsoid_metadata")
    write_ellmeta(suit_mod, emodfile)

    s_name <- paste0(suit_fol, "/suitability_current", ftype)
    raster::writeRaster(suit_layer, filename = s_name, format = out_format)

    l_name <- paste0(suit_fol, "/suitability_lgm", ftype)
    raster::writeRaster(suit_lgm, filename = l_name, format = out_format)

    ## names for later
    sp_name <- normalizePath(s_name)
    lp_name <- normalizePath(l_name)

    ## preparing interpolation cycles
    message("\nPreparing interpolations:")
    int_vals <- interpolation_values(simulation_period, transition_to_lgm,
                                     stable_lgm, lgm_to_current,
                                     stable_current, scenario_span)

    ## interpolations and suitability layer projections
    suit_name <- interpolation(suit_mod, suitability_threshold, int_vals,
                               barr, current_variables,
                               projection_variables, sp_name, lp_name,
                               out_format, suit_fol)
  }

  # --------
  # occurrences in suitable areas, starting scenario of simulation
  occ_suit <- suit_mod[[1]][, 1:2]
  suit_bar <- raster::extract(raster::raster(suit_name[1]), occ_suit)
  occ_suit <- occ_suit[suit_bar > 0, ]

  ## records
  oca <- data.frame(Species = sp_nam, occ_suit)
  oca_nam <- paste0(output_directory, "/occ_simulation.csv")
  write.csv(oca, oca_nam, row.names = FALSE)

  # --------
  # figure of niche centroid model in E space
  save_nicheplot(suit_mod, suitability_threshold, current_variables,
                 size_proportion = 0.55, suit_fol)

  # --------
  # running simulation
  out_dir <- normalizePath(output_directory)

  ## script
  message("")
  res <- dispersal_simulationR(data = oca, suit_layers = suit_name,
                               starting_proportion = starting_proportion,
                               proportion_to_disperse = 1,
                               sampling_rule = sampling_rule,
                               dispersal_kernel = dispersal_kernel,
                               kernel_spread = kernel_spread,
                               max_dispersers = max_dispersers,
                               dispersal_events = dispersal_events,
                               replicates = replicates,
                               threshold = access_threshold,
                               results_by = "scenario", set_seed = set_seed,
                               return = "accessed", write_to_directory = TRUE,
                               write_all = write_all_scenarios,
                               raster_format = out_format,
                               output_directory = out_dir)

  # --------
  # preparing, writing, and outputs
  ## M
  message("\nPreparing M as a spatial polygon...")
  m <- M_preparation(output_directory, res$A, raster_format = out_format)
  m_poly <- m[[2]]
  m <- m[[1]]

  ## plot and save figure of M in geographic space
  save_Mplot(suit_mod, suit_layer, m_poly, size_proportion = 0.55,
             output_directory)

  message("\nM simulation finished\nCheck results in:  ", out_dir, "\n")

  # return
  return(list(Simulation_occurrences = oca, Simulation_scenarios = suit_name,
              Summary = res$Summary, A_raster = m, A_polygon = m_poly,
              A_mean = res$A_mean, A_var = res$A_var, Barriers = barriers))
}
