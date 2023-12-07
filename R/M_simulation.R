#' Simulation of species accessible areas (M) Python version
#'
#' @description M_simulation generates an area that has been potentially
#' accessible to a species based on simulations of dispersal events determined
#' by environmental suitability and user-defined parameters. The dispersal
#' simulation is performed using Python (see details).
#'
#' @param data (character) name of the csv file with all the occurrences used to
#' run the simulation; columns must be: species, longitude, latitude.
#' @param current_variables (character) name of the folder where environmental
#' variables representing current conditions are. T least two variables are
#' needed and they must be in ascii format.
#' @param barriers (character) "optional" name of the raster file representing
#' barriers for the species dispersal. Barriers must have the same projection,
#' format, and extent than variables. The only values allowed in this layer are
#' 1 and NA; 1 = barrier. Default = NULL.
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
#' @param projection_variables (character) name of the folder where environmental
#' variables representing the "Last Glacial Maximum" scenario. Variable names,
#' projection, and extent of these layers must be the same than those in
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
#' @param output_directory (character) name of the output directory to be created
#' in which all results will be written.
#' @param overwrite (logical) whether or not to overwrite the
#' \code{output_directory} if it already exists. Default = FALSE.
#'
#' @return
#' The complete set of results derived from data preparation and the simulation
#' is written in \code{output_directory}. These results include:
#' - occurrences found in suitable areas in the scenario where the simulation
#' started.
#' - a folder containing results from the PCA performed
#' - a folder containing results from the preparation of suitability layer(s)
#' - accessible areas as raster layers (value 1 = accessed)
#' - other raster layers representing statistics of accessibility:
#' mean and variance
#' - accessible areas as a shapefile (only accessed areas)
#' - a plot representing the accessible areas and the occurrences
#' - a simple report from the simulation process
#'
#' @usage
#' M_simulation(data, current_variables, barriers = NULL, project = FALSE,
#'              projection_variables, scale = TRUE, center = TRUE,
#'              dispersal_kernel = "normal", kernel_spread = 1,
#'              max_dispersers = 4, suitability_threshold = 5,
#'              replicates = 10, dispersal_events = 20,
#'              access_threshold = 5, simulation_period = 50,
#'              stable_lgm = 7, transition_to_lgm = 100, lgm_to_current = 7,
#'              stable_current = 13, scenario_span = 1, output_directory,
#'              overwrite = FALSE)
#'
#' @export
#' @importFrom terra writeRaster rast trim as.polygons crs vect plot writeVector
#' @importFrom ellipse ellipse
#'
#' @details
#' An external dependency for this function is Python >= 3.6 and some of its
#' libraries: os, numpy, linecache, csv, copy, math, time. We recommend to
#' install Python 3 via Anaconda which will include all the required libraries.
#'
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
#' \dontrun{
#' # preparing data and directories for examples
#' ## directories
#' tempdir <- file.path(tempdir(), "msim")
#' dir.create(tempdir)
#'
#' cvariables <- paste0(tempdir, "/variables")
#' dir.create(cvariables)
#'
#' lgmvariables <- paste0(tempdir, "/LGM")
#' dir.create(lgmvariables)
#'
#' ## data
#' data("records", package = "grinnell")
#' variables <- terra::rast(system.file("extdata/variables.tif",
#'                                      package = "grinnell"))
#' variables_lgm <- terra::rast(system.file("extdata/variables_lgm.tif",
#'                                          package = "grinnell"))
#' names(variables_lgm) <- names(variables)
#' barrier <- terra::rast(system.file("extdata/barrier.tif",
#'                                    package = "grinnell"))
#'
#' ## writing data in temporal directories
#' occ <- paste0(tempdir, "/records1.csv")
#' write.csv(records, occ, row.names = FALSE)
#'
#' barr <- paste0(tempdir, "/barrier1.asc")
#' terra::writeRaster(barrier, filename = barr)
#'
#' vnam <- paste0(cvariables, "/var_", 1:6, ".asc")
#' terra::writeRaster(variables, filename = vnam)
#'
#' vnam <- paste0(lgmvariables, "/var_" 1:6, ".asc")
#' terra::writeRaster(variables_lgm, filename = vnam)
#'
#' odir1 <- paste0(tempdir, "/eg_msim1")
#' odir2 <- paste0(tempdir, "/eg_msim2")
#' odir3 <- paste0(tempdir, "/eg_msim3")
#'
#' # simulations
#' ## example in current scenario
#' M_simulation(data = occ, current_variables = cvariables,
#'              max_dispersers = 2, replicates = 3, dispersal_events = 5,
#'              output_directory = odir1)
#'
#' ## example under changing climatic conditions (starting from the past)
#' M_simulation(data = occ, current_variables = cvariables,
#'              project = TRUE, projection_variables = lgmvariables,
#'              kernel_spread = 2, max_dispersers = 2,
#'              replicates = 3, dispersal_events = 25,
#'              simulation_period = 25, stable_lgm = 7,
#'              transition_to_lgm = 3, lgm_to_current = 3,
#'              stable_current = 7, scenario_span = 3,
#'              output_directory = odir2)
#'
#' ## example under changing conditions, considering dispersal barriers
#' M_simulation(data = occ, current_variables = cvariables,
#'              barriers = barr, project = TRUE,
#'              projection_variables = lgmvariables,
#'              kernel_spread = 2, max_dispersers = 2,
#'              replicates = 3, dispersal_events = 25,
#'              simulation_period = 25, stable_lgm = 7,
#'              transition_to_lgm = 3, lgm_to_current = 3,
#'              stable_current = 7, scenario_span = 3,
#'              output_directory = odir3)
#' }

M_simulation <- function(data, current_variables, barriers = NULL, project = FALSE,
                         projection_variables, scale = TRUE, center = TRUE,
                         dispersal_kernel = "normal", kernel_spread = 1,
                         max_dispersers = 4, suitability_threshold = 5,
                         replicates = 10, dispersal_events = 20,
                         access_threshold = 5, simulation_period = 50,
                         stable_lgm = 7, transition_to_lgm = 100, lgm_to_current = 7,
                         stable_current = 13, scenario_span = 1,
                         output_directory, overwrite = FALSE) {

  # --------
  # testing for initial requirements
  ## existing directpry
  if (missing(output_directory)) {
    stop("Argument 'output_directory' must be defined")
  }
  if (overwrite == FALSE & dir.exists(output_directory)) {
    stop("'output_directory' already exists, to replace it use overwrite = TRUE")
  }
  if (overwrite == TRUE & dir.exists(output_directory)) {
    unlink(x = output_directory, recursive = TRUE, force = TRUE)
  }

  if (missing(data)) {
    stop("Argument 'data' must be defined")
  }

  ## python 3.6 or superior
  message("Checking dependencies")
  if (.Platform$OS.type == "unix") {
    py <- system("python3 -V", intern = TRUE)
  } else {
    py <- system("python", intern = TRUE)
  }
  ncl <- gregexpr("Python 3.\\d", py)
  ncla <- regmatches(py, ncl)
  pyver <- unlist(ncla)

  if (length(pyver) == 0) {
    stop("Python version >= 3.6 is needed to run simulations.\n",
         "  Installing Anaconda will install Python and all libraries needed:\n",
         "  https://www.anaconda.com/products/individual#Downloads")
  }

  version <- as.numeric(strsplit(pyver, " ")[[1]][2])

  if (version < 3.6) {
    stop("Python version >= 3.6 is needed to run simulations.\n",
         "  Installing Anaconda will install Python and all libraries needed:\n",
         "  https://www.anaconda.com/products/individual#Downloads")
  }

  ## other arguments
  if (missing(current_variables)) {
    stop("Argument 'current_variables' must be defined")
  }
  var <- list.files(current_variables, pattern = ".asc$", full.names = TRUE)
  n <- length(var)
  varss <- terra::rast(var)

  if (n < 2) {
    stop("At least 2 variables are needed in 'current_variables'")
  }
  if (project == TRUE) {
    if (missing(projection_variables)) {
      stop("If 'project' = TRUE, argument 'projection_variables' must be defined")
    }
    pvar <- list.files(projection_variables, pattern = ".asc$", full.names = TRUE)
    lgmss <- terra::rast(pvar)
    if (!all(terra:ext(varss) == terra::ext(lgmss))) {
      stop("'projection_variables' and 'current_variables' must have the same extent")
    }
    if (!all(names(varss) == names(lgmss))) {
      stop("Variable names in 'projection_variables' and 'current_variables' must bee the same")
    }
  }

  # --------
  # preparation of data in R
  message("Preparing data to run simulation...")

  # output directory and slash type for directories
  dir.create(output_directory)

  # occurrences
  data <- read.csv(data)
  sp_nam <- as.character(data[1, 1])
  data <- data[, 2:3]

  # --------
  # principal components; interpolations if needed; and suitability layer(s)
  ## folder for suitability layers
  suit_fol <- paste0(output_directory, "/Initial_suitability")
  dir.create(suit_fol)

  npcs <- ifelse(n > 3, 3, n) # number of pcs
  opca_fol <- paste0(output_directory, "/PCA_results") # PCA folder

  if (project == FALSE) {
    ## pca
    variables <- pca_raster(variables = varss, scale = scale, center = center,
                            n_pcs = npcs, project = project,
                            write_to_directory = TRUE,
                            output_directory = opca_fol)[[2]]

    ## suitability
    suit_mod <- ellipsoid_suitability(data, variables, suitability_threshold)
    suit_layer <- suit_mod[[5]]

    ## barrier consideration
    if (!is.null(barriers)) {
      message("  Correcting suitability layer using barriers")
      barriers <- terra::rast(barriers)
      barr <- is.na(barriers)
      suit_layer <- suit_layer * barr
    }

    ### write suitability layer
    emodfile <- paste0(suit_fol, "/Ellipsoid_metadata")
    write_ellmeta(suit_mod, emodfile)

    s_name <- paste0(suit_fol, "/suitability.asc")
    terra::writeRaster(suit_layer, s_name)

    suit_name <- paste0("\"", normalizePath(s_name), "\"")
    suit_name <- gsub("\\\\", "/", suit_name)

  } else {
    pcs <- pca_raster(variables = varss, scale = scale, center = center,
                      n_pcs = npcs, project = project,
                      projection_variables = lgmss, return_projection = TRUE,
                      write_to_directory = TRUE,
                      output_directory = opca_fol)[c(2, 3)]

    variables <- pcs[[1]]
    lgm <- pcs[[2]][[1]]

    ## initial suitabilities
    suit_mod <- ellipsoid_suitability(data, variables, suitability_threshold,
                                      project = TRUE, projection_variables = lgm)
    suit_layer <- suit_mod[[5]][[1]]
    suit_lgm <- suit_mod[[5]][[2]]

    ## barrier consideration
    if (!is.null(barriers)) {
      message("  Suitability layers will be corrected using barriers")
      barriers <- terra::rast(barriers)
      barr <- is.na(barriers)
      suit_layer <- suit_layer * barr
      suit_lgm <- suit_lgm * barr
    } else {
      barr <- barriers
    }

    ### write suitability layer current and lgm
    emodfile <- paste0(suit_fol, "/Ellipsoid_metadata")
    write_ellmeta(suit_mod, emodfile)

    s_name <- paste0(suit_fol, "/suitability_current.asc")
    terra::writeRaster(suit_layer, s_name)

    l_name <- paste0(suit_fol, "/suitability_lgm.asc")
    terra::writeRaster(suit_lgm, l_name)

    ### names for python
    sp_name <- paste0("\"", normalizePath(s_name), "\"")
    sp_name <- gsub("\\\\", "/", sp_name)
    lp_name <- paste0("\"", normalizePath(l_name), "\"")
    lp_name <- gsub("\\\\", "/", lp_name)

    ## preparing interpolation cycles
    int_vals <- interpolation_values(simulation_period, transition_to_lgm,
                                     stable_lgm, lgm_to_current,
                                     stable_current, scenario_span)

    ## interpolations and suitability layer projections
    suit_name <- interpolation(suit_mod, suitability_threshold, int_vals,
                               barr, variables, lgm, sp_name, lp_name,
                               out_format = "ascii", suit_fol)
    suit_name <- gsub("\"", "", suit_name)
    suit_name <- gsub("\\\\", "/", suit_name)
    suit_name <- paste0("\"", suit_name, "\"")
  }

  # --------
  # occurrences in suitable areas
  occ_suit <- as.matrix(suit_mod[[1]][, 1:2])
  suit_bar <- terra::extract(terra::rast(gsub("\"", "", suit_name[1])),
                              occ_suit)
  occ_suit <- occ_suit[suit_bar > 0, ]

  ## records
  oca <- data.frame(Species = sp_nam, occ_suit)
  oca_nam <- paste0(output_directory, "/occ_simulation.csv")
  write.csv(oca, oca_nam, row.names = FALSE)

  occ_name <- paste0("\"", normalizePath(oca_nam),"\"")
  occ_name <- gsub("\\\\", "/", occ_name)

  # --------
  # figure of niche centroid model in E space
  save_nicheplot(suit_mod, suitability_threshold, variables,
                 size_proportion = 0.55, suit_fol)

  # --------
  # python script preparation and execution
  out_dir <- paste0("\"", normalizePath(output_directory), "\"")
  out_dir <- gsub("\\\\", "/", out_dir)

  ## script
  dispersal_simulation(data = occ_name, suit_layers = suit_name,
                       dispersal_kernel = dispersal_kernel,
                       kernel_spread = kernel_spread,
                       max_dispersers = max_dispersers,
                       replicates = replicates,
                       dispersal_events = dispersal_events,
                       access_threshold = access_threshold,
                       output_directory = out_dir)

  ## execution
  message("Running simulation...")
  py.script_r <- normalizePath(paste0(output_directory, "/M_simulation.py"))
  if (.Platform$OS.type == "unix") {
    system(paste("python3", py.script_r))
  } else {
    system(paste("python", py.script_r))
  }

  # --------
  # preparing, writing, and outputs
  ## M
  message("Preparing M as a spatial polygon...")
  afile <- paste0("A_S", length(suit_name),".asc")
  m <- M_preparation(output_directory, A_name = afile, raster_format = "ascii")
  m_poly <- m[[2]]
  m <- m[[1]]

  ## plot and save figure of M in geographic space
  save_Mplot(suit_mod, suit_layer, m_poly, size_proportion = 0.55,
             output_directory)

  message("\nM simulation finished\nCheck results in:  ",
          gsub("\"", "", out_dir))
}
