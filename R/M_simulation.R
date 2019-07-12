#' Simulation of species accessible areas (M in BAM)
#'
#' @description M_simulation generates an area that has been potentially
#' accessible for a given species during relevant periods of time, derived from
#' a dispersal simulation process.
#'
#' @param data (character) name of the csv file with all the occurrences used to
#' run the simulation; columns must be: species, longitude, latitude.
#' @param current_variables (character) name of the folder where environmental
#' variables representing current conditions are. Variables must be in ascii
#' format and at least two. If \code{project} = TRUE, variables in
#' \code{projection_variables} must have the same extent than variables in this
#' folder.
#' @param barriers (character) "optional" name of the folder where the barriers
#' for the species dispersal are. Barriers must be the same projection than
#' variables (WGS84 and not planar projection is recommended). This barriers can
#' be in two formats: ascii format with the same resolution than variables,
#' or shapefile.
#' @param barrier_format (character) if \code{barriers} are defined, format of
#' the layers to be used. Options are: "ascii" and "shp".
#' @param dispersal_kernel (character) dispersal kernel (dispersal function)
#' used to simulate the movement of the species. Options are: "normal",
#' "log_normal". Default = "normal".
#' @param kernel_spread (numeric) standar deviation for the
#' \code{dispersal_kernel}. Default = 2.
#' @param max_dispersers (numeric) maximum number of dispersors that depart from
#' each colonized pixel. Depending on suitability this number will decrease in
#' less suitable areas.
#' @param suitability_threshold (numeric) percentage of omission error to define
#' what is suitable and whai is not. Default = 5.
#' @param replicates (numeric) number of times that the complete simulation
#' process will be repeated.
#' @param dispersal_events (numeric) number of dispersal events to take place
#' during the entire period of simulation if \code{project} = FALSE, or each
#' \code{scenario_span} if \code{project} = TRUE.
#' @param project (logical) whether or not to perform the simulation using
#' suitabilities that change among scenarios that are interpolated based in
#' current and past conditions. If TRUE, arguments \code{projection_variables},
#' \code{simulation_period}, \code{stable_current}, \code{stable_lgm},
#' \code{transition_to_lgm}, \code{lgm_to_current}, and \code{scenario_span},
#' need to be defined. Default = FALSE.
#' @param projection_variables (character) name of the folder where environmental
#' variables for representing the Last Glacial Maximum scenario are. Variables
#' must be in ascii format and they must correspond with variables in
#' \code{current_variables} (i.e., their name and extent must be the same but
#' conditions should represent the LGM).
#' @param simulation_period (numeric) time in thousands of years for the complete
#' simulation.
#' @param stable_current (numeric) time in thousands of years in which climate
#' is assumed to be relatively stable in the current scenario.
#' @param stable_lgm (numeric) time in thousands of years in which climate is
#' assumed to be relatively stable in the last maximum glacial scenario.
#' @param transition_to_lgm (numeric)
#' @param lgm_to_current (numeric)
#' @param scenario_span (numeric) time in thousands of years that have to pass
#' for changing the scenario in which the simulation will be performed.
#' Default = 1 (one thousend years).
#' @param access_threshold (numeric) percentage of frequency values to be
#' considered as unprobably visited during the process of dispersal. Default = 5.
#' @param write_replicates (logical) whether or not write M ascii files resulted
#' of each replicate. Default = FALSE
#' @param output_directory (character) name of the output directory to be created
#' in which subdirectories containing all results will be written. Default =
#' "M_simulation".
#' @param plot (logical) whether or not to plot the species' niche ellipsoid and
#' the resultant M on the environmental suitability based on the species' niche,
#' the species occurrences are included as well. Default = TRUE.
#' @param mask_variables (logical) whether or not to mask the variables to the
#' simulated accessible area (M). Default = FALSE
#' @param directory_masked (character) name of the folder to be created to save
#' the masked variables. Default = "M_variables".
#'
#' @return Folder \code{output_directory} conatining the results of the simulation.
#' These results include: A plot of the M and the occurrences on an environmental
#' layer map that will return to a graphic device if \code{plot.ellipse} = TRUE;
#' the M as a shapefile and as a raster layer in ascii format; a folder with ....
#'
#' If \code{mask_variables} = TRUE, the environmental variables will be masked
#' to the M and the masked layers will be written in the \code{directory_masked}
#' folder.
#'
#' @details
#' When more than three variables are present in the \code{var.folder} directory,
#' a principal component analysis of them is performed. The the three first
#' principal components are used then, instead of the complete set of variables,
#' when creating the suitability layer over which the dispersal simulation is
#' performed.
#'
#' If barriers are used, the suitability values in the areas where barriers exist
#' turn into cero. This is, populations cannot stablish there and dispersal will
#' be truncated unless the dispersal kernel defined by argumens
#' \code{dispersal_kernel} and \code{kernel_spread}, allow the species to overpass
#' the barriers.
#'
#' @export

M_simulation <- function(data, current_variables, project = FALSE,
                         projection_variables, scale = TRUE,
                         dispersal_kernel = "normal", kernel_spread = 1,
                         max_dispersers = 4, suitability_threshold = 5,
                         replicates = 10, dispersal_events = 20,
                         access_threshold = 5, simulation_period = 50,
                         stable_lgm = 7, transition_to_lgm = 100, lgm_to_current = 7,
                         stable_current = 13, scenario_span = 1,
                         barriers = NULL, barrier_format = NULL,
                         write_replicates = FALSE,
                         output_directory = "M_simulation", plot = TRUE,
                         mask_variables = FALSE, directory_masked = "Variables_M") {

  # --------
  # testing for initial requirements
  ## python 3.6 or superior
  cat("\nChecking dependencies:\n")
  py <- system("python", intern = TRUE)
  ncl <- gregexpr("Python 3.\\d", py)
  ncla <- regmatches(py, ncl)
  pyver <- unlist(ncla)

  if (length(pyver) == 0) {
    stop(paste("Python version >= 3.6 is needed to run the simulations. Please download and install it:\n",
               "https://www.python.org/downloads/release/python-366/", sep = ""))
  }

  version <- as.numeric(strsplit(pyver, " ")[[1]][2])

  if (version < 3.6) {
    stop(paste("Python version >= 3.6 is needed to run the simulations. Please download and install it:\n",
               "https://www.python.org/downloads/release/python-366/", sep = ""))
  }else{
    cat("Presence of Python 3.6 or superior confirmed.\n")
  }

  ## other arguments
  if (missing(current_variables)) {
    stop("Argument current_variables must be defined. See functions help.")
  }
  var <- list.files(current_variables, pattern = ".asc$", full.names = TRUE) # initial variables
  n <- length(var)

  if (n < 2) {
    stop("At least 2 variables are needed in current_variables. See functions help.")
  }
  if (project == TRUE) {
    if (missing(projection_variables)) {
      stop("If projections are needed, argument projection_variables must be defined. See functions help.")
    }
    varss <- raster::raster(var[1])

    lgmss <- raster::raster(list.files(projection_variables, pattern = ".asc$", full.names = TRUE)[1])
    if (raster::extent(varss) != raster::extent(lgmss)) {
      stop("Layers in projection_variables must have the same extent than layers in current_variables.")
    }
  }


  # --------
  # preparation of data in R
  cat("\nPreparing data for running simulation, please wait...\n")

  # output directory and slash type for directories
  dir.create(output_directory)
  sl <- "/"

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
    variables <- pca_raster(variables = current_variables, scale = scale,
                            project = project, return_back = TRUE,
                            n_pcs = npcs, output_directory = opca_fol)[[3]]

    ## suitability
    suit_mod <- ellipsoid_suitability(data, variables, suitability_threshold)
    suit_layer <- suit_mod[[5]]

    ### write suitability layer
    s_name <- (paste(suit_fol, "suitability.asc", sep = "/"))
    raster::writeRaster(suit_layer, s_name, format = "ascii")

    suit_name <- gsub("/", sl, paste("\"", paste(getwd(), s_name, sep = sl),
                                     "\"", sep = ""))

  }else {
    pvar <- list.files(projection_variables, pattern = ".asc$", full.names = TRUE)
    p_stack <- raster::stack(pvar)
    pcs <- pca_raster(variables = current_variables, scale = scale, project = project,
                      projection_variables = p_stack, return_back = TRUE,
                      n_pcs = npcs, output_directory = opca_fol)[c(3, 4)]

    variables <- pcs[[1]]
    lgm <- pcs[[2]]

    ## initial suitabilities
    suit_mod <- ellipsoid_suitability(data, variables, suitability_threshold,
                                      project = TRUE, projection_variables = lgm)

    suit_layer <- suit_mod[[5]]
    suit_lgm <- suit_mod[[6]]

    ### write suitability layer current and lgm
    s_name <- (paste(suit_fol, "suitability_current.asc", sep = "/"))
    raster::writeRaster(suit_layer, s_name, format = "ascii")

    sink(paste(suit_fol, "ellipsoid_metadata.txt", sep = "/"))
    print(suit_mod[c(2:4)])
    sink()

    l_name <- (paste(suit_fol, "suitability_lgm.asc", sep = "/"))
    raster::writeRaster(suit_lgm, l_name, format = "ascii")

    ### names for python
    sp_name <- paste0("\"", paste(getwd(), s_name, sep = sl), "\"")
    lp_name <- paste0("\"", paste(getwd(), l_name, sep = sl), "\"")

    ## preparing interpolation cycles
    int_vals <- interpolation_values(transition_to_lgm, stable_lgm, lgm_to_current,
                                     stable_current, simulation_period, scenario_span)
    types_clim <- int_vals[[1]]
    pos_scenarios <- int_vals[[2]]

    ## interpolations and sutiability layer projections
    suit_name <- interpolation(suit_mod, types_clim, pos_scenarios, variables,
                               lgm, sp_name, lp_name, sl, opca_fol, suit_fol)
  }

  # --------
  # occurences in suitable areas
  occ_suit <- suit_mod[[1]][, 1:2]
  suit_bar <- raster::extract(raster::raster(gsub("\"", "", suit_name[1])), occ_suit)
  occ_suit <- occ_suit[suit_bar > 0, ]

  ## records
  oca <- data.frame(Species = sp_nam, occ_suit)
  oca_nam <- paste(output_directory, "occ_simulation.csv", sep = "/")
  write.csv(oca, oca_nam, row.names = FALSE)

  occ_name <- gsub("/", sl, paste0("\"", paste(getwd(), oca_nam, sep = sl),"\""))

  # --------
  # figure of niche centroid model in E space
  save_nicheplot(suit_mod, suitability_threshold, variables,
                 size_proportion = 0.55, suit_fol, plot)

  # --------
  # other python needs
  ## folder for simulation's replicates
  if (write_replicates == TRUE) {
    rep_fol <- paste0(output_directory, "/Replicates")
    dir.create(rep_fol)
  }

  # --------
  # python script preparation and execution
  out_dir <- gsub("/", sl,
                  paste0("\"", paste(getwd(), output_directory, sep = sl), "\""))

  ## script
  dispersal_simulation(data = occ_name, suit_layers = suit_name,
                       dispersal_kernel = dispersal_kernel,
                       kernel_spread = kernel_spread,
                       max_dispersers = max_dispersers,
                       replicates = replicates, dispersal_events = dispersal_events,
                       access_threshold = access_threshold,
                       write_replicates = write_replicates, output_directory = out_dir)

  ## execution
  cat("\nRunning simulation, please wait...\n")

  if(.Platform$OS.type == "unix") {
    system("python -m pip install numpy", show.output.on.console = FALSE)
    py.script_r <- paste(getwd(), output_directory, "M_simulation", sep = "/")
    system(paste("python", paste(py.script_r, ".py", sep = "")))

  } else{
    system("python -m pip install numpy", show.output.on.console = FALSE)
    sil <- "\\\\"
    py.script_r <- paste(gsub("/", sil, paste0(getwd(), "/", output_directory)),
                         "M_simulation", sep = "\\")
    system(paste("python", paste(py.script_r, ".py", sep = "")))
  }

  # --------
  # preparing, writing, and outputs
  ## M
  cat("\nPreparing M in shapefile format...\n")
  m <- M_preparation(output_directory, pattern = "A_S\\d.*c$")
  m_poly <- m[[2]]
  m <- m[[1]]

  ## variables masked to M, if asked
  if (mask_variables == TRUE) {
    cat("\nMasking variables to M and writing them in", directory_masked,
        "    Please wait...\n")
    m_variables <- raster::mask(raster::crop(variab, m), m)

    var_names <- names(m_variables)
    var_names <- paste0(paste(directory_masked, var_names, sep = "/"), ".asc")
    dir.create(directory_masked)

    for (i in 1:length(raster::unstack(m_variables))) {
      raster::writeRaster(m_variables[[i]], filename = var_names[i], format = "ascii")
    }
  }

  ## plot and save figure of M in geographic space
  save_Mplot(suit_mod, variables, suit_layer, m_poly, size_proportion = 0.55,
             output_directory, plot)

  cat(paste("\nM simulation finished.\n",
            "Check your working directory:\t", getwd(), sep = ""))
}
