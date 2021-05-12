#' Helper function to interpolate and prepare layer names
#'
#' @param ellipsoid_model object produced with the function
#' \code{ellipsoid_suitability}.
#' @param suitability_threshold (numeric) percentage of omission error to define
#' what is not suitable. Default = 5.
#' @param interpolation_values list of results from using the function
#' \code{\link{interpolation_values}}.
#' @param barriers RasterLayer representing dispersal barriers for the species.
#' This layer must have the same extent and projection than
#' \code{current_variables}. The only values allowed in this layer are 0 and 1;
#' 0 = barrier, 1 = all other areas. Default = NULL.
#' @param current_variables RasterStack of environmental variables representing
#' "current" conditions (interglacial).
#' @param projection_variables RasterStack of environmental variables
#' representing the "Last Glacial Maximum" scenario. Variable names, projection
#' and extent of these layers must be the same than those in
#' \code{current_variables}.
#' @param current_suitability (character) name of the file representing values
#' of suitability for the interglacial (current) scenario.
#' @param lgm_suitability (character) name of the file representing values of
#' suitability for the glacial (LGM) scenario.
#' @param out_format (character) format of layers to be written in
#' \code{output_directory}. Options are "ascii", "GTiff", and "EHdr" = bil.
#' Default = "GTiff".
#' @param output_directory (character) name of the folder where results will be
#' written.
#'
#' @export
#'
#' @usage
#' interpolation(ellipsoid_model, suitability_threshold = 5,
#'               interpolation_values, barriers = NULL,
#'               current_variables, projection_variables,
#'               current_suitability, lgm_suitability,
#'               out_format = "GTiff", output_directory)
#'
#' @return
#' A vector with names of files corresponding to each interpolation scenario.
#'
#' All layers resulting from interpolation will be written in
#' \code{output_directory}.

interpolation <- function(ellipsoid_model, suitability_threshold = 5,
                          interpolation_values, barriers = NULL,
                          current_variables, projection_variables,
                          current_suitability, lgm_suitability,
                          out_format = "GTiff", output_directory) {
  # initial tests
  if (missing(current_variables)) {
    stop("Argument 'current_variables' must be defined")
  }

  # interpolation values
  types_clim <- interpolation_values[[1]]
  pos_scenarios <- interpolation_values[[2]]

  # processing
  suit_name <- vector()
  for (i in 1:length(pos_scenarios)) {
    spot_val <- types_clim[pos_scenarios[i]]

    if (!spot_val %in% c(0, 1)) {
      ## interpolation of PCs
      pc_inter <- list()
      for (j in 1:3) {
        pc_inter[[j]] <- (current_variables[[j]] * (1 - spot_val)) +
          (projection_variables[[j]] * spot_val)
      }
      pc_inter <- do.call(raster::stack, pc_inter)

      ## suitability projections
      suit_p <- predict_esuitability(ellipsoid_model = ellipsoid_model,
                                     variables = pc_inter,
                                     suitability_threshold = suitability_threshold)

      ## barrier consideration
      if (!is.null(barriers)) {
        suit_p <- suit_p$suitability_layer * barriers
      }

      ## write suitability layer other scenarios
      ip_name <- paste0(output_directory, "/suitability_interpolation", i,
                        rformat_type(out_format))
      raster::writeRaster(suit_p, ip_name, format = out_format)

      suit_name[i] <- paste0("\"", normalizePath(ip_name), "\"")

      message("\tInterpolation ", i, " of ", length(pos_scenarios))
    } else {
      suit_name[i] <- ifelse(spot_val == 0, lgm_suitability,
                             current_suitability)

      message("\tInterpolation not needed, using glacial or interglacial layers")
    }
  }

  suit_name <- suit_name[rev(seq(1:length(suit_name)))]

  return(suit_name)
}



#' Helper function to generate interpolation values
#'
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
#' for changing the scenario. Default = 3 (three thousand years).
#'
#' @export
#'
#' @usage
#' interpolation_values(simulation_period = 140, transition_to_lgm = 13,
#'                      stable_lgm = 25, lgm_to_current = 7,
#'                      stable_current = 25, scenario_span = 3)
#'
#' @return
#' A list containing:
#' - a vector with proportions of interglacial conditions corresponding to
#' each interpolation
#' - a vector represenitng the position of each interpolation in the previous
#' vector.

interpolation_values <- function(simulation_period = 140, transition_to_lgm = 13,
                                 stable_lgm = 25, lgm_to_current = 7,
                                 stable_current = 25, scenario_span = 3) {

  # types of clime considering current
  types_clim <- rev(c(rev(seq(0, 1, 1 / (transition_to_lgm * 1000))),
                      rep(0, (stable_lgm * 1000) - 1),
                      seq(0, 1, 1 / (lgm_to_current * 1000)),
                      rep(1, (stable_current * 1000) - 1)))


  if ((simulation_period * 1000) - 1 > length(types_clim)) {
    timesbig <- ceiling(((simulation_period * 1000) -1) / length(types_clim))
    types_clim <- rep(types_clim, timesbig)
  }

  # position of each scenario in types of clime
  pos_scenarios <- seq(0, (simulation_period * 1000) - 1, scenario_span * 1000)
  pos_scenarios[1] <- 1

  return(list(climate_prop = types_clim, scenario_position = pos_scenarios))
}
