#' Suitability model based on Mahalanobis distance to niche centroid
#'
#' @description ellipsoid_suitability produces a model of environmental
#' suitability based on the environmental Mahalanobis distance to the centriod
#' of the environmental values characterized by species occurrences.
#'
#' @param data a numerical matrix containing geographic coordinates of species
#' occurrences to be used; columns must be: longitude and latitude.
#' @param variables a RasterStack of variables representing the environmental
#' conditions in the area of interest.
#' @param suitability_threshold (numeric) value (percentage) to be used as
#' threshold; default = 5.
#' @param project (logical) whether or not to project the model to other
#' scenario(s). If TRUE, argument \code{projection_variables} must be defined.
#' Default = FALSE.
#' @param projection_variables a RasterStack (if only one scenario) or named
#' list of RasterStacks (if more than one scenario) with variables representing
#' the environmental conditions of scenarios to transfer the model to.
#' Variable names must match between initial and projection scenarios.
#' @param tolerance the tolerance for detecting linear dependencies.
#' Default = 1e-60.
#'
#' @return
#' A list containing:
#' - occurrences (geographic coordinates and environmental values)
#' - niche centroid
#' - covariance matrix
#' - proportion of suitable areas in environmental and geographic space
#' - suitability RasterLayer if \code{projection_variables} is a RasterStack or
#' list of suitability RasterLayer(s) if \code{projection_variables} is a list
#'
#' @details Distance used for creating the model is Mahalanobis distance. All
#' values outside the ellipsoid produced using the centroid and covariance
#' matrix derived from environmental values in \code{data}, and the value in
#' \code{suitability_threshold} will be zero.
#'
#' Values in maps (from 0 to 1) can be interpreted as suitability.
#'
#' @export
#' @importFrom raster extract
#' @importFrom stats cov mahalanobis na.omit qchisq
#'
#' @usage
#' ellipsoid_suitability(data, variables, suitability_threshold = 5,
#'                       project = FALSE, projection_variables,
#'                       tolerance = 1e-60)

ellipsoid_suitability <- function(data, variables, suitability_threshold = 5,
                                  project = FALSE, projection_variables,
                                  tolerance = 1e-60) {
  # conditions
  if (missing(data)) {
    stop("Argument 'data' must be defined")
  }
  if (missing(variables)) {
    stop("Argument 'variables' must be defined")
  }
  if (project == TRUE) {
    if (missing(projection_variables)) {
      stop("If projections are needed, argument 'projection_variables' must be defined")
    }

    if (!class(projection_variables)[1] %in% c("RasterStack", "list")) {
      stop("Argument 'projection_variables' must be of class RasterStack or list")
    }
  }

  # raster and data processing
  occ_data <- na.omit(raster::extract(variables, data))
  centroid <- colMeans(occ_data)

  covariance <- cov(occ_data)

  if (project == FALSE) {
    suit_model <- predict_esuitability(centroid = centroid,
                                       covariance_matrix = covariance,
                                       variables = variables,
                                       suitability_threshold = suitability_threshold,
                                       tolerance = tolerance)

    # occurrences for results
    suit <- na.omit(raster::extract(suit_model[[4]], data))
    occ_comp <- cbind(occ_data, suit)
    colnames(occ_comp) <- c("Longitude", "Latitude", names(variables), "Suitability")

    results <- c(list(occurrences = occ_comp), suit_model)

  }else {
    # raster and data processing for projections
    not_suitable <- list()
    suit_layer <- list()

    if (class(projection_variables)[1] == "RasterStack") {
      projection_variables <- list(variables, projection_variables)
      proj_names <- "projection"

    }else {
      if (is.null(names(projection_variables))) {
        proj_names <- paste0("projection_", 1:length(projection_variables))
      } else {
        proj_names <- paste0("projection_", names(projection_variables))
      }

      projection_variables <- c(variables, projection_variables)
    }

    suit_names <- c("suitability_layer", proj_names)

    # running in loop
    for (i in 1:length(projection_variables)) {
      suit_model <- predict_esuitability(centroid = centroid,
                                         covariance_matrix = covariance,
                                         variables = projection_variables[[i]],
                                         suitability_threshold = suitability_threshold,
                                         tolerance = tolerance)

      not_suitable[[i]] <- suit_model[[3]]

      # raster
      suit_layer[[i]] <- suit_model[[4]]
    }

    scenarios <- rep(c("Initial", paste("Transfer area",
                                        1:(length(projection_variables) - 1))),
                     each = 2)
    not_suitable <- data.frame(Scenario = scenarios, do.call(rbind, not_suitable))

    names(suit_layer) <- suit_names

    # occurrences for results
    suit <- na.omit(raster::extract(suit_layer[[1]], data))
    occ_comp <- cbind(occ_data, suit)
    colnames(occ_comp) <- c("Longitude", "Latitude", names(variables), "Suitability")

    # returning results
    results <- list(occurrences = occ_comp, centroid = suit_model[[1]],
                    covariance_matrix = suit_model[[2]],
                    suitable_area_prop = not_suitable,
                    suitability_layer = suit_layer)
  }

  return(results)
}


#' Predict suitability based on Mahalanobis distance to niche centroid
#'
#' @description predict_esuitability predicts a RasterLayer of environmental
#' suitability based on Mahalanobis distances to a niche centroid and a
#' covariance matrix.
#'
#' @param ellipsoid_model object produced with the function
#' \code{ellipsoid_suitability}. If defined, arguments \code{centroid} and
#' \code{covariance_matrix} are not used. Default = NULL. This argument can be
#' overlooked if \code{centroid} and \code{covariance_matrix} are defined.
#' @param centroid centroid from which distances are measured in \code{variables}.
#' The length and variables of centroid must correspond with those of layers in
#' \code{variables}. Ignored if \code{ellipsoid_model} is defined. Default = NULL.
#' @param covariance_matrix covariance matrix to be used when measuring the
#' Mahalanobis distances. Ignored if \code{ellipsoid_model} is defined.
#' Default = NULL.
#' @param variables a RasterStack of variables representing the environmental
#' conditions in the area where the model will be transferred to.
#' @param suitability_threshold (numeric) value (percentage) to be used as
#' threshold; default = 5.
#' @param tolerance the tolerance for detecting linear dependencies.
#' Default = 1e-60.
#'
#' @export
#'
#' @usage
#' predict_esuitability(ellipsoid_model = NULL, centroid = NULL,
#'                      covariance_matrix = NULL, variables,
#'                      suitability_threshold = 5, tolerance = 1e-60)
#'
#' @return
#' A list containing:
#' - niche centroid
#' - covariance matrix
#' - proportion of suitable areas in environmental and geographic space
#' - suitability RasterLayer

predict_esuitability <- function(ellipsoid_model = NULL, centroid = NULL,
                                 covariance_matrix = NULL, variables,
                                 suitability_threshold = 5, tolerance = 1e-60) {

  # preparing data from ellipsoid_model if given
  if (!missing(ellipsoid_model)) {
    occ <- ellipsoid_model[[1]]
    centroid <- ellipsoid_model[[2]]
    covariance_matrix <- ellipsoid_model[[3]]
  }else {
    if (missing(centroid) | missing(covariance_matrix)) {
      stop("Argument 'ellipsoid_model' is missing, centroid and covariance_matrix must be defined")
    }
  }

  # raster data
  back <- na.omit(variables[])

  # calculate de Mahalanobis distance
  maha <- mahalanobis(x = back, center = centroid, cov = covariance_matrix,
                      tol = tolerance)

  # defining what is inside or outside the ellipsoid
  alpha <- (100 - suitability_threshold) / 100
  chi_sq <- qchisq(alpha, ncol(back))

  # distances to suitabilities
  suitability <- exp(-0.5 * maha)
  suitability <- ifelse(maha / chi_sq <= 1, suitability, 0) # inside only

  # counting not suitable areas (G and E)
  p_no_suit_g <- sum(suitability != 0) / length(suitability)
  u_suit <- suitability[!duplicated(back)]
  p_no_suit_e <- sum(u_suit != 0) / length(u_suit)

  suitable <- data.frame(c("Geographic_space", "Environmental_space"),
                             c(p_no_suit_g, p_no_suit_e))
  colnames(suitable) <- c("Space", "Proportion_suitable")

  # raster generation
  suit_layer <- variables[[1]]
  suit_layer[!is.na(suit_layer[])] <- suitability

  results <- list(centroid = centroid, covariance_matrix = covariance_matrix,
                  suitable_area_prop = suitable, suitability_layer = suit_layer)

  return(results)
}
