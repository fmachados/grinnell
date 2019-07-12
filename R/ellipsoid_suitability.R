#' Suitability model based on niche centroid distance
#'
#' @description ellipsoid_suitability produces a model of environmental
#' suitability based on the Mahalanobis distance of each environmental value
#' in a given set of variables, to the centriod of the environmental values
#' characterized by the occurrences of a given species.
#'
#' @param data a numerical matrix containing geographic coordinates of the species
#' occurrences to be used; columns must be: longitude and latitude.
#' @param variables a RasterStack of variables representing the environmental
#' conditions in the area on which the models will be created.
#' @param suitability_threshold (numeric) value from 0 to 100 that will be used as
#' threshold (E); default = 5.
#' @param project (logical) whether or not to project the species niche to other
#' scenario(s). If TRUE, argument \code{projection_variables} needs to be defined.
#' Default = FALSE.
#' @param projection_variables a RasterStack (if only one scenario) or named list of
#' RasterStacks (if more than one scenario) of variables representing the
#' environmental conditions of the scenario at which the species niche will be
#' transferred. Variables names must correspond between initial and projection
#' scenarios.
#'
#' @return
#' A list containing: the occurrences (geographic coordinates and environmental
#' values), the corrdinates of the niche centroid, the covariance matrix, the
#' proportion of non-suitable areas as measured in the environmental (unique
#' environmental conditions) and the geographic space (all conditions), and a
#' suitability RasterLayer (if \code{projection_variables} is a RasterStack) or
#' a suitability RasterStack (if \code{projection_variables} is a list of
#' RasterStacks).
#'
#' @details Distance used for creating the model is Mahalanobis distance. All
#' values outside the ellipsoid produced by the distance measured from the
#' centroid limit of the occurrences exluding a percentage given by the
#' \code{suitability_threshold}, will have cero suitability values.
#'
#' Values in maps (from 0 to 1) are interpreted as suitability. For better
#' visualization use a color palette based on percetage of data with distinct
#' values.
#'
#' @export

ellipsoid_suitability <- function(data, variables, suitability_threshold = 5,
                                  project = FALSE, projection_variables,
                                  tolerance = 1e-60) {
  # conditions
  if (missing(data)) {
    stop("Argument data must be defined. See functions help.")
  }
  if (missing(variables)) {
    stop("Argument variables must be defined. See functions help.")
  }
  if (project == TRUE) {
    if (missing(projection_variables)) {
      stop("If projections are needed, argument projection_variables must be defined.\nSee functions help.")
    }

    if (class(projection_variables)[1] == "RasterStack" |
        class(projection_variables)[1] == "list") {
    }else {
      stop("If projections are needed, projection_variables must be a RasterStack or a list of RasterStacks.\nSee functions help.")
    }
  }

  # raster and data processing
  occ_data <- na.omit(raster::extract(variables, data))
  centroid <- apply(occ_data, 2, mean)

  covariance <- cov(occ_data)

  if (project == FALSE) {
    suit_model <- predict_esuitability(centroid = centroid,
                                       covariance_matrix = covariance,
                                       variables = variables,
                                       suitability_threshold = suitability_threshold,
                                       tolerance = tolerance)

    # occurrences for results
    occ_comp <- cbind(data[!is.na(raster::extract(variables[[1]], data)), ],
                      occ_data, na.omit(raster::extract(suit_model[[4]], data)))
    names(occ_comp) <- c("Longitude", "Latitude", names(variables), "Suitability")

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

    for (i in 1:length(projection_variables)) {
      suit_model <- predict_esuitability(centroid = centroid,
                                         covariance_matrix = covariance,
                                         variables = projection_variables[[i]],
                                         suitability_threshold = suitability_threshold,
                                         tolerance = tolerance)

      not_suitable[[i]] <- suit_model[[3]]

      # raster generation
      suit_layer[[i]] <- suit_model[[4]]
    }

    scenarios <- rep(c("Initial", paste("Transfer area",
                                        1:(length(projection_variables) - 1))),
                     each = 2)
    not_suitable <- data.frame(Scenario = scenarios, do.call(rbind, not_suitable))

    names(suit_layer) <- suit_names

    # occurrences for results
    occ_comp <- cbind(data[!is.na(raster::extract(variables[[1]], data)), ],
                      occ_data, na.omit(raster::extract(suit_layer[[1]], data)))
    names(occ_comp) <- c("Longitude", "Latitude", names(variables), "Suitability")

    # returning results
    results <- append(suit_layer, list(occurrences = occ_comp,
                                       centroid = suit_model[[1]],
                                       covariance_matrix = suit_model[[2]],
                                       non_suitable_area = not_suitable), after = 0)
  }

  return(results)
}


#' Predict suitability based on niche centroid distances
#'
#' @description predict_esuitability predicts a RasterLayer of environmental
#' suitability based on Mahalanobis distances to a niche centroid and a
#' covariance matrix.
#'
#' @param ellipsoid_model object produced with the function
#' \code{ellipsoid_suitability}. If defined, arguments \code{centroid} and
#' \code{covariance_matrix} are not necessary.
#' @param centroid centroid to wich distance will be measured in \code{variables}.
#' Length of centroid must correspond with number of layers in \code{variables}.
#' Ignored if \code{ellipsoid_model} is defined.
#' @param covariance_matrix covariance matrix to be used in measuring the
#' Mahalanobis distances. Ignored if \code{ellipsoid_model} is defined.
#' @param variables a RasterStack of variables representing the environmental
#' conditions in the area on which the models will be created.
#' @param suitability_threshold (numeric) value from 0 to 100 that will be used
#' as threshold (E); default = 5.
#' @param tolerance the tolerance for detecting linear dependencies.
#' Default = 1e-60.
#'
#' @export

predict_esuitability <- function(ellipsoid_model, centroid, covariance_matrix,
                                 variables, suitability_threshold = 5,
                                 tolerance = 1e-60) {

  # preparing data from ellipsoid_model if given
  if (!missing(ellipsoid_model)) {
    occ <- ellipsoid_model[[1]]
    centroid <- ellipsoid_model[[2]]
    covariance_matrix <- ellipsoid_model[[3]]
  }else {
    if (missing(centroid) | missing(covariance_matrix)) {
      stop("If argument ellipsoid_model is missing, centroid and covariance_matrix must be defined.")
    }
  }

  # raster data
  back <- na.omit(raster::values(variables))

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
  p_no_suit_g <- length(back[suitability != 0, 1]) / length(back[, 1])
  u_back <- unique(back)
  u_suit <- suitability[!duplicated(back)]
  p_no_suit_e <- length(unique(u_back[u_suit != 0, ])[, 1]) / length(u_back[, 1])

  not_suitable <- data.frame(c("Geographic_space", "Environmental_space"),
                             c(p_no_suit_g, p_no_suit_e))
  names(not_suitable) <- c("Space", "Proportion_not_suitable")

  # raster generation
  suit_layer <- variables[[1]]
  suit_layer[!is.na(raster::values(suit_layer))] <- suitability

  results <- list(centroid = centroid, covariance_matrix = covariance_matrix,
                  non_suitable_area = not_suitable, suitability_layer = suit_layer)

  return(results)
}
