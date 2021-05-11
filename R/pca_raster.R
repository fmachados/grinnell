#' Principal component analysis for raster layers
#'
#' @description pca_raster performs a principal component analysis with a set of
#' variables and returns results including PCs as raster layers. If needed,
#' principal components can be projected to other scenarios.
#'
#' @param variables RasterStack or character. If RasterStack, stack of layers to
#' be used; if character, name of the folder where environmental variables
#' are.
#' @param in_format (character) if class of \code{variables} = "character",
#' format of variables in such a folder. Options are "ascii", "GTiff", and
#' "EHdr" = bil. Default = NULL.
#' @param scale (logical) whether or not to scale variables while performing
#' principal component analyses.
#' @param center (logical) whether or not to center variables while performing
#' principal component analyses.
#' @param n_pcs (numeric) number of principal components to be returned as
#' raster layers. By default all principal components are returned.
#' @param project (logical) whether or not to project the species niche to other
#' scenario(s). If TRUE, argument \code{projection_variables} must be defined.
#' Default = FALSE.
#' @param projection_variables character or RasterStack. If character, name of
#' the folder where one or more sub-directories containing variables of distinct
#' environmental scenarios are (useful if multiple projections are needed).
#' If RasterStack, stack of variables representing one scenario for projection.
#' Variable names must correspond with those of \code{variables} (i.e., their
#' names must match).
#' @param return_projection (logical) whether to return raster PCs for
#' projection scenario(s) as part of the resulting list. Default = FALSE.
#' @param write_to_directory (logical) whether to write results in
#' \code{output_directory}. Default = FALSE.
#' @param out_format (character) format of layers to be written in
#' \code{output_directory}. Options are "ascii", "GTiff", and "EHdr" = bil.
#' Default = "GTiff".
#' @param output_directory (character) name of the folder to be created to save
#' results.
#'
#' @return
#' A list containing:
#' - PCA results
#' - RasterStack of principal components for \code{variables}
#' - if return_projection = TRUE, list of RasterStack objects with principal
#' components for other scenarios. NULL if \code{project} = FALSE.
#'
#' @export
#' @importFrom raster nlayers
#'
#' @usage
#' pca_raster(variables, in_format = NULL, scale = TRUE, center = TRUE,
#'            n_pcs = NULL, project = FALSE, projection_variables,
#'            return_projection = FALSE, write_to_directory = FALSE,
#'            out_format = "GTiff", output_directory)

pca_raster <- function(variables, in_format = NULL, scale = TRUE, center = TRUE,
                       n_pcs = NULL, project = FALSE, projection_variables,
                       return_projection = FALSE, write_to_directory = FALSE,
                       out_format = "GTiff", output_directory) {

  # checking for potential errors
  if (missing(variables)) {
    stop("Argument 'variables' must be defined")
  }
  clsv <- class(variables)[1]
  if (!clsv %in% c("character", "RasterStack", "RasterBrick")) {
    stop("'variables' must be of class 'character' or 'RasterStack'")
  }
  if (clsv == "character") {
    if (is.null(in_format)) {
      stop("If class of 'variables' = 'character', 'in_format' must be defined")
    }
  }
  if (project == TRUE) {
    if (missing(projection_variables)) {
      stop("If 'project' = TRUE, argument 'projection_variables' must be defined")
    }
    if (return_projection == FALSE & write_to_directory == FALSE) {
      message("Setting 'return_projection' = TRUE to keep results from projections")
      return_projection <- TRUE
    }
  }
  if (write_to_directory) {
    stop("If 'write_to_directory' = TRUE, 'output_directory' must be defined")
  }

  # preparing things
  patt1 <- rformat_type(out_format)

  if (write_to_directory == TRUE) {
    dir.create(output_directory)
    pca_fol <- paste0(output_directory, "/Initial")
    dir.create(pca_fol)
  }

  # reading variables
  if (clsv == "character") {
    var <- list.files(variables, pattern = patt, full.names = TRUE)
    variables <- raster::stack(var)
  }
  var_points <- na.omit(variables[])
  pcra <- variables[[1]]
  to_fill <- !is.na(pcra[])

  if (is.null(n_pcs)) {n_pcs <- ncol(var_points)}

  # pca analyses
  pca <- prcomp(var_points, center = center, scale = scale)
  pcras <- list()

  message("Preparing raster PCs...")
  for (i in 1:n_pcs) {
    pcra[to_fill] <- pca$x[, i]
    pcras[[i]] <- pcra

    if (write_to_directory == TRUE) {
      filenam <- paste0(pca_fol, "/PC", i, patt1)
      raster::writeRaster(pcra, filenam, format = out_format)
    }
  }

  pcras <- do.call(raster::stack, pcras)
  names(pcras) <- paste0("PC", 1:(raster::nlayers(pcras)))

  SumPCAMat <- summary(pca)$importance

  if (write_to_directory == TRUE) {
    txtfile <- paste0(pca_fol, "/PCA_summary.txt")
    cat("Principal component analysis:\n\nPCA summary\n", file = txtfile)
    suppressWarnings(write.table(SumPCAMat, sep = "\t", file = txtfile,
                                 append = TRUE, quote = FALSE))

    cat("\n\nPCA loadings\n", file = txtfile, append = TRUE)
    suppressWarnings(write.table(pca$rotation, sep = "\t", file = txtfile,
                                 append = TRUE, quote = FALSE))

    rd_file <- paste0(pca_fol, "/PCA_results.RData")
    save(pca, file = rd_file)
  }

  # projecting PCs
  if (project == TRUE) {
    clspr <- class(projection_variables)[1]
    if (return_projection == TRUE) {ppcrass <- list()} ########

    message("Projecting PCs to distinct scenarios...")
    if (clspr == "character") {
      proj_dirs <- list.dirs(projection_variables, recursive = FALSE)
      proj_names <- list.dirs(projection_variables, recursive = FALSE,
                              full.names = FALSE)
      fol_names <- paste(output_directory, proj_names, sep = "/")
    } else {
      proj_dirs <- "projection"
      proj_names <- "Projected_PCs"
      fol_names <- paste(output_directory, proj_names, sep = "/")
    }

    for (h in 1:length(proj_dirs)) {
      if (clspr == "character") {
        pvar <- list.files(proj_dirs[h], pattern = patt, full.names = TRUE)
        projection_variables <- raster::stack(pvar)
      }
      pcra <- projection_variables[[1]]
      to_fillp <- !is.na(pcra[])

      if (write_to_directory == TRUE) {
        dir.create(fol_names[h])
      }

      if (return_projection == TRUE) {ppcras <- list()}

      p_stackp <- na.omit(projection_variables[])
      p_pcs <- predict(pca, newdata = p_stackp)

      for (i in 1:n_pcs) {
        pcra[to_fillp] <- p_pcs[, i]
        if (return_projection == TRUE) {ppcras[[i]] <- pcra}

        if (write_to_directory == TRUE) {
          filenam <- paste0(fol_names[h], "/PC", i, patt1)
          raster::writeRaster(pcra, filenam, format = out_format)
        }
      }

      if (return_projection == TRUE) {
        ppcrass[[h]] <- do.call(raster::stack, ppcras)
        names(ppcrass[[h]]) <- paste0("PC", 1:(raster::nlayers(ppcrass[[h]])))
      }
    }

    if (return_projection == TRUE) {
      names(ppcrass) <- paste0("PCRaster_", proj_names)
    }
  }

  if (project == TRUE & return_projection == TRUE) {
    results <- list(PCA_results = pca, PCRaster_initial = pcras,
                    PCRaster_projection = ppcrass)
  } else {
    results <- list(PCA_results = pca, PCRaster_initial = pcras,
                    PCRaster_projection = NULL)
  }

  message("Raster PCA finished")
  if (write_to_directory == TRUE) {
    message("Check results in:  ", normalizePath(output_directory))
  }

  return(results)
}
