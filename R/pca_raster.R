#' Principal componens for raster layers and projections
#'
#' @description pca_raster performs a principal component analysis with a set of
#' variables and produces raster layers of them. If needed the pricipal
#' components are projected to other scenarios.
#'
#' @param variables (character) name of the folder where environmental variables
#' are.
#' @param in_format (character) format of variables in \code{variables}. Options
#' are "ascii", "GTiff", and "EHdr" = bil. Default = "ascii".
#' @param out_format (character) format of variables to be written in distinct
#' sets inside \code{output_directory}. Options are "ascii", "GTiff", and
#' "EHdr" = bil. Default = "ascii".
#' @param scale (logical) wheter or not to scale variables while performing
#' principal component analyses.
#' @param project (logical) whether or not to project the species niche to other
#' scenario(s). If TRUE, argument \code{projection_variables} needs to be defined.
#' Default = FALSE.
#' @param projection_variables (character or RasterStack) if character, name of
#' the folder where subfolders with environmental variables of scenarios for
#' projections are (useful if multiple projections are needed). If RasterStack,
#' object containing stacked variables of only one projection scenario.
#' Variables must correspond with variables in \code{variables} (i.e., their name
#' must correspond but they should represent conditions in other scenario).
#' @param return_back (logical) whether or not return raster layers of principal
#' components to the R environment in a list with other results.
#' @param n_pcs (numeric) number of principal components to be returned as rasters.
#' By default all principal components are returned as RasterLayers.
#' @param output_directory (character) name of the folder to be created to save
#' the results of the analyses. Default = "PCA_results".
#'
#' @return
#' A list containing PCA summary and PCA loadings as matrices;
#' if \code{return_back} = TRUE, one or multiple (if projected)
#' RasterStacks of principal components are returned additionally.
#'
#' All results are written in \code{output_directory}.
#'
#' @details
#' If \code{scale} = TRUE, variables are centered to cero and scaled using
#' \code{\link[base]{scale}}.
#'
#' @export


pca_raster <- function(variables, in_format = "ascii", out_format = "ascii",
                       scale = TRUE, project = FALSE, projection_variables,
                       return_back = FALSE, n_pcs, output_directory = "PCA_results") {

  # checking for potential errors
  if (missing(variables)) {
    stop("Argument variables must be defined. See functions help.")
  }
  if (project == TRUE) {
    if (missing(projection_variables)) {
      stop("If projections are needed, argument projection_variables must be defined. See functions help.")
    }
  }

  # preparing raster formats
  if (in_format == "ascii") {
    patt <- ".asc$"
  }
  if (in_format == "GTiff") {
    patt <- ".tif$"
  }
  if (in_format == "EHdr") {
    patt <- ".bil$"
  }
  if (out_format == "ascii") {
    patt1 <- ".asc"
  }
  if (out_format == "GTiff") {
    patt1 <- ".tif"
  }
  if (out_format == "EHdr") {
    patt1 <- ".bil"
  }

  # reading variables
  var <- list.files(variables, pattern = patt, full.names = TRUE)
  variab <- raster::stack(var)
  var_points <- na.omit(raster::values(variab))

  # pca analyses
  if (scale == TRUE) {
    pca <- prcomp(var_points, center = TRUE, scale = TRUE)
  } else {
    pca <- prcomp(var_points, center = TRUE, scale = FALSE)
  }

  dir.create(output_directory)
  pca_fol <- paste(output_directory, "Initial", sep = "/")
  dir.create(pca_fol)

  if (missing(n_pcs)) {n_pcs <- length(var)}
  if (return_back == TRUE) {pcras <- list()}

  cat("\nWriting raster PCs in Output folder, please wait...\n")
  scores <- pca$x
  for (i in 1:n_pcs) {
    pcra <- variab[[1]]
    pcra[!is.na(raster::values(pcra))] <- scores[, i]
    filenam <- paste(pca_fol, "/pc_", i, patt1, sep = "")
    raster::writeRaster(pcra, filenam, format = out_format)

    if (return_back == TRUE) {pcras[[i]] <- pcra}
  }

  if (return_back == TRUE) {
    pcras <- do.call(raster::stack, pcras)
    names(pcras) <- paste0("pc_", 1:dim(pcras)[3])
  }

  StdDev <- pca$sdev
  VarExp <- pca$sdev^2/sum(pca$sdev^2)
  CumVar <- cumsum(VarExp)
  SumPCAMat <- rbind(StdDev, VarExp, CumVar)
  colnames(SumPCAMat) <- paste("PC", seq(1, length(StdDev)), sep = "")
  row.names(SumPCAMat) <- c("Standard deviation", "Proportion of Variance",
                            "Cumulative Proportion")

  sink(paste(paste(pca_fol, "pca_results.txt", sep = "/")))
  cat("Principal component analysis results\n")
  cat("\nPCA loadings\n")
  print(pca$rotation)

  cat("\n\nPCA summary\n")
  print(SumPCAMat)
  sink()

  # pca results to be returned
  loadings <- pca$rotation
  respca <- SumPCAMat

  # projecting PCs
  if (project == TRUE) {
    if (return_back == TRUE) {ppcrass <- list()}

    cat("\nProjecting and writing projected raster PCs in Output folder, please wait...\n")
    if (class(projection_variables)[1] == "character") {
      proj_dirs <- list.dirs(projection_variables, recursive = FALSE)
      proj_names <- list.dirs(projection_variables, recursive = FALSE,
                              full.names = FALSE)
      fol_names <- paste(output_directory, proj_names, sep = "/")
    }
    if (class(projection_variables)[1] %in% c("RasterStack", "RasterBrick")) {
      proj_dirs <- "projection"
      proj_names <- "Projected_PCs"
      fol_names <- paste(output_directory, proj_names, sep = "/")
    }

    for (h in 1:length(proj_dirs)) {
      if (class(projection_variables)[1] == "character") {
        pvar <- list.files(proj_dirs[h], pattern = patt, full.names = TRUE)
        p_stack <- raster::stack(pvar)
      }
      if (class(projection_variables)[1] %in% c("RasterStack", "RasterBrick")) {
        p_stack <- projection_variables
      }
      dir.create(fol_names[h])

      if (return_back == TRUE) {ppcras <- list()}
      p_stackp <- na.omit(raster::values(p_stack))
      names(p_stackp) <- names(pca[[4]])
      p_pcs <- predict(pca, newdata = p_stackp)

      for (i in 1:n_pcs) {
        pcra <- p_stack[[1]]
        pcra[!is.na(raster::values(pcra))] <- p_pcs[, i]
        filenam <- paste(fol_names[h], "/pc_", i, patt1, sep = "")
        raster::writeRaster(pcra, filenam, format = out_format)

        if (return_back == TRUE) {ppcras[[i]] <- pcra}
      }

      if (return_back == TRUE) {
        ppcrass[[h]] <- do.call(raster::stack, ppcras)
        names(ppcrass[[h]]) <- paste0("pc_", 1:dim(ppcrass[[h]])[3])
      }
    }

    if (return_back == TRUE) {
      names(ppcrass) <- paste("PCRasters", proj_names, sep = "_")
    }
  }

  if (return_back == TRUE) {
    if (project == TRUE) {
      results <- c(list(loadings, respca, pcras), ppcrass)
      names(results)[1:3] <- c("PCA_loadings", "PCA_results", "PCRasters_initial")
    }else {
      results <- list(loadings, respca, pcras)
      names(results) <- c("PCA_loadings", "PCA_results", "PCRasters_initial")
    }
  }else {
    results <- list(loadings, respca)
    names(results) <- c("PCA_loadings", "PCA_results")
  }

  cat("\nRaster PCA finished. Check your output directory",
      paste(getwd(), output_directory, sep = "/"), "\n")

  return(results)
}
