# writes ellipsoid metadata
write_ellmeta <- function(ellipsoid_model, file) {
  file <- paste0(file, ".txt")
  cat("Ellipsoid metadata:", file = file)
  cat("\n\nCentroid:\n", file = file, append = TRUE)
  suppressWarnings(write.table(t(as.matrix(ellipsoid_model[[2]])),
                               sep = "\t", file = file,
                               append = TRUE, quote = FALSE, row.names = FALSE))
  cat("\n\nCovariance matrix:\n", file = file, append = TRUE)
  cat(c(" ", colnames(ellipsoid_model[[3]]), "\n"), sep = "\t",
      file = file, append = TRUE)
  suppressWarnings(write.table(ellipsoid_model[[3]], sep = "\t", file = file,
                               append = TRUE, quote = FALSE, col.names = FALSE))
  cat("\n\nSuitability proportions:\n", file = file, append = TRUE)
  suppressWarnings(write.table(ellipsoid_model[[4]], sep = "\t", file = file,
                               append = TRUE, quote = FALSE, row.names = FALSE))

  save(ellipsoid_model, file = gsub("txt$", "RData", file))
}



# saves plot of ellipsoid niche
save_nicheplot <- function(ellipsoid_model, suitability_threshold, variables,
                           size_proportion = 0.55, output_directory) {
  el <- ellipse::ellipse(x = ellipsoid_model[[3]], centre = ellipsoid_model[[2]],
                         level = (100 - suitability_threshold) / 100)

  background <- terra::as.data.frame(variables[[1:2]])
  nr <- nrow(background)
  if (nr > 5000) {background <- background[sample(nr, size = 5000), ]}

  xlim <- range(range(el[, 1]), range(background[, 1]))
  ylim <- range(range(el[, 2]), range(background[, 2]))
  ylim <- c(ylim[1], ylim[2] + (diff(ylim) * 0.1))

  ## saving figure
  png(paste0(output_directory, "/niche_ellipsoid.png"), width = 80, height = 80,
      units = "mm", res = 600)
  par(mar = c(4.5, 4, 0.5, 0.5), cex = size_proportion)
  plot(background, xlim = xlim, ylim = ylim, pch = 1, col = "grey75",
       xlab = "PC 1", ylab = "PC 2", cex = 0.55)
  lines(el, col = "black", lwd = 1)
  points(ellipsoid_model[[1]][, 3:4], pch = 16, col = "black", cex = 0.65)
  legend("topright", legend = c("Background", "Occurrences", "Niche ellipsoid"),
         lty = c(NA, NA, 1), lwd = c(NA, NA, 1), pch = c(1, 16, NA),
         col = c("grey75", "black", "black"), bty = "n", horiz = TRUE,
         cex = 0.85)
  dev.off()
}


# metadata from base raster layer
layer_metadata <- function(layer) {
  return(list(
    NW_vertex = terra::ext(layer)[c(4, 1)],
    cell_size = terra::res(layer)[1],
    layer_dim = c(terra::nrow(layer), terra::ncol(layer))
  ))
}


# prepare S matrix from base raster
base_matrix <- function(layer) {
  return(terra::as.matrix(layer, wide = TRUE))
}


# binarize matrices according to a threshold percentage, threshold value, or is NA
binarize_matrix <- function(m, threshold_percentage = NULL,
                            threshold_value = NULL) {
  if (is.null(threshold_value) & is.null(threshold_percentage)) {
    m[!is.na(m)] <- 1
  } else {
    if (!is.null(threshold_percentage) & is.null(threshold_value)) {
      threshold_value <- quantile(m[m > 0], probs = threshold_percentage / 100)
    }
    m <- (m > threshold_value) * 1
  }

  return(m)
}


# matrix to raster layer
matrix_to_rlayer <- function(m, layer, name = NULL) {
  layer[] <- c(t(m))
  names(layer) <- ifelse(is.null(name), "mlayer", name)

  return(layer)
}


# calculates statistics among simulation replicates contained in a list
replicate_stats <- function(list_replicates, base_matrix, layer, threshold) {
  if (length(list_replicates) > 1) {
    list_replicates <- do.call(cbind, list_replicates)
    mean_A <- rowMeans(list_replicates)
    var_A <- rowSums((list_replicates - mean_A)^2) / (ncol(list_replicates) - 1)
  } else {
    mean_A <- unlist(list_replicates)
    var_A <- rep(0, length(mean_A))
  }

  all_meanA <- mean_A[mean_A > 0]
  threshold <- quantile(all_meanA, probs = threshold / 100)
  if (threshold == max(all_meanA)) {
    A_bin <- as.numeric(mean_A >= threshold)
  } else {
    A_bin <- as.numeric(mean_A > threshold)
  }

  layer_dim <- dim(base_matrix)

  mean_A <- matrix(mean_A, nrow = layer_dim[1], ncol = layer_dim[2])
  var_A <- matrix(var_A, nrow = layer_dim[1], ncol = layer_dim[2])
  A_bin <- matrix(A_bin, nrow = layer_dim[1], ncol = layer_dim[2])

  Sbin <- binarize_matrix(base_matrix)
  mean_A <- mean_A * Sbin
  var_A <- var_A * Sbin
  A_bin <- A_bin * Sbin

  mean_A <- matrix_to_rlayer(mean_A, layer, name = "mean")
  var_A <- matrix_to_rlayer(var_A, layer, name = "var")
  A_bin <- matrix_to_rlayer(A_bin, layer, name = "binary")

  return(list(mean = mean_A, var = var_A, bin = A_bin))
}


# finds raster format type according to format name
rformat_type <- function(format) {
  if (missing(format)) {stop("Argument 'format' needs to be defined")}
  if (format == "GTiff") {format1 <- ".tif"}
  if (format == "ascii") {format1 <- ".asc"}
  return(format1)
}


# write raster layers of simulation results
write_stats <- function(rep_stat_list, names, format, directory) {
  format1 <- rformat_type(format)
  terra::writeRaster(
    rep_stat_list[[1]], filename = paste0(directory, "/", names[1], format1)
  )
  terra::writeRaster(
    rep_stat_list[[2]], filename = paste0(directory, "/", names[2], format1)
  )
  terra::writeRaster(
    rep_stat_list[[3]], filename = paste0(directory, "/", names[3], format1)
  )
}



# transforms raster M in spatial polygon and saves it
M_preparation <- function(directory, A_bin = NULL, A_name = NULL,
                          raster_format = "GTiff") {
  if (is.null(A_bin) & !is.null(A_name)) {
    A_bin <- terra::rast(paste0(directory, "/", A_name))
  }

  A_bin[A_bin[] == 0] <- NA
  A_bin <- terra::trim(A_bin)
  shpm <- terra::as.polygons(A_bin, dissolve = TRUE)

  terra::writeVector(shpm,
                     filename = paste0(directory, "/accessible_area_M.shp"))
  m_name <- paste0(directory, "/accessible_area_M", rformat_type(raster_format))
  terra::writeRaster(A_bin, filename = m_name)

  return(list(A_bin, shpm))
}


# transparent colors
t_col <- function(col, alpha = 1, names = NULL) {
  rgb_col <- col2rgb(col)
  t_col <- rgb(rgb_col[1, ], rgb_col[2, ], rgb_col[3, ],
               alpha = (alpha * 100) * 255 / 100,
               names = names, maxColorValue = 255)

  return(t_col)
}



# saves plot of accessed areas (M) after simulation is done
save_Mplot <- function(ellipsoid_model, suitability_layer, M_polygon,
                       size_proportion = 0.55, output_directory) {

  png(paste0(output_directory, "/Accessible_area_M.png"), width = 80,
      height = 80, units = "mm", res = 600)
  par(cex = size_proportion)

  terra::plot(suitability_layer, legend = FALSE,
              mar = c(1.2, 1.2, 0.5, 0.5))
  terra::plot(M_polygon, lwd = 1, border = "blue4", add = TRUE)
  points(ellipsoid_model[[1]][, 1:2], pch = 16, cex = 0.6)

  dev.off()
}



# saves plot results from event - wise simulation
save_event_plot <- function(data, event_simulation_results, barriers = NULL,
                            size_proportion = 0.55, output_directory) {

  cola <- t_col("black", 0.5)
  mar1 <- c(1.2, 1.2, 1.5, 3.5)

  png(paste0(output_directory, "/Accessed_colonized_per_event.png"),
      width = 166, height = 80, units = "mm", res = 600)

  par(mfrow = c(1, 2))
  par(cex = size_proportion)
  terra::plot(event_simulation_results$A_events, main = "Accessed areas",
              legend = FALSE, mar = mar1)
  if (!is.null(barriers)) {
    terra::plot(barriers, col = cola, axes = FALSE, legend = FALSE, add = TRUE)
  }
  points(data[, 2:3], pch = 16, cex = 0.3)

  par(cex = size_proportion)
  terra::plot(event_simulation_results$C_events, main = "Colonized areas",
              mar = mar1)
  if (!is.null(barriers)) {
    terra::plot(barriers, col = cola, axes = FALSE, legend = FALSE, add = TRUE)
  }
  points(data[, 2:3], pch = 16, cex = 0.3)

  dev.off()
}
