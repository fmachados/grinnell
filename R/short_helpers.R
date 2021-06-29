# writes ellipsoid metadata
write_ellmeta <- function(emodel, file) {
  file <- paste0(file, ".txt")
  cat("Ellipsoid metadata:", file = file)
  cat("\n\nCentroid:\n", file = file, append = TRUE)
  suppressWarnings(write.table(t(as.matrix(emodel[[2]])), sep = "\t", file = file,
                               append = TRUE, quote = FALSE, row.names = FALSE))
  cat("\n\nCovariance matrix:\n", file = file, append = TRUE)
  cat(c(" ", colnames(emodel[[3]]), "\n"), sep = "\t", file = file, append = TRUE)
  suppressWarnings(write.table(emodel[[3]], sep = "\t", file = file,
                               append = TRUE, quote = FALSE, col.names = FALSE))
  cat("\n\nSuitability proportions:\n", file = file, append = TRUE)
  suppressWarnings(write.table(emodel[[4]], sep = "\t", file = file,
                               append = TRUE, quote = FALSE, row.names = FALSE))

  save(emodel, file = gsub("txt$", "RData", file))
}



# saves plot of ellipsoid niche
save_nicheplot <- function(ellipsoid_model, suitability_threshold, variables,
                           size_proportion = 0.55, output_directory) {
  el <- ellipse::ellipse(x = ellipsoid_model[[3]], centre = ellipsoid_model[[2]],
                         level = (100 - suitability_threshold) / 100)

  background <- na.omit(raster::values(variables[[1:2]]))
  nr <- nrow(background)
  if (nr > 5000) {background <- background[sample(nr, size = 5000), ]}

  xlim <- range(range(el[, 1]), range(background[, 1]))
  ylim <- range(range(el[, 2]), range(background[, 2]))

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
         col = c("grey75", "black", "black"), bty = "n", horiz = TRUE, cex = 0.85)
  dev.off()
}



# calculates statistics among simulation replicates contained in a list
stats_rep <- function(list_rep, layer_dim, Sbin, layer, threshold) {
  list_rep <- do.call(cbind, list_rep)
  mean_A <- rowMeans(list_rep)
  var_A <- rowSums((list_rep - mean_A)^2) / (ncol(list_rep) - 1)

  all_meanA <- mean_A[mean_A > 0]
  threshold <- quantile(all_meanA, probs = threshold / 100)
  A_bin <- as.numeric(mean_A > threshold)

  mean_A <- matrix(mean_A, nrow = layer_dim[1], ncol = layer_dim[2])
  var_A <- matrix(var_A, nrow = layer_dim[1], ncol = layer_dim[2])
  A_bin <- matrix(A_bin, nrow = layer_dim[1], ncol = layer_dim[2])

  mean_A <- mean_A * Sbin
  var_A <- var_A * Sbin
  A_bin <- A_bin * Sbin

  layer[] <- c(t(mean_A)); mean_A <- layer
  layer[] <- c(t(var_A)); var_A <- layer
  layer[] <- c(t(A_bin)); A_bin <- layer

  return(list(mean_A, var_A, A_bin))
}


# finds raster format type according to format name
rformat_type <- function(format) {
  if (missing(format)) {stop("Argument 'format' needs to be defined")}
  if (format == "GTiff") {format1 <- ".tif"}
  if (format == "EHdr") {format1 <- ".bil"}
  if (format == "ascii") {format1 <- ".asc"}
  return(format1)
}


# write raster layers of simulation results
write_stats <- function(rep_stat_list, names, format, directory) {
  format1 <- rformat_type(format)
  raster::writeRaster(
    rep_stat_list[[1]], filename = paste0(directory, "/", names[1], format1),
    format = format
  )
  raster::writeRaster(
    rep_stat_list[[2]], filename = paste0(directory, "/", names[2], format1),
    format = format
  )
  raster::writeRaster(
    rep_stat_list[[3]], filename = paste0(directory, "/", names[3], format1),
    format = format
  )
}



# transforms raster M in spatial polygon and saves it
M_preparation <- function(directory, A_bin = NULL, A_name = NULL,
                          raster_format = "GTiff") {
  if (is.null(A_bin) & !is.null(A_name)) {
    A_bin <- raster::raster(paste0(directory, "/", A_name))
  }

  A_bin[A_bin[] == 0] <- NA
  A_bin <- raster::trim(A_bin)
  shpm <- raster::rasterToPolygons(A_bin, dissolve = TRUE)
  sp::proj4string(shpm) <- A_bin@crs

  rgdal::writeOGR(shpm, directory, "accessible_area_M", driver = "ESRI Shapefile")
  m_name <- paste0(directory, "/accessible_area_M", rformat_type(raster_format))
  raster::writeRaster(A_bin, filename = m_name, format = raster_format)

  return(list(A_bin, shpm))
}



# saves plot of accessed areas (M) after simulation is done
save_Mplot <- function(ellipsoid_model, suitability_layer, M_polygon,
                       size_proportion = 0.55, output_directory) {
  boxpam <- t(matrix(suitability_layer@extent, 2, byrow = T))
  boxpam <- sp::SpatialPointsDataFrame(boxpam, data.frame(boxpam),
                                       proj4string = M_polygon@proj4string)

  png(paste0(output_directory, "/Accessible_area_M.png"), width = 80, height = 80,
      units = "mm", res = 600)
  par(mar = c(0.5, 0.5, 0.5, 0.5), cex = size_proportion)
  sp::plot(boxpam, col = NA)
  raster::image(suitability_layer, col = rev(terrain.colors(255)),
                axes = FALSE, add = TRUE)
  sp::plot(M_polygon, lwd = 1, border = "blue4", add = TRUE)
  points(ellipsoid_model[[1]][, 1:2], pch = 16, cex = 0.6)
  box()
  dev.off()
}
