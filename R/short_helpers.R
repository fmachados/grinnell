# define number of dispersers according to suitability and rules
nd_sval <- function(S_values, S_list, Nd_list) {
  sl <- c(0, S_list)
  nd <- sapply(1:(length(sl)-1), function(i) {
    ifelse(S_values > sl[i] & S_values <= sl[i+1], Nd_list[i], 0)
  })
  apply(nd, 1, sum)
}


# finds locations in the matrix according to coordinates and matrix features
setLoc <- function(coord, NW_vertex, cell_size) {
  loc_val <- cbind(row = NW_vertex[1] - coord[, 2],
                   column = coord[, 1] - NW_vertex[2])
  round(loc_val / cell_size)
}


# creates a matrix where all cells with records are set as 1 and the others as 0
setPop <- function(spdata, NW_vertex, layer_dim, cell_size,
                   proportion = 0.5, set_seed = 1) {
  if (proportion == 1) {
    samp <- spdata[, 2:3]
  } else {
    set.seed(set_seed)
    samp <- spdata[sample(nrow(spdata), round(nrow(spdata) * proportion),
                          replace = FALSE), 2:3]
  }
  cells <- setLoc(samp, NW_vertex, cell_size)
  C <- matrix(data = 0, nrow = layer_dim[1], ncol = layer_dim[2])
  for (i in 1:nrow(samp)) {
    C[cells[i, 1], cells[i, 2]] <- 1
  }
  return(C)
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
  if (missing(format)) {stop("Argument 'format' needs to be defined.")}
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
