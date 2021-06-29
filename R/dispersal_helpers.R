# defines number of dispersers according to suitability and rules
nd_sval <- function(S_values, S_list, Nd_list) {
  sl <- c(0, S_list)
  nd <- sapply(1:(length(sl)-1), function(i) {
    ifelse(S_values > sl[i] & S_values <= sl[i+1], Nd_list[i], 0)
  })
  apply(nd, 1, sum)
}


# finds locations in the matrix according to coordinates and matrix features
set_loc <- function(coord, NW_vertex, cell_size) {
  loc_val <- cbind(row = NW_vertex[1] - coord[, 2],
                   column = coord[, 1] - NW_vertex[2])
  round(loc_val / cell_size)
}


# creates a matrix where all cells with records are set as 1 and the others as 0
set_pop <- function(spdata, NW_vertex, layer_dim, cell_size,
                    proportion = 0.5, set_seed = 1) {
  if (proportion == 1) {
    samp <- spdata[, 2:3]
  } else {
    set.seed(set_seed)
    samp <- spdata[sample(nrow(spdata), round(nrow(spdata) * proportion),
                          replace = FALSE), 2:3]
  }
  cells <- set_loc(samp, NW_vertex, cell_size)
  C <- matrix(data = 0, nrow = layer_dim[1], ncol = layer_dim[2])
  for (i in 1:nrow(samp)) {
    C[cells[i, 1], cells[i, 2]] <- 1
  }
  return(C)
}

# runs simulation steps to update accessed areas
dispersal_steps <- function(access_matrix, colonized_cells, n_dispersers,
                            dispersal_kernel = "normal", kernel_spread = 1,
                            set_seed = 1) {

  # initial tests
  if (missing(access_matrix)) {
    stop("Argument 'access_matrix' must be defined")
  }
  if (missing(colonized_cells)) {
    stop("Argument 'colonized_cells' must be defined")
  }
  if (missing(n_dispersers)) {
    stop("Argument 'n_dispersers' must be defined")
  }
  if (!dispersal_kernel %in% c("normal", "log_normal")) {
    stop("Argument 'dispersal_kernel' not valid")
  }

  # matrix dimensions
  layer_dim <- dim(access_matrix)

  # maximum number of dispersers
  Ndmax <- max(n_dispersers)

  # kernel logic
  fat <- ifelse(dispersal_kernel == "log_normal", TRUE, FALSE)

  # loops to run simulation
  for (l in 1:Ndmax) {
    rl <- n_dispersers >= l
    nv <- sum(rl)

    ## angle
    set.seed(set_seed)
    theta <- runif(n = nv, min = 0, max = 2) * pi

    ## rad
    if (fat == FALSE) {
      rad <- rnorm(n = nv, mean = 0, sd = kernel_spread)
    } else {
      rad <- rlnorm(n = nv, meanlog = 0, sdlog = kernel_spread)
    }

    ## new coordinates using old ones and rad
    d_lat <- round(rad * sin(theta))
    d_lon <- round(rad * cos(theta))

    cpos <- colonized_cells[rl, ]
    if (nv == 1) {
      cpos <- matrix(cpos, 1)
    }
    cond <- (cpos[, 1] + d_lat) >= 1 & (cpos[, 1] + d_lat) < layer_dim[1] &
      (cpos[, 2] + d_lon) >= 1 & (cpos[, 2] + d_lon) < layer_dim[2]

    rc_now <- cbind(
      row <- ifelse(cond, cpos[, 1] + d_lat, cpos[, 1]),
      col <- ifelse(cond, cpos[, 2] + d_lon, cpos[, 2])
    )

    access_matrix[rc_now] <- access_matrix[rc_now] + 1
  }

  return(access_matrix)
}
