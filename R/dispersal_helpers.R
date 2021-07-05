# dispersers according to values from suitability layer
disperser_rules <- function(layer, max_dispersers) {
  maxsuit <- max(layer[], na.rm = TRUE)
  Nd_list <- seq(1, max_dispersers)
  incNd <- maxsuit / max_dispersers
  S_list <- seq(incNd, maxsuit, incNd)

  return(list(Nd_list = Nd_list, S_list = S_list))
}

# defines number of dispersers according to suitability and rules
nd_sval <- function(S_values, disperser_rules) {
  S_list <- disperser_rules$S_list
  Nd_list <- disperser_rules$Nd_list
  sl <- c(0, S_list)
  nd <- sapply(1:(length(sl)-1), function(i) {
    ifelse(S_values > sl[i] & S_values <= sl[i+1], Nd_list[i], 0)
  })

  return(apply(nd, 1, sum))
}


# finds locations in the matrix according to coordinates and matrix features
set_loc <- function(coord, NW_vertex, cell_size) {
  loc_val <- cbind(row = NW_vertex[1] - coord[, 2],
                   column = coord[, 1] - NW_vertex[2])

  return(round(loc_val / cell_size))
}


# creates a matrix where all cells with records are set as 1 and the others as 0
set_pop <- function(data, NW_vertex, layer_dim, cell_size,
                    proportion = 1, rule = "random", set_seed = 1) {
  if (proportion == 1) {
    samp <- data[, 2:3]
  } else {
    nr <- nrow(data)
    n <- round(nr * proportion)

    set.seed(set_seed)
    if (rule == "random") {
      samp <- data[sample(nr, n), 2:3]
    } else {
      samp <- data[sample(nr, n, prob = data[, 4]), 2:3]
    }
  }

  cells <- set_loc(samp, NW_vertex, cell_size)
  C <- matrix(data = 0, nrow = layer_dim[1], ncol = layer_dim[2])

  for (i in 1:nrow(samp)) {
    C[cells[i, 1], cells[i, 2]] <- 1
  }

  return(C)
}


# generates an initial matrix with colonized cells according to a layer and points
initial_colonized <- function(data, base_layer, proportion = 1, rule = "random",
                              set_seed = 1) {
  l_meta <- layer_metadata(base_layer)

  C <- set_pop(data, l_meta$NW_vertex, l_meta$layer_dim, l_meta$cell_size,
               proportion, rule, set_seed)

  return(C)
}



# sample cells from suitable areas
suitable_cells <- function(suit_layer, data = NULL) {
  if (is.null(data)) {
    noz <- which((suit_layer[] > 0))
    suit_bar <- suit_layer[noz]
    noz <- raster::xyFromCell(suit_layer, noz)
  } else {
    suit_bar <- raster::extract(suit_layer, data[, 2:3])
    tokeep <- suit_bar > 0 & !is.na(suit_bar)
    noz <- data[tokeep, 2:3]
    suit_bar <- suit_bar[tokeep]
  }
  sp <- ifelse(is.null(data), "Species", as.character(data[1, 1]))


  return(data.frame(species = sp, longitude = noz[, 1], latitude = noz[, 2],
                    suitability = suit_bar))
}


# id cells that have been colonized in a matrix
which_colonized <- function(colonized_matrix, suitability_matrix = NULL,
                            proportion = 1, rule = "random", set_seed = 1) {
  if (!is.null(suitability_matrix)) {
    colonized_cells <- which(colonized_matrix >= 1 & suitability_matrix > 0,
                             arr.ind = TRUE)
  } else {
    colonized_cells <- which(colonized_matrix >= 1, arr.ind = TRUE)
  }

  if (proportion < 1) {
    n <- nrow(colonized_cells)
    ns <- ceiling(n * proportion)

    set.seed(set_seed)
    colonized_cells <- colonized_cells[sample(n, ns), ]
    if (rule == "random") {
      colonized_cells <- colonized_cells[sample(n, ns), ]
    } else {
      suit <- suitability_matrix[colonized_cells]
      colonized_cells <- colonized_cells[sample(n, ns, prob = suit), ]
    }

    if (ns == 1) {
      colonized_cells <- matrix(colonized_cells, ncol = 2)
      colnames(colonized_cells) <- c("row", "col")
    }
  }

  return(colonized_cells)
}


# derive angle and distance of a movement based on a dispersal kernel
angle_distance <- function(n = 1, dispersal_kernel = "normal", kernel_spread = 1,
                           set_seed = 1) {
  set.seed(set_seed)
  theta <- runif(n = n, min = 0, max = 2) * pi

  if (dispersal_kernel == "normal") {
    rad <- rnorm(n = n, mean = 0, sd = kernel_spread)
  }
  if (dispersal_kernel == "log_normal") {
    rad <- rlnorm(n = n, meanlog = 0, sdlog = kernel_spread)
  }

  return(list(angle = theta, distance = rad))
}


# get accessed cells based on angle_distance and previously
which_accessed <- function(colonized_cells, angle_distance, layer_dim) {
  d_lat <- round(angle_distance$distance * sin(angle_distance$angle))
  d_lon <- round(angle_distance$distance * cos(angle_distance$angle))

  cond <- (colonized_cells[, 1] + d_lat) >= 1 &
    (colonized_cells[, 1] + d_lat) < layer_dim[1] &
    (colonized_cells[, 2] + d_lon) >= 1 &
    (colonized_cells[, 2] + d_lon) < layer_dim[2]

  rc_now <- cbind(
    row <- ifelse(cond, colonized_cells[, 1] + d_lat, colonized_cells[, 1]),
    col <- ifelse(cond, colonized_cells[, 2] + d_lon, colonized_cells[, 2])
  )

  return(rc_now)
}


# runs simulation steps to update accessed areas
dispersal_steps <- function(colonized_matrix, suitability_matrix, disperser_rules,
                            proportion_to_disperse = 1, sampling_rule = "random",
                            dispersal_kernel = "normal", kernel_spread = 1,
                            set_seed = 1) {

  # initial tests
  if (missing(colonized_matrix)) {
    stop("Argument 'colonized_matrix' must be defined")
  }
  if (missing(suitability_matrix)) {
    stop("Argument 'suitability_matrix' must be defined")
  }
  if (missing(disperser_rules)) {
    stop("Argument 'disperser_rules' must be defined")
  }
  if (!dispersal_kernel %in% c("normal", "log_normal")) {
    stop("Argument 'dispersal_kernel' not valid")
  }
  if (!sampling_rule %in% c("random", "suitability")) {
    stop("Argument 'sampling_rule' is not valid")
  }

  # data preparation
  colonized_cells <- which_colonized(colonized_matrix, suitability_matrix,
                                     proportion_to_disperse, sampling_rule,
                                     set_seed)
  Sv <- suitability_matrix[colonized_cells]
  n_dispersers <- nd_sval(Sv, disperser_rules)

  ## matrix dimensions
  layer_dim <- dim(colonized_matrix)

  ## access matrix
  access_matrix <- matrix(0, nrow = layer_dim[1], ncol = layer_dim[2])

  ## maximum number of dispersers
  Ndmax <- max(n_dispersers)

  # loop to run simulation
  for (l in 1:Ndmax) {
    rl <- n_dispersers >= l
    nv <- sum(rl)

    ## angle and distance of movement
    set_seed1 <- set_seed + l - 1
    ad <- angle_distance(nv, dispersal_kernel, kernel_spread, set_seed1)

    ## correcting C position according to n of dispersers
    cpos <- colonized_cells[rl, ]
    if (nv == 1) {
      cpos <- matrix(cpos, 1)
    }

    ## new coordinates using old ones, angles and directions
    rc_now <- which_accessed(cpos, ad, layer_dim)

    access_matrix[rc_now] <- access_matrix[rc_now] + 1
  }

  return(access_matrix)
}



# update A after dispersal step
update_accessed <- function(accessed_matrix, accessed_matrix_now) {
  whichA <- which(accessed_matrix_now > 0, arr.ind = TRUE)
  accessed_matrix[whichA] <- accessed_matrix[whichA] + accessed_matrix_now[whichA]

  return(accessed_matrix)
}


# update C after dispersal step
update_colonized <- function(colonized_matrix, accessed_matrix,
                             suitability_matrix) {
  whichC <- which(accessed_matrix > 0 & suitability_matrix > 0, arr.ind = TRUE)
  colonized_matrix[whichC] <- colonized_matrix[whichC] + accessed_matrix[whichC]

  return(colonized_matrix)
}
