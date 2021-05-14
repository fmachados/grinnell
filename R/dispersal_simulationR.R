#' Helper function to simulate species dispersal processes
#'
#' @description dispersal_simulationR performs a multistep simulation of species
#' dispersal to reconstruct areas that have been accessed and/or colonized
#' based on environmental suitability and user-defined dispersal parameters.
#'
#' @param data data.frame containing geographic coordinates of occurrences of
#' the species of interest. Columns must be: species, longitude, latitude, in
#' that order.
#' @param suit_layers (character) vector of names of suitability layers to be
#' used as distinct scenarios. If more than one, the layer names must be ordered
#' starting from the oldest scenario. Layer names should include parent
#' directory if needed.
#' @param starting_porportion (numeric) proportion of \code{data} to be used as
#' starting points for the simulation. Default = 0.5. All data is used if a
#' value of 1 is defined.
#' @param dispersal_kernel (character) dispersal kernel (dispersal function)
#' used to simulate the movement of the species. Options are: "normal",
#' "log_normal". Default = "normal".
#' @param kernel_spread (numeric) standard deviation for the
#' \code{dispersal_kernel}. Default = 1.
#' @param max_dispersers (numeric) maximum number of dispersers that depart from
#' each colonized pixel. Depending on suitability this number will automatically
#' decrease in areas with low suitability values. Default = 4.
#' @param replicates (numeric) number of times to repeat the simulation
#' per scenario. Default = 10.
#' @param dispersal_events (numeric) number of dispersal events to happen per
#' scenario. Default = 25.
#' @param threshold (numeric) percentage to be considered when excluding
#' accessed or colonized cells with lower values. Default = 5.
#' @param set_seed (numeric) a seed to be used when sampling \code{data}
#' according to \code{starting_porportion}. Default = 1.
#' @param return (character) the results to be returned or written. Options are:
#' "all", "accessed", "colonized". Default = "all"
#' @param write_to_directory (logical) whether to write results in
#' \code{output_directory}. Default = FALSE.
#' @param write_all_scenarios (logical) whether or not to write results
#' for all scenarios. The default, FALSE, writes only the final results if
#' \code{write_to_directory} = TRUE.
#' @param raster_format (character) format to use for raster layers to be
#' written. Options are: "GTiff", "EHdr", and "ascii", Default = "GTiff.
#' @param output_directory (character) name of the output directory where
#' results should be written. If this directory does not exist, it will be
#' created.
#'
#' @export
#' @importFrom raster raster extent res nrow ncol writeRaster
#' @importFrom stats quantile runif rnorm rlnorm
#' @importFrom utils read.csv write.csv write.table
#'
#' @usage
#' dispersal_simulationR(data, suit_layers, starting_porportion = 0.5,
#'                       dispersal_kernel = "normal",
#'                       kernel_spread = 1, max_dispersers = 4,
#'                       replicates = 10, dispersal_events = 25,
#'                       threshold = 5, set_seed = 1,
#'                       return = "all", write_to_directory = FALSE,
#'                       write_all_scenarios = FALSE,
#'                       raster_format = "GTiff", output_directory)
#'
#' @return
#' If \code{return} = "all', all elements described below will be returned as a
#' list, if "accessed" or "colonized" are chosen instead, only the elements
#' corresponding to either "accessed" or "colonized" areas will be returned.
#'
#' The list returned contains:
#' - a list with a summary of scenarios and parameters used
#' - a binary RasterLayer showing accessed = 1 and non-accessed areas = 0
#' - a RasterLayer representing mean values of accessibility frequency among
#' replicates
#' - a RasterLayer representing variance among values of accessibility frequency
#' of all replicates
#' - a binary RasterLayer showing colonized = 1 and non-colonized areas = 0
#' - a RasterLayer representing mean values of frequency of colonization among
#' replicates
#' - a RasterLayer representing variance among values of frequency of
#' colonization of all replicates
#'
#' If \code{write_to_directory} is set to TRUE, the elements described above
#' and raster layers corresponding to all scenarios
#' (if \code{write_all_scenarios} = TRUE) will be written in
#' \code{output_directory}.

dispersal_simulationR <- function(data, suit_layers, starting_porportion = 0.5,
                                  dispersal_kernel = "normal",
                                  kernel_spread = 1, max_dispersers = 4,
                                  replicates = 10, dispersal_events = 25,
                                  threshold = 5, set_seed = 1,
                                  return = "all", write_to_directory = FALSE,
                                  write_all_scenarios = FALSE,
                                  raster_format = "GTiff", output_directory) {

  # initial tests
  if (missing(data)) {
    stop("Argument 'data' must be defined")
  }
  if (missing(suit_layers)) {
    stop("Argument 'suit_layers' must be defined")
  }
  if (write_to_directory & missing(output_directory)) {
    stop("If 'write_to_directory' = TRUE, 'output_directory' must be defined")
  }


  # data
  list_suit <- suit_layers
  spdata <- data

  # parameters
  proportion <- starting_porportion
  rep <- replicates
  steps <- dispersal_events
  spread <- kernel_spread
  NdMax <- max_dispersers
  fat <- ifelse(dispersal_kernel == "LogNormal", TRUE, FALSE)

  # initial part of report
  if (write_to_directory == TRUE) {
    if (dir.exists(output_directory) == FALSE) {
      dir.create(output_directory)
    }
    outText <- paste0(output_directory, "/report.txt")
    if (file.exists(outText)) {invisible(file.remove(outText))}
    cat(
      "Simulation parameters\n\n",
      "   Suitability scenarios: ", length(list_suit), "\n",
      "   Replicates: ", rep, "\n",
      "   Dispersal events: ", steps, "\n",
      "   Dispersal kernel: ", dispersal_kernel, "\n",
      "   Kernel spread (SD): ", spread, "\n",
      "   Maximum number of dispersers: ", NdMax, "\n",
      file = outText,
      append = TRUE
    )
  }

  # initial values
  cur_layer <- raster::raster(list_suit[length(list_suit)])
  NW_vertex <- raster::extent(cur_layer)[c(4, 1)]
  cell_size <- raster::res(cur_layer)[1]
  layer_dim <- c(raster::nrow(cur_layer), raster::ncol(cur_layer))
  maxsuit <- max(cur_layer[], na.rm = TRUE)

  # other parameters
  Nd_list <- seq(1, NdMax)
  incNd <- maxsuit / NdMax
  S_list <- seq(incNd, maxsuit, incNd)

  # running
  ## start time
  start <- Sys.time()

  ## running in loop all analysis
  message("Running simulation:")
  for (i in 1:length(list_suit)) {
    s <- raster::raster(list_suit[i])
    S <- matrix(s[], nrow = layer_dim[1], ncol = layer_dim[2], byrow = TRUE)
    Sbin <- S
    Sbin[!is.na(Sbin)] <- 1

    ## lists depending on what is needed
    if(return == "all") {
      list_rep <- list(); list_col <- list()
    } else {
      if(return == "accessed") {
        list_rep <- list()
      } else {
        list_col <- list()
      }
    }

    ## loop for all replicates
    message("   Replicates:", appendLF = FALSE)

    for (j in 1:rep) {
      ## preparing C matrix
      if (i == 1) {
        C <- setPop(spdata, NW_vertex, layer_dim, cell_size, proportion,
                    set_seed)
        A <- C
      } else {
        C[S == 0] <- 0
      }

      ### simulation steps
      for (k in 1:steps) {
        A_now <- matrix(0, nrow = layer_dim[1], ncol = layer_dim[2], byrow = TRUE)
        Cpos <- which(C >= 1 & S > 0, arr.ind = TRUE)
        Sv <- S[Cpos]
        Nd <- nd_sval(Sv, S_list, Nd_list)
        Ndmax <- max(Nd)

        #### running steps according to dispersers
        for (l in 1:Ndmax) {
          rl <- Nd >= l
          nv <- sum(rl)

          ##### angle
          theta <- runif(n = nv, min = 0, max = 2) * pi

          ##### rad
          if (fat == FALSE) {
            rad <- rnorm(n = nv, mean = 0, sd = spread)
          } else {
            rad <- rlnorm(n = nv, meanlog = 0, sdlog = spread)
          }

          ##### new coordinates using old ones and rad
          d_lat <- round(rad * sin(theta))
          d_lon <- round(rad * cos(theta))

          cpos <- Cpos[rl, ]
          if (nv == 1) {
            cpos <- matrix(cpos, 1)
          }
          cond <- (cpos[, 1] + d_lat) >= 1 & (cpos[, 1] + d_lat) < layer_dim[1] &
            (cpos[, 2] + d_lon) >= 1 & (cpos[, 2] + d_lon) < layer_dim[2]

          rc_now <- cbind(
            row <- ifelse(cond, cpos[, 1] + d_lat, cpos[, 1]),
            col <- ifelse(cond, cpos[, 2] + d_lon, cpos[, 2])
          )

          A_now[rc_now] <- A_now[rc_now] + 1
        }

        #### updating A and C
        whichA <- which(A_now > 0, arr.ind = TRUE)
        A[whichA] <- A[whichA] + A_now[whichA]

        whichC <- which(A_now > 0 & S > 0, arr.ind = TRUE)
        C[whichC] <- C[whichC] + A_now[whichC]
      }

      ### keeping replicates of accessed areas
      if(return == "all") {
        list_rep[[j]] <- c(A)
        list_col[[j]] <- c(C)
      } else {
        if(return == "accessed") {
          list_rep[[j]] <- c(A)
        } else {
          list_col[[j]] <- c(C)
        }
      }

      message(" ", j, appendLF = FALSE)
    }

    ## statistics
    if(return == "all") {
      Amvb <- stats_rep(list_rep, layer_dim, Sbin, s, threshold)
      Cmvb <- stats_rep(list_col, layer_dim, Sbin, s, threshold)
    } else {
      if(return == "accessed") {
        Amvb <- stats_rep(list_rep, layer_dim, Sbin, s, threshold)
      } else {
        Cmvb <- stats_rep(list_col, layer_dim, Sbin, s, threshold)
      }
    }

    ## writing results
    if (write_to_directory == TRUE) {
      if (write_all_scenarios == TRUE) {
        if(return == "all") {
          namesA <- paste0("A_", c("mean", "var", "bin"), "_", i)
          write_stats(Amvb, namesA, raster_format, output_directory)
          namesC <- paste0("C_", c("mean", "var", "bin"), "_", i)
          write_stats(Cmvb, namesC, raster_format, output_directory)
        } else {
          if(return == "accessed") {
            namesA <- paste0("A_", c("mean", "var", "bin"), "_", i)
            write_stats(Amvb, namesA, raster_format, output_directory)
          } else {
            namesC <- paste0("C_", c("mean", "var", "bin"), "_", i)
            write_stats(Cmvb, namesC, raster_format, output_directory)
          }
        }
      } else {
        if(i == length(list_suit)) {
          if(return == "all") {
            namesA <- paste0("A_", c("mean", "var", "bin"))
            write_stats(Amvb, namesA, raster_format, output_directory)
            namesC <- paste0("C_", c("mean", "var", "bin"))
            write_stats(Cmvb, namesC, raster_format, output_directory)
          } else {
            if(return == "accessed") {
              namesA <- paste0("A_", c("mean", "var", "bin"))
              write_stats(Amvb, namesA, raster_format, output_directory)
            } else {
              namesC <- paste0("C_", c("mean", "var", "bin"))
              write_stats(Cmvb, namesC, raster_format, output_directory)
            }
          }
        }
      }
    }

    message(""); message("Scenario ", i, " of ", length(list_suit))
  }

  ## preparing layers if needed

  ## end time
  end <- Sys.time()

  # last parts of report and preparing layers if needed
  if (write_to_directory == TRUE) {
    timetot <- end - start
    cat("\nSimulation time\n\n",
        "   Start date|time: ", format(start,usetz = TRUE), "\n",
        "   Running time: ", timetot, attr(timetot, "units"),
        file = outText, append = TRUE)
  }

  # returning results
  summ <- list(Scenarios = length(list_suit),
               Starting_porportion = starting_porportion,
               Replicates = rep, Dispersal_events = steps,
               Dispersal_kernel = dispersal_kernel,
               Kernel_spread_SD = spread, Max_dispersers = NdMax)

  if(return == "all") {
    res <- list(Summary = summ, A = Amvb[[3]], A_mean = Amvb[[1]],
                A_var = Amvb[[2]], C = Cmvb[[3]], C_mean = Cmvb[[1]],
                C_var = Cmvb[[2]])
  } else {
    if(return == "accessed") {
      res <- list(Summary = summ, A = Amvb[[3]], A_mean = Amvb[[1]],
                  A_var = Amvb[[2]], C = NULL, C_mean = NULL, C_var = NULL)
    } else {
      res <- list(Summary = summ, A = NULL, A_mean = NULL, A_var = NULL,
                  C = Cmvb[[3]], C_mean = Cmvb[[1]], C_var = Cmvb[[2]])
    }
  }

  return(res)
}
