#' Helper function to generate interpolation values
#'
#' @export

interpolation_values <- function(transition_to_lgm, stable_lgm, lgm_to_current,
                                 stable_current, simulation_period, scenario_span) {
  ### types of clime considering current
  types_clim <- rev(c(rev(seq(0, 1, 1 / (transition_to_lgm * 1000))),
                      rep(0, (stable_lgm * 1000) - 1),
                      seq(0, 1, 1 / (lgm_to_current * 1000)),
                      rep(1, (stable_current * 1000) - 1)))


  if ((simulation_period * 1000) - 1 > length(types_clim)) {
    timesbig <- ceiling(((simulation_period * 1000) -1) / length(types_clim))
    types_clim <- rep(types_clim, timesbig)
  }

  ### position of each scenario in types of clime
  pos_scenarios <- seq(0, (simulation_period * 1000) - 1, scenario_span * 1000)
  pos_scenarios[1] <- 1

  return(list(types_clim, pos_scenarios))
}

#' Helper function to interpolate and prepare layer names
#'
#' @export

interpolation <- function(ellipsoid_model, suitability_threshold, types_clim,
                          pos_scenarios, current, lgm, current_suitability,
                          lgm_suitability, slash_type, pca_directory,
                          suitability_directory) {
  suit_name <- vector()
  for (i in 1:length(pos_scenarios)) {
    spot_val <- types_clim[pos_scenarios[i]]

    if (!spot_val %in% c(0, 1)) {
      if (!dir.exists(pca_directory)) {dir.create(pca_directory)}
      inter_fol <- paste0(pca_directory, "/", paste0("Interpolation_", i))
      dir.create(inter_fol)

      ### interpolation of PCs
      pc_inter <- list()
      for (j in 1:3) {
        pc_inter[[j]] <- (current[[j]] * (1 - spot_val)) + (lgm[[j]] * spot_val)
        filenam <- paste0(inter_fol, "/", paste0(names(current[[j]]), ".asc"))
        raster::writeRaster(pc_inter[[j]], filenam, format = "ascii")
      }
      pc_inter <- do.call(raster::stack, pc_inter)

      ### suitability projections
      suit_preds <- predict_esuitability(ellipsoid_model = ellipsoid_model,
                                         variables = pc_inter,
                                         suitability_threshold = suitability_threshold,
                                         tolerance = 1e-60)

      ### write suitability layer other scenarios
      ip_name <- paste0(suitability_directory, "/suitability_interpolation", i, ".asc")
      raster::writeRaster(suit_preds$suitability_layer, ip_name, format = "ascii")

      suit_name[i] <- gsub("/", slash_type,
                           paste0("\"", paste0(getwd(), slash_type, ip_name), "\""))

      ## truncation with barriers
      #if (!is.null(barriers)) {
      #  barrier_processing
      #}

    } else {
      if (spot_val == 0){
        suit_name[i] <- lgm_suitability
      } else {
        suit_name[i] <- current_suitability
      }
    }

    cat("\tInterpolation", i, "of", length(pos_scenarios), "has finished\n")
  }

  suit_name <- suit_name[rev(seq(1:length(suit_name)))]

  return(suit_name)
}


#' Helper function to prepare M files
#'
#' @param directory (character)
#' @param pattern (character)
#' @param crs (character)
#'
#' @export

M_preparation <- function(directory, pattern = "A_S\\d.*c$",
                          crs = "+proj=longlat +datum=WGS84 +no_defs") {
  # reading the last A_S
  mss <- list.files(directory, pattern = pattern)
  pl <- gregexpr("S\\d*", mss)
  pla <- regmatches(mss, pl)
  place <- as.numeric(gsub("S", "", unlist(pla)))

  m <- raster::raster(paste0(directory, "/", paste0("A_S", max(place), ".asc")))
  m[m[] == 0] <- NA
  m <- raster::trim(m)
  shpm <- raster::rasterToPolygons(m, dissolve = TRUE)

  sp::proj4string(shpm) <- sp::CRS(crs)

  # writing shapefile and raster of M
  rgdal::writeOGR(shpm, directory, "M", driver = "ESRI Shapefile")

  m_name <- paste(directory, "M.asc", sep = "/")
  raster::writeRaster(m, filename = m_name, format = "ascii")

  return(list(m, shpm))
}


#' Helper function to save and plot niche elllipsoid figure
#'
#' @export

save_nicheplot <- function(ellipsoid_model, suitability_threshold, variables,
                           size_proportion = 0.55, output_directory, plot) {
  el <- ellipse::ellipse(x = ellipsoid_model[[3]], centre = ellipsoid_model[[2]],
                         level = (100 - suitability_threshold) / 100)

  background <- na.omit(raster::values(variables[[1:2]]))
  background <- background[sample(nrow(background), size = 10000), ]

  xlim <- range(range(el[, 1]), range(background[, 1]))
  ylim <- range(range(el[, 2]), range(background[, 2]))

  if (plot == TRUE) {
    if (.Platform$OS.type == "unix") {
      quartz()
    } else {
      x11()
    }
    par(mar = c(4.5, 4, 0.5, 0.5), cex = 0.9)
    plot(background, xlim = xlim, ylim = ylim, pch = 1, col = "grey75",
         xlab = "PC 1", ylab = "PC 2")
    lines(el, col = "black", lwd = 2)
    points(ellipsoid_model[[1]][, 3:4], pch = 16, col = "black", cex = 0.85)
    legend("topright", legend = c("Background", "Occurrences", "Niche ellipsoid"),
           lty = c(NA, NA, 1), lwd = c(NA, NA, 2), pch = c(1, 16, NA),
           col = c("grey75", "black", "black"), bty = "n", horiz = TRUE)
  }

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


#' Helper function to save and plot M figure
#'
#' @export

save_Mplot <- function(ellipsoid_model, suitability_layer, M_polygon,
                       size_proportion = 0.55, output_directory, plot) {
  xlims <- raster::extent(suit)[1:2]
  ylims <- raster::extent(suit)[3:4]

  if (plot == TRUE) {
    if (.Platform$OS.type == "unix") {
      quartz()
    } else {
      x11()
    }

    par(mar = c(0.5, 0.5, 0.5, 0.5))
    sp::plot(M_polygon, border = "transparent", xlim = xlims, ylim = ylims)
    raster::image(suitability_layer, col = rev(terrain.colors(255)),
                  axes = FALSE, add = TRUE)
    sp::plot(M_polygon, lwd = 2, border = "blue4", add = TRUE)
    points(ellipsoid_model[[1]][, 1:2], pch = 16, cex = 0.7)
    box()

    Sys.sleep(2)
  }

  png(paste0(output_directory, "/M_in_geography.png"), width = 80, height = 80,
      units = "mm", res = 600)
  sp::plot(M_polygon, border = "transparent", xlim = xlims, ylim = ylims)
  raster::image(suitability_layer, col = rev(terrain.colors(255)),
                axes = FALSE, add = TRUE)
  sp::plot(M_polygon, lwd = 1.2, border = "blue4", add = TRUE)
  points(ellipsoid_model[[1]][, 1:2], pch = 16, cex = 0.7)
  box()
  dev.off()
}
