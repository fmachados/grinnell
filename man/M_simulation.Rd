% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/M_simulation.R
\name{M_simulation}
\alias{M_simulation}
\title{Simulation of species accessible areas (M) Python version}
\usage{
M_simulation(data, current_variables, barriers = NULL, project = FALSE,
             projection_variables, scale = TRUE, center = TRUE,
             dispersal_kernel = "normal", kernel_spread = 1,
             max_dispersers = 4, suitability_threshold = 5,
             replicates = 10, dispersal_events = 20,
             access_threshold = 5, simulation_period = 50,
             stable_lgm = 7, transition_to_lgm = 100, lgm_to_current = 7,
             stable_current = 13, scenario_span = 1, output_directory,
             overwrite = FALSE)
}
\arguments{
\item{data}{(character) name of the csv file with all the occurrences used to
run the simulation; columns must be: species, longitude, latitude.}

\item{current_variables}{(character) name of the folder where environmental
variables representing current conditions are. T least two variables are
needed and they must be in ascii format.}

\item{barriers}{(character) "optional" name of the raster file representing
barriers for the species dispersal. Barriers must have the same projection,
format, and extent than variables. The only values allowed in this layer are
1 and NA; 1 = barrier. Default = NULL.}

\item{project}{(logical) whether or not to project environmental suitability
to past scenarios. The projection is done to the scenario defined by
\code{projection_variables} and to any other scenario resulting from
interpolations between current and past conditions. If TRUE,
arguments \code{projection_variables}, \code{simulation_period},
\code{stable_current}, \code{stable_lgm}, \code{transition_to_lgm},
\code{lgm_to_current}, and \code{scenario_span} need to be defined.
Default = FALSE.}

\item{projection_variables}{(character) name of the folder where environmental
variables representing the "Last Glacial Maximum" scenario. Variable names,
projection, and extent of these layers must be the same than those in
\code{current_variables}.}

\item{scale}{(logical) whether or not to scale variables while performing
principal component analyses.}

\item{center}{(logical) whether or not to center variables while performing
principal component analyses.}

\item{dispersal_kernel}{(character) dispersal kernel (dispersal function)
used to simulate the movement of the species. Options are: "normal",
"log_normal". Default = "normal".}

\item{kernel_spread}{(numeric) standard deviation for the
\code{dispersal_kernel}. Default = 2.}

\item{max_dispersers}{(numeric) maximum number of dispersers that depart from
each colonized pixel. Depending on suitability this number will decrease in
less suitable areas.}

\item{suitability_threshold}{value (percentage) to be used as threshold
for suitability; default = 5. Below this value environments are considered
unsuitable.}

\item{replicates}{(numeric) number of times to repeat the simulation
per scenario. Default = 10.}

\item{dispersal_events}{(numeric) number of dispersal events to happen per
scenario. Default = 25.}

\item{access_threshold}{(numeric) percentage of frequency values to be
considered as highly unlikely to have been visited during process of
dispersal. Default = 5.}

\item{simulation_period}{(numeric) time in thousands of years for the
complete period of simulation.}

\item{stable_lgm}{(numeric) time in thousands of years for the period when
glacial (LGM) conditions are assumed to be relatively stable.}

\item{transition_to_lgm}{(numeric) time in thousands of years for the
transition period from current-like (interglacial) to glacial (LGM) climatic
conditions.}

\item{lgm_to_current}{(numeric) time in thousands of years for the
transition period from glacial (LGM) to current-like (interglacial) climatic
conditions.}

\item{stable_current}{(numeric) time in thousands of years for the period when
current-like (interglacial) conditions are assumed to be relatively stable.}

\item{scenario_span}{(numeric) time in thousands of years that have to pass
for changing the scenario. Default = 1 (one thousand years).}

\item{output_directory}{(character) name of the output directory to be created
in which all results will be written.}

\item{overwrite}{(logical) whether or not to overwrite the
\code{output_directory} if it already exists. Default = FALSE.}
}
\value{
The complete set of results derived from data preparation and the simulation
is written in \code{output_directory}. These results include:
\itemize{
\item occurrences found in suitable areas in the scenario where the simulation
started.
\item a folder containing results from the PCA performed
\item a folder containing results from the preparation of suitability layer(s)
\item accessible areas as raster layers (value 1 = accessed)
\item other raster layers representing statistics of accessibility:
mean and variance
\item accessible areas as a shapefile (only accessed areas)
\item a plot representing the accessible areas and the occurrences
\item a simple report from the simulation process
}
}
\description{
M_simulation generates an area that has been potentially
accessible to a species based on simulations of dispersal events determined
by environmental suitability and user-defined parameters. The dispersal
simulation is performed using Python (see details).
}
\details{
An external dependency for this function is Python >= 3.6 and some of its
libraries: os, numpy, linecache, csv, copy, math, time. We recommend to
install Python 3 via Anaconda which will include all the required libraries.

A principal component analysis is performed with \code{current_variables}.
Then the three first principal components are used to calculate the
suitability layer used in dispersal simulations. Values of suitability are
derived from an ellipsoid envelope model created based on occurrence records
and principal components. The ellipsoid model is used because it is a simple
yet reliable representation of a species ecological niche that does not
require a background or pseudo-absences.

If \code{barriers} are used, suitability values in the areas where barriers
exist become zero. This is, populations cannot establish there and dispersal
will be truncated unless dispersal abilities defined by arguments
\code{dispersal_kernel} and \code{kernel_spread}, allow the species to
overpass the barriers.

If \code{project} = TRUE, the simulation will run on a set of scenarios
representing glacial-interglacial climate conditions. This set of scenarios
are constructed based on interpolations between environmental conditions in
\code{current_variables} and \code{projection_variables}. The later set of
variables must represent Last Glacial Maximum conditions. Interpolations
are linear and depend on the distance between the two initial set of
conditions and other parameter defined in \code{simulation_period},
\code{stable_current}, \code{stable_lgm}, \code{transition_to_lgm},
\code{lgm_to_current}, and \code{scenario_span}.
}
\examples{
\dontrun{
# preparing data and directories for examples
## directories
tempdir <- file.path(tempdir(), "msim")
dir.create(tempdir)

cvariables <- paste0(tempdir, "/variables")
dir.create(cvariables)

lgmvariables <- paste0(tempdir, "/LGM")
dir.create(lgmvariables)

## data
data("records", package = "grinnell")
variables <- terra::rast(system.file("extdata/variables.tif",
                                     package = "grinnell"))
variables_lgm <- terra::rast(system.file("extdata/variables_lgm.tif",
                                         package = "grinnell"))
names(variables_lgm) <- names(variables)
barrier <- terra::rast(system.file("extdata/barrier.tif",
                                   package = "grinnell"))

## writing data in temporal directories
occ <- paste0(tempdir, "/records1.csv")
write.csv(records, occ, row.names = FALSE)

barr <- paste0(tempdir, "/barrier1.asc")
terra::writeRaster(barrier, filename = barr)

vnam <- paste0(cvariables, "/var_", 1:6, ".asc")
terra::writeRaster(variables, filename = vnam)

vnam <- paste0(lgmvariables, "/var_" 1:6, ".asc")
terra::writeRaster(variables_lgm, filename = vnam)

odir1 <- paste0(tempdir, "/eg_msim1")
odir2 <- paste0(tempdir, "/eg_msim2")
odir3 <- paste0(tempdir, "/eg_msim3")

# simulations
## example in current scenario
M_simulation(data = occ, current_variables = cvariables,
             max_dispersers = 2, replicates = 3, dispersal_events = 5,
             output_directory = odir1)

## example under changing climatic conditions (starting from the past)
M_simulation(data = occ, current_variables = cvariables,
             project = TRUE, projection_variables = lgmvariables,
             kernel_spread = 2, max_dispersers = 2,
             replicates = 3, dispersal_events = 25,
             simulation_period = 25, stable_lgm = 7,
             transition_to_lgm = 3, lgm_to_current = 3,
             stable_current = 7, scenario_span = 3,
             output_directory = odir2)

## example under changing conditions, considering dispersal barriers
M_simulation(data = occ, current_variables = cvariables,
             barriers = barr, project = TRUE,
             projection_variables = lgmvariables,
             kernel_spread = 2, max_dispersers = 2,
             replicates = 3, dispersal_events = 25,
             simulation_period = 25, stable_lgm = 7,
             transition_to_lgm = 3, lgm_to_current = 3,
             stable_current = 7, scenario_span = 3,
             output_directory = odir3)
}
}
