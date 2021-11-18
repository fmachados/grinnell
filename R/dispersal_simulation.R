#' Helper function to simulate species dispersal processes
#'
#' @description dispersal_simulation generates a python script to perform a
#' multistep simulation of species dispersal to reconstruct areas that have been
#' accessed based on environmental suitability and user-defined dispersal parameters.
#'
#' @param data (character) name of the csv file with all the species occurrences,
#' columns must be: species, longitude, latitude.
#' @param suit_layers (character) vector of names of suitability layers to be
#' used as distinct scenarios. If more than one, the layer names must be ordered
#' from past to current. Layer names should include parent directory if needed.
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
#' @param access_threshold (numeric) percentage to be considered when excluding
#' accessed cells with lower values. Default = 5.
#' @param output_directory (character) name of the output directory where
#' results should be written. If this directory does not exist, it will be
#' created.
#'
#' @export
#'
#' @usage
#' dispersal_simulation(data, suit_layers, dispersal_kernel = "normal",
#'                      kernel_spread = 1, max_dispersers = 4,
#'                      replicates = 10, dispersal_events = 25,
#'                      access_threshold = 5, output_directory)

dispersal_simulation <- function(data, suit_layers, dispersal_kernel = "normal",
                                 kernel_spread = 1, max_dispersers = 4,
                                 replicates = 10, dispersal_events = 25,
                                 access_threshold = 5, output_directory) {

  sl <- "/"
  dl <- "/"

  out <- gsub("\"", "", output_directory)
  report <- paste0("\"", paste0(out, "/report.txt"), "\"")
  ms <- paste0("\"", paste0(out, "/A_S"), "\"")
  msc <- paste0("\"", paste0(out, "/A_mean_S"), "\"")
  vsc <- paste0("\"", paste0(out, "/A_var_S"), "\"")

  replicates_folder <- paste0("\"", paste0(out, "/M_replicates/M_S"), "\"")

  sink(paste(out, "M_simulation.py", sep = sl))

  cat("#!usr/bin/env python\n")

  cat("# -*- coding: utf-8 -*-\n")

  cat("import os
import numpy
import linecache
import csv
import copy
import math
import time\n
os.chdir(")

  cat(output_directory)

  cat(")
outText = open(")

  cat(report)

  cat(",\"w\")\n
metadata = []
HeaderMeta = []
for j in range(6):
    metadata.append(linecache.getline(")

  cat(suit_layers[length(suit_layers)])

  cat(", j+1))
    HeaderMeta.append(metadata[j])
    metadata[j] = metadata[j].split()
    metadata[j][1] = float(metadata[j][1])

Ncol = int(metadata[0][1])
Nrow = int(metadata[1][1])

list_asc = [")


  cat(paste(suit_layers, collapse = ",\n"))

  cat("]\nwith open(")

  cat(data)

  cat(", \"r\") as f:
    reader = csv.reader(f)
    spdata = list(reader)
del spdata[0]

rep = ")

  cat(replicates)

  cat("\nsteps = ")

  cat(dispersal_events)

  cat("\nspread = ")

  cat(kernel_spread)

  if (dispersal_kernel == "normal") {
    cat("\nfat = False")
  }
  if (dispersal_kernel == "log_normal") {
    cat("\nfat = True")
  }

  cat("\nthres = ")

  cat(access_threshold)

  cat("\nNdMax = ")

  cat(max_dispersers)

  cat("\nNd_list = list(range(1, NdMax + 1))

Dk = \"Normal\"
if (fat):
    Dk = \"LogNormal\"

outText.write(\"***Parameters:\" + \"\\n\" + \"\\n\" +
              \"Suitability scenarios: \" + str(len(list_asc)) + \"\\n\" +
              \"Replicates: \" + str(rep) + \"\\n\" +
              \"Dispersal events: \" + str(steps) + \"\\n\" +
              \"Dispersal kernel: \" + Dk + \"\\n\" +
              \"Kernel spread (SD): \" + str(spread) + \"\\n\" +
              \"Maximum number of dispersers: \" + str(NdMax) + \"\\n\")

incNd = 1 / NdMax

S_list = list(numpy.arange(incNd, 1.0 + incNd, incNd))\n")

  cat("\ndef setLoc(coord):
    NWlon, NWlat = round(metadata[2][1], 12), round((metadata[1][1]*metadata[4][1] + metadata[3][1]), 12)
    NWvertex = [NWlat, NWlon]
    loc_val, loc = [], []
    loc_val.append(NWvertex[0] - coord[1])
    loc_val.append(coord[0] -NWvertex[1])

    for i in range(2):
        loc.append(round(loc_val[i] / metadata[4][1]))
    return loc

def setPop(spdata):
    C = [[0]*Ncol for i in range(Nrow)]
    point = []
    size_l = round(len(spdata)*0.50)
    rowset = numpy.random.choice(range(1, len(spdata)), size = size_l, replace = False)
    for i in range(len(rowset)):
        point = setLoc([float(spdata[rowset[i]][1]), float(spdata[rowset[i]][2])])
        C[point[0]][point[1]] = 1
    return C

def updateA(A, A_now):
    for row in range(Nrow):
        for col in range (Ncol):
            if (A_now[row][col] != 0):
                A[row][col] = A[row][col] + A_now[row][col]
    return A

def updateC(C, A_now, S):
    for row in range(Nrow):
        for col in range (Ncol):
            if (A_now[row][col] != 0):
                if (S[row][col] > 0):
                    C[row][col] = C[row][col] + A_now[row][col]
    return C

def newS(C, S):
    for row in range(Nrow):
        for col in range (Ncol):
            if (S[row][col] <= 0):
                    C[row][col] = 0
    return C")

  cat("\n\nstart = time.time()
\nlist_rep = [[] for r in range(rep)]
\nlist_rep_C = [[] for r in range(rep)]\n
for s in range(len(list_asc)):
    S_raw = numpy.loadtxt(list_asc[s], skiprows = 6)
    S_data = []
    for row in range(Nrow):
        for col in range(Ncol):
            S_data.append(S_raw[row][col])
    S =  [S_data[Ncol * i : Ncol * (i + 1)] for i in range(Nrow)]

    mean_A = [[0]*Ncol for i in range(Nrow)]
    var_A = [[0]*Ncol for i in range(Nrow)]
    bin_A = [[0]*Ncol for i in range(Nrow)]

    for r in range(rep):
        if (s == 0):
            C = setPop(spdata)
            A = [[0]*Ncol for i in range(Nrow)]
            for row in range(Nrow):
                for col in range(Ncol):
                    if (C[row][col] == 1):
                        A[row][col] = 1
        else:
            A = list_rep[r]
            C = newS(list_rep_C[r], S)

        for t in range(steps):
            A_now = [[0]*Ncol for i in range(Nrow)]
            for row in range(Nrow):
                for col in range(Ncol):
                    if (C[row][col] >= 1):
                        if (NdMax != 1):
                            unkD = True
                            while unkD:
                                for i in range(len(S_list)):
                                    if (S[row][col] <= S_list[i]):
                                        Nd = Nd_list[i]
                                        unkD = False
                                        break
                        else:
                            Nd = NdMax
                        for i in range(Nd):
                            outside = True
                            while outside:
                                theta = numpy.random.uniform(low = 0., high = 2.) * math.pi
                                if (fat):
                                    rad = numpy.random.lognormal(mean = 0, sigma = spread)
                                    d_lat = round(rad * math.sin(theta))
                                    d_lon = round(rad * math.cos(theta))
                                    if (0 <= row + d_lat < Nrow and 0 <= col + d_lon < Ncol):
                                        row_now = int(row + d_lat)
                                        col_now = int(col + d_lon)
                                        outside = False
                                else:
                                    rad = numpy.abs(numpy.random.normal(loc = 0, scale = spread))
                                    d_lat = round(rad * math.sin(theta))
                                    d_lon = round(rad * math.cos(theta))
                                    if (0 <= row + d_lat < Nrow and 0 <= col + d_lon < Ncol):
                                        row_now = int(row + d_lat)
                                        col_now = int(col + d_lon)
                                        outside = False

                            A_now[row_now][col_now] += 1

            A = updateA(A, A_now)
            C = updateC(C, A_now, S)
        list_rep_C[r] = copy.deepcopy(C)
        list_rep[r] = copy.deepcopy(A)

        for row in range(Nrow):
            for col in range(Ncol):
                mean_A[row][col] = mean_A[row][col] + list_rep[r][row][col]
    for row in range(Nrow):
        for col in range(Ncol):
            if (S[row][col] == -9999):
                mean_A[row][col] = -9999
            else:
                mean_A[row][col] = mean_A[row][col] / rep

    MeanOut = open(")

  cat(msc)

  cat(" + str(s+1) + \".asc\",\"w\")
    for j in range(6):
        MeanOut.write(HeaderMeta[j])
    for row in range(Nrow):
        for col in range(Ncol):
            MeanOut.write(str(mean_A[row][col]) + \" \")
        if (row != Nrow - 1):
                MeanOut.write(\"\\n\")
    MeanOut.close()

    list_sq = [[[0]*Ncol for i in range(Nrow)] for r in range(rep)]

    for r in range(rep):
        for row in range(Nrow):
            for col in range(Ncol):
                list_sq[r][row][col] = (list_rep[r][row][col] - mean_A[row][col])**2

    for r in range(rep):
        for row in range(Nrow):
            for col in range(Ncol):
                var_A[row][col] = var_A[row][col] + list_sq[r][row][col]

    for row in range(Nrow):
        for col in range(Ncol):
            if (S[row][col] == -9999):
                var_A[row][col] = -9999
            else:
                var_A[row][col] = var_A[row][col] / (rep - 1)

    VarOut = open(")

  cat(vsc)

  cat(" + str(s+1) + \".asc\",\"w\")
    for j in range(6):
        VarOut.write(HeaderMeta[j])
    for row in range(Nrow):
        for col in range(Ncol):
            VarOut.write(str(var_A[row][col]) + \" \")
        if (row != Nrow - 1):
                VarOut.write(\"\\n\")
    VarOut.close()

    Abin_list = []
    for row in range(Nrow):
        for col in range(Ncol):
            Abin_list.append(mean_A[row][col])

    clean_Abin_list = []
    for i in range(len(Abin_list)):
        if (Abin_list[i] > 0):
            clean_Abin_list.append(Abin_list[i])

    p5 = numpy.percentile(clean_Abin_list, q = thres)

    for row in range(Nrow):
        for col in range(Ncol):
            if (mean_A[row][col] < 0):
                bin_A[row][col] = -9999
            elif (mean_A[row][col] < p5 and mean_A[row][col] >= 0):
                bin_A[row][col] = 0
            else:
                bin_A[row][col] = 1

    AbinOut = open(")

  cat(ms)

  cat(" + str(s+1) + \".asc\",\"w\")
    for j in range(6):
        AbinOut.write(HeaderMeta[j])
    for row in range(Nrow):
        for col in range(Ncol):
            AbinOut.write(str(bin_A[row][col]) + \" \")
        if (row != Nrow - 1):
                AbinOut.write(\"\\n\")
    AbinOut.close()\n")

  cat("\nend = time.time()


outText.write(\"\\n***Simulation time: \" + str(round(end - start, 3)) + \" sec\")
outText.close()")

  sink()
}
