### Plot the scalograms from cross-validation results
### CAVEAT: The code below only works when there are 4 methods in the
### cross-validation procedure.

### This program is part of RoLDSIS
###
### Copyright (C) 2020  Rafael Laboissière
### Copyright (C) 2020  Adrielle de Carvalho Santana
### Copyright (C) 2020  Hani Camille Yehia
###
### This program is free software: you can redistribute it and/or modify it
### under the terms of the GNU General Public License as published by the
### Free Software Foundation, either version 3 of the License, or (at your
### option) any later version.
###
### This program is distributed in the hope that it will be useful, but
### WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
### General Public License for more details.
###
### You should have received a copy of the GNU General Public License along
### with this program.  If not, see <http://www.gnu.org/licenses/>.

### * Load the local libraries
source ("paths.r")
source ("dwt-lib.r")
source ("compare-methods.r")
source ("scalogram.r")

### * Load the system library
load.pkgs ("Cairo")

output <- c ("phy", "psy")
title <- list (phy = "Φ", psy = "Ψ")

### * Open the PDF file
cairo_pdf (file = file.path (figures.dir, "cv-scalograms.pdf"),
           width = 8, height = 4)

## ** Specify the panels
layout (matrix (seq (1, 10), nrow = 2, byrow = TRUE),
        widths = c (1.26, rep (1, 3), 0.2), heights = c (0.91, 1))

### * Counter for panels
panel <- 0

### * Loop over output types
for (out in output) {

    ## ** Loop over the methods
    for (method in names (methods)) {

        ## *** Increase counter
        panel <- panel + 1

        ## ** Load the results of the cross-validation procedure
        load (file.path (results.dir, sprintf ("cross-validation-%s.dat", out)))

        ## *** Cumulate the coefficients values per wavelet
        cf <- 0
        for (subj in cohort) {
            ## **** For RoLDSIS use regression coeeficients for reduced data sets
            ## For the other methods, use the coefficoents obtained from the
            ## cross-validation with 3 folds.
            if (method == "RoLDSIS") {
                cf.subj <- cv.results [[subj]] [[1]] [[method]] $ coef.full
                cf.subj <- cf.subj / sum (cf.subj ^ 2)
                cf <- cf + cf.subj ^ 2
            } else {
                cf.subj <- cv.results [[subj]] [[1]] [[method]] $ coef.cv
                cf.subj <- cf.subj / sum (cf.subj ^ 2)
                cf <- cf + cf.subj ^ 2
            }
        }

        cf <- sqrt (cf / length (cohort))

        ## *** Plot panel
        plot.scalogram (vec.to.dwt (cf, dwt.length),
                        main = ifelse (panel < 5, method, NA),
                        x.axis = ifelse (panel > 4, TRUE, FALSE),
                        y.axis = ifelse (panel == 1 | panel == 5, TRUE, FALSE))

    }

    par (mar = c (ifelse (panel > 4, 5, 0), 0, ifelse (panel > 4, 0, 3), 0))
    plot (0, 0, type = "n", bty = "n", xlab = "", ylab = "",
          xaxt = "n", yaxt = "n")
    text (0, 0, title [[out]], adj = c (0.5, 0.5), cex = 2)

} # out

### * Close the PDF file
dummy <- dev.off ()
