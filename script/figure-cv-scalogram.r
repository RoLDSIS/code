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

### * Load the results of the cross-validation procedure
load (file.path (results.dir, "cross-validation.dat"))

### * Open the PDF file
pdf (file = file.path (figures.dir, "cv-scalograms.pdf"),
     width = 6, height = 4)

### * Specify the panels
layout (matrix (c (1, 2, 3, 4), nrow = 2),
        widths = c (1, 0.82), heights = c (0.68, 1))

### * Counter for panels
panel <- 1

### * Loop over the methods
for (method in names (methods)) {

    ## ** Cumulate the coefficients values per wavelet
    cf <- 0
    for (subj in cohort) {
        ## *** For RoLDSIS use regression coeeficients for reduced data sets
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

    ## ** Plot panel
    plot.scalogram (vec.to.dwt (cf, dwt.length) , main = method,
                    x.axis = ifelse (panel == 1 | panel == 3, FALSE, TRUE),
                    y.axis = ifelse (panel > 2, FALSE, TRUE))

    ## ** Increase counter
    panel <- panel + 1

}

### * Close the PDF file
dummy <- dev.off ()
