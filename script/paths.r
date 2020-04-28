### Support functions for paths an package loading

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

### * Load local library
source ("dwt-parameters.r")

### * Force creation of directories
force.dir.create <- function (dir)
    suppressWarnings (dir.create (dir, recursive = TRUE))

### * Needed directories
data.dir <- "../data"
figures.dir <- "../figures"
force.dir.create (figures.dir)
results.dir <- "../results"
force.dir.create (results.dir)
exp.dir <- file.path (data.dir, "Experimento")
id.dir <- file.path (data.dir, "identification")

dwt.dir <- file.path (data.dir,
                      sprintf ("%s-%d-%d",
                               dwt.filter,
                               dwt.length,
                               dwt.nb.levels))
force.dir.create (dwt.dir)

axes.dir <- file.path (dwt.dir, sprintf ("start-W%d", dwt.start.level))
force.dir.create (axes.dir)

### * Get name of subject specific directory
subject.dirname <- function (exp.feat, exp.type, subject) {
    feat.dir <- file.path (dwt.dir, exp.feat)
    type.dir <- file.path (feat.dir, exp.type)
    subj.label <- sprintf ("S%02d", subject)
    subj.dir <- file.path (type.dir, subj.label)
    force.dir.create (subj.dir)
    return (subj.dir)
}

### * Get name of the DWT results directory
dwt.coefs.stem <- "dwt-coefs"
dwt.coefs.filename <- function (exp.feat, exp.type, subject)
    file.path (subject.dirname (exp.feat, exp.type, subject),
               sprintf ("%s.dat", dwt.coefs.stem))

### * Get name of the cross-validation directory
cv.stem <- "cross-validation"
cv.filename <- function (exp.feat, exp.type, subject)
    file.path (subject.dirname (exp.feat, exp.type, subject),
               sprintf ("%s.dat", cv.stem))

### * Utility function for Loading R pakcages, installing them if necessary
load.pkgs <- function (list.of.pkgs)
    for (pkg in list.of.pkgs)
        if (! require (pkg, character.only = TRUE)) {
            install.packages (pkg)
            require (pkg, character.only = TRUE)
        }
