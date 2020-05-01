### Remove spikes from the averaged evoked repsonses, based on the
### elimination of outlier DWT components.

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

### * Load the local library
source ("paths.r")

### * Load the system packages
load.pkgs (c ("EnvStats", "wavelets"))

### * Find outliers using Rosner's generalized extreme Studentized deviate test
which.outliers <- function (x) {
    ## ** Start with one possible outlier
    k <- 1
    ## ** Grow set until the k-th sample is not an outlier
    while (TRUE) {
        t <- rosnerTest (x, k = k)
        if (! t$all.stats$Outlier [k])
            break
        k <- k + 1
    }
    ## ** Return value
    if (k == 1)
        return (NULL)
    else
        return (t$all.stats$Obs.Num [1 : (k - 1)])
}

### * Remove spikes and return "clean" signal
remove.spikes <- function (x, max.level = NULL) {
    ## ** Comoute the DWT
    d <- dwt (x)
    ## ** Use all available W bands of the DWT, if max level is not specified
    if (is.null (max.level))
        max.level <- length (d@W)
    ## ** Iterate over the W bands
    for (i in seq (1, max.level)) {
        ## *** Replace outliers with the median value of the band
        w <- d@W [[i]]
        idx <- which.outliers (w)
        w [idx] <- median (w)
        d@W [[i]] <- w
    }
    ## ** Return value
    return (idwt (d))
}
