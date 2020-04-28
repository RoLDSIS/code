### Support functions for DWT objects

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

### * Load the system library
load.pkgs ("wavelets")

### * Transform a DWT object into a vector of coeficients object

dwt.to.coefs <- function (d) {
    ## ** Accumulate the W bands
    coefs <- c ()
    w <- names (d@W)
    for (n in w)
        coefs <- c (coefs, d@W [[n]])
    ## ** Add the first V band
    v <- tail (names (d@V), n = 1)
    coefs <- c (coefs, d@V [[v]])
    ## ** Return a list containing the wavelet coefficients and the list of
    ## levels that were used to generate the coefficients
    list (coefs = coefs, levels = c (w, v))
}

### * Transform a coefficients object into the equivalent a DWT object

coefs.to.dwt <- function (coefs) {
    ## ** Initialize an empty (i.e. all coefficients equal to zero) DWT object
    n <- length (coefs$coefs)
    nb.levels <- length (coefs$levels) - 1
    d <- dwt (rep (0, n), n.levels = nb.levels)
    ## ** Extract the W bands
    start <- 1
    for (i in seq (1,  nb.levels)) {
        n <- n / 2
        w <- coefs$levels [i]
        d@W [[w]] <- matrix (coefs$coefs [start : (start + n - 1)], ncol = 1)
        start <- start + n
    }
    ## ** Extract the V band
    v <- coefs$levels [nb.levels + 1]
    d@V [[v]] <- matrix (coefs$coefs [start : (start + n - 1)], ncol = 1)
    d
}

### *** Transform a shortened vector of coefficients into the equivalent a
### DWT object.  This function assumes that the input vector contains only
### the low-frequency components of the wavelet transform.  The argument n
### specifies the desired length of the output DWT object.

vec.to.dwt <- function (x, n = length (x)) {
    ## ** Hi-frequency components are at the beginning of the vector.  Assume
    ## that they have value equal to zero.
    x <- c (rep (0, n - length (x)), x)
    ## ** Generate the appropriate argument for function coefs.to.dwt
    d <- list (coefs = x,
               levels = dwt.to.coefs (dwt (x, n.levels = dwt.nb.levels))$levels)
    coefs.to.dwt (d)
}

### *** Transform a shortened vector of coefficients into the equivalent
### temporal series.  The argument n (a power of 2) specifies the length of
### the signal.

vec.to.signal <- function (x, n = length (x))
    idwt (vec.to.dwt (x, n))
