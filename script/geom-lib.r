### Support functions for geometric manipulations

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

### * Load system library
load.pkgs ("pracma") # For nullspace()

### * Given a N×M matrix x, with N points in an M-dimensional space, return
### N-1 orthonormal vectors that define the sub-space spanned by those N
### points.

sub.space <- function (x) {
    nb.points <- nrow (x)
    nb.dims <- ncol (x)
    vectors <- matrix (NA, nrow = nb.points - 1, ncol = nb.dims)
    vectors [1, ] <- normalize (x [2, ] - x [1, ])
    if (nb.points > 2)
        for (i in seq (3, nb.points)) {
            ns <- nullspace (vectors [1 : (i - 2), , drop = FALSE])
            v <- x [i, , drop = FALSE] - x [1, , drop = FALSE]
            v.ns <- ns %*% t (v %*% ns)
            vectors [i - 1, ] <- normalize (v.ns)
        }
    t (vectors)
}

### * Utility function for normalizing a vector

normalize <- function (v) {
    v / sqrt (sum (v ^ 2))
}

### * Project a series of points P [n×m matrix] onto a line defined by an
### origin M [m vector] and a direction vector V [m vector].  n is the
### number of points and m is the dimension of the space.
### Function sweep subtracts the mean M, line by line of P

project <- function (P, M, V)
    sweep (P, 2, M) %*% V

