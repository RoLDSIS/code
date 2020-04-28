### RoLDSIS functions

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
source ("geom-lib.r")

### * Regression function

roldsis <- function (data.points, response) {

    ## ** Get the subspace embeding the mean responses for each stimulus
    axes <- sub.space (data.points)

    ## ** Compute the subject's grand mean
    M <- colMeans (data.points)

    B <- cbind (project (data.points, M, axes), rep (1, nrow (data.points)))

    sol <- solve (B, response)
    dir <- axes %*% sol [1 : (length (sol) - 1)]
    dir <- dir / sqrt (sum (dir ^ 2))

    pred <- B %*% sol
    proj <- (data.points %*% dir) %*% t (dir)

    for (i in seq (1, nrow (proj)))
        proj [i, ] <- M + proj [i, ]

    list (direction = dir,
          projection = proj,
          prediction = pred,
          solution = sol,
          mean = M,
          axes = axes)

}

### * Prediction function

predict.roldsis <- function (sol, newx) {

    B <- cbind (project (newx, sol$mean, sol$axes), rep (1, nrow (newx)))
    return (B %*% sol$solution)

}

### * Coefficients function

coefficients.roldsis <- function (sol)
    return (sol$direction)
