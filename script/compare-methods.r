### Definition of regression methods

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
source ("cross-validation.r")
source ("roldsis.r")

### * Load the system libraries
load.pkgs (c ("glmnet", "spls"))

### * Auxiliary functions for glmnet-base LASSO and Ridge Regression

glmnet.predict.fct <- function (fm, newx) {
    n <- length (fm$lambda)
    p <- predict (fm, newx)
    return (p [, n])
}

glmnet.coefficients.fct <- function (fm) {
    n <- length (fm$lambda)
    coef <- coefficients (fm)
    return (unname (coef [2 : dim (coef) [1], n]))
}

### * Methods

methods <- list (

    RoLDSIS = list (
        regression = function (x, y, params, k)
            return (roldsis (x, y)),
        predict = predict.roldsis,
        coefficients = coefficients.roldsis,
        init.params = NULL,
        lb.params = NULL,
        ub.params = NULL,
        k.max.fct = function (n) NULL
    ),

    LASSO = list (
        regression = function (x, y, params, k)
            return (glmnet (x, y, alpha = 1, lambda = params,
                            standardize = FALSE)),
        predict = glmnet.predict.fct,
        coefficients = glmnet.coefficients.fct,
        init.params = 0,
        lb.params = 0,
        ub.params = Inf,
        k.max.fct = function (n) 1
    ),

    Ridge = list (
        regression = function (x, y, params, k)
            return (glmnet (x, y, alpha = 0, lambda = params,
                            standardize = FALSE)),
        predict = glmnet.predict.fct,
        coefficients = glmnet.coefficients.fct,
        init.params = 0,
        lb.params = 0,
        ub.params = Inf,
        k.max.fct = function (n) 1
    ),

    SPLS = list (
        regression = function (x, y, params, k) {
            fm <- spls (x, y, eta = params, K = k)
            cat (sprintf ("\rk = %d  eta = %f  err = %f",
                          k, params,
                          sum ((y - predict (fm)) ^ 2)))
            return (fm)
        },
        predict = predict,
        coefficients = coefficients,
        init.params = 0.5,
        lb.params = 1e-5,
        ub.params = 1 - 1e-5,
        k.max.fct = function (n) {(n - 1) * 5 - 1}
    )

    ## spcr = list (
    ##     regression = function (x, y, params, k) {
    ##         fm <- spcr (x, y, k,
    ##                     lambda.B = params [1],
    ##                     lambda.gamma = params [2],
    ##                     w = params [3],
    ##                     xi = params [4])
    ##         cat (sprintf ("\rk = %d  lambda.B = %f  lambda.gamma = %f  w = %f  xi = %f  err = %f",
    ##                       k, params [1], params [2], params [3], params [4],
    ##                       sum ((y - fm$gamma0 - x %*% fm$loadings.B %*% fm$gamma) ^ 2)))
    ##         flush (stdout ())
    ##         return (fm)
    ##     },
    ##     predict = function (fm, newx)
    ##         return (fm$gamma0 + newx %*% fm$loadings.B %*% fm$gamma),
    ##     coefficients = function (fm)
    ##         return (fm$loadings.B %*% fm$gamma),
    ##     init.params = c (0, 0, 0.1, 0.01),
    ##     lb.params = c (0, 0, 0, 0),
    ##     ub.params = c (Inf, Inf, 1, 1),
    ##     k.max.fct = function (n) {(n - 1) * 5 - 1}
    ## )

)

compare.methods <- function (x, y, nb.folds) {

    retval <- list ()

    for (m in names (methods)) {

        cat (sprintf ("\nRunning method %s\n", m))
        flush (stdout ())

        meth <- methods [[m]]

        retval [[m]] <- cross.validation (x, y,
                                          nb.folds,
                                          meth$regression,
                                          meth$predict,
                                          meth$coefficients,
                                          meth$init.params,
                                          meth$lb.params,
                                          meth$ub.params,
                                          meth$k.max.fct)

    }

    return (retval)
}
