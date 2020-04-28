### Cross-validation functions

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

### * Load the nloptr library
load.pkgs ("nloptr")

### * Folds the output values into k groups
###
### Arguments:
### x: matrix with the input vectors for the regression, one line per trial
### y: vector with the output values for the regression.  The folds are
###    organized in terms of the unique values present in y.
### nb.folds: number of folds
###
### Return values:
### x: averaged trials
### y: associated output values
### id: associated fold id

k.folds <- function (x, y, nb.folds) {

    ## ** Get unique values in the y vector
    y.vals <- sort (unique (y))
    nb.vals <- length (y.vals)

    ## ** Initialize return variables
    retval <- list (x = c (), y = c (), id = c ())

    ## ** Loop over the output values
    for (i in seq (1, nb.vals)) {

        ## *** Get corresponding indices in y vector
        idx.val <- which (y == y.vals [i])

        ## *** Initialize loop variables
        amount <- length (idx.val)
        start <- 1

        ## *** Loop over the folds
        for (j in seq (1, nb.folds)) {

            ## **** Get the number of samples in the ith group
            count <- round (amount / (nb.folds - j + 1))

            ## **** Compute the associated row indices
            idx.fold <- seq (start, start + count - 1)

            ## **** Add to the x and y output variables
            retval$x <- rbind (retval$x, colMeans (x [idx.val [idx.fold], ]))
            retval$y <- c (retval$y, y.vals [i])
            retval$id <- c (retval$id, j)

            ## **** Update the loop variables
            amount <- amount - count
            start <- start + count
        }

    }

    ## ** Return value
    return (retval)

}

### * Compute the regression erros for the cross validation
###
### Arguments:
### params: the list of continuous regularization parameters
### k: Integer regularization parameter
### folds: data folds (returned by k.folds)
### regression.fct: Regression function for the specific method.
###                 input: (x, y, params, k)
###                 output: fitted model object
### predict.fct: Prediction function
###              input: (fm, newx)
###              output: vector with predicted y
###
### Return value:
### Matrix [n,m] with prediction errors on test folds, columns correspond
### to the folds

cv.regression.errors <- function (params,
                                  k,
                                  folds,
                                  regression.fct,
                                  predict.fct) {

    ## ** Initialize output matrix
    err <- c ()

    ## ** Loop over folds
    for (i in unique (folds$id)) {

        ## *** Get test fold
        idx <- which (folds$id == i)
        x.test <- folds$x [idx, ]
        y.test <- folds$y [idx]

        ## *** Get train fold
        idx <- which (folds$id != i)
        xs <- folds$x [idx, ]
        ys <- folds$y [idx]

        y.train <- x.train <- c ()
        for (s in unique (ys)) {
            y.train <- c (y.train, s)
            idx <- which (ys == s)
            if (length (idx) == 1)
                x.train <- rbind (x.train, xs [idx,])
            else
                x.train <- rbind (x.train, colMeans (xs [idx,]))
        }

        ## *** Get fitted model
        fm <- regression.fct (x.train, y.train, params, k)

        ## *** Cumulate regression error for current test fold
        err <- rbind (err, as.vector ((y.test - predict.fct (fm, x.test)) ^ 2))

    }

    ## ** Return value
    return (err)

}


### * Objective function for the optimization procedure

objective.fct <- function (params,
                           k,
                           folds,
                           regression.fct,
                           predict.fct) {

    err <- cv.regression.errors (params,
                                 k,
                                 folds,
                                 regression.fct,
                                 predict.fct)

    return (sum (err))

}

### * Do the cross validation
###
### Arguemnts:
### x: input matrix. trials in rows and features in columns
### y: output vector
### nb.folds: number of folds
### regression.fct: Function that returns the fitted object.
###              Must have signature regression.fct (x, y)
### predict.fct: Function that returns the predicted values.
###              Must have signature predict.fct (fm, newx),
###              where fm is the fitted object
### coefficients.fct: Function that returns the fitted coefficients.
###              Must have signature coefficients.fct (fm)
### init.params: Initial values for the continuous regularization parameters
### lb.params: lower bounds for the continuous regularization parameters
### ub.params: upper bounds for the continuous regularization parameters
### k.max.fct: Function for determining the maximum value of k, the
###            integer regularization parameter.  The number of folds is given
###            as input.  Should return NULL if there is no integer parameter.
###
### Values:
### params: optimal regularization parameters
### k: optimal interger regularization parameter
### coef: regression coefficients of optimal model
### sse: sum of squared errors of regression for all folds
### sse.test: squared errors of regression for the test folds

cross.validation <- function (x, y,
                              nb.folds,
                              regression.fct,
                              predict.fct,
                              coefficients.fct,
                              init.params,
                              lb.params,
                              ub.params,
                              k.max.fct) {

    ## ** Get the folds
    folds <- k.folds (x, y, nb.folds)

    ## ** Initalize loop parameters
    opt.sse <- opt.params <- c ()
    k.max <- k.max.fct (nb.folds)

    if (is.null (k.max)) {

        ## *** There is no integer parameter k
        fm <- regression.fct (folds$x, folds$y, NULL, NULL)
        opt.k <- opt.params <- NULL

    } else {

        ## *** Initialize variables
        opt.sse <- Inf
        opt.params <- init.params

        ## *** Loop over the integer parameter k
        ## k should not exceed 4, because, in some casesm, there will be
        ## only 5 points for doing the regression.
        for (k in seq (1, min (4, k.max))) {

            ## **** Optimize contiuous regularization parameters
            sol <- nloptr (x0 = opt.params,
                           eval_f = objective.fct,
                           opts = list (algorithm = "NLOPT_LN_COBYLA",
                                        xtol_rel = 1e-4,
                                        ftol_rel = 1e-6,
                                        maxeval = 1000000),
                           lb = lb.params,
                           ub = ub.params,
                           k = k,
                           folds = folds,
                           regression.fct = regression.fct,
                           predict.fct = predict.fct)

            ## **** Clean progress meter
            cat ("\n")

            ## **** Get SSE and end loop if it increased
            sse <- sol$objective
            if (sse > opt.sse)
                break

            ## **** Update loop variables
            opt.sse <- sse
            opt.params <- sol$solution
            opt.k <- k

        }

        ## *** Get optimal fitted model
        fm <- regression.fct (folds$x, folds$y, opt.params, opt.k)

    }

    ## ** Get the optimal coefficients
    coef.cv <- as.vector (coefficients.fct (fm))

    ## ** Compute total SSE
    sse.train <- sum ((folds$y - predict.fct (fm, folds$x)) ^ 2)

    ## ** Compute SS errors for the test folds
    sse.test <- cv.regression.errors (opt.params, opt.k, folds,
                                      regression.fct, predict.fct)

    ## ** Apply optimal model to whole data with a single fold
    folds <- k.folds (folds$x, folds$y, 1)
    fm <- regression.fct (folds$x, folds$y, opt.params, opt.k)
    coef.full <- as.vector (coefficients.fct (fm))
    sse.full <- sum ((folds$y - predict.fct (fm, folds$x)) ^ 2)

    ## ** Return value
    return (list (params = opt.params,
                  k = opt.k,
                  coef.cv = coef.cv,
                  sse.train = sse.train,
                  sse.test = sse.test,
                  coef.full = coef.full,
                  sse.full = sse.full))

}
