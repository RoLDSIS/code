### Find the optimal lambda that transforms a chi-squared distributed
### variable X, such that the transformed variable X^lambda is optimally
### normal.
###
### Reference:
### Hawkins DM, Wixley RAJ (1986) A note on the transformation of chi-squared
### variables to normality.  Am Stat 40(4):296. DOI:
### 10.2307/2684608. (http://dx.doi.org/10.2307/2684608).

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

### * Returns the exponent of the tranformation
chisq.to.normal <- function (df) {

    ## ** High df should use lambda = 1/3
    ## if (df >= 30)
    ##     return (1 / 3)

    ## ** Find optimum lambda
    ret <- optimize (
        ## *** Implement the transformation x^lambda
        function (lambda, df) {
            ## *** Find the empirical initial values of the mean and the stderr
            ## of the transformed normal variable
            m <- qchisq (0.5, df) ^ lambda
            s <- mean (c (m - qchisq (pnorm (-1), df)  ^ lambda,
                          qchisq (pnorm (1), df) ^ lambda - m))
            ## *** Find the optimum mean and the stderr
            ret <- optim (c (m, s),
                ## **** Compute the Kullback-Leibner divergence
                function (par, df) {
                    ## ***** Extract the parameters
                    m <- par [1]
                    s <- par [2]
                    ## ***** Implement integral equation
                    ret <- integrate (function (x) {
                        ret <- rep (0, length (x))
                        ## ****** Compute prob. distribution of xhi^2 variable
                        dx <- dchisq (x, df)
                        ## ****** Idem for normal, transformed variable
                        ## This uses the formula for transformation of
                        ## probability densities, with the derivative
                        ## of the transformation function.
                        dn <- (dnorm (x ^ lambda, mean = m, sd = s)
                               * (lambda * x ^ (lambda - 1)))
                        ## ****** Avoid singularities for zero values
                        i <- which (dx > 0 & dn > 0)
                        ret [i] <- dx [i] * (log (dx [i]) - log (dn [i]))
                        return (ret)
                    }, 0, 5 * df,
                    ## ****** Avoid integration convergence problems
                    ## This seems to be harmful to the computations,
                    ## at least until lambda = 30.
                    stop.on.error = FALSE)
                    return (ret$value)
                }, df = df)
            return (ret$value)
        },
        ## ** Initial interval for searching lambda
        ## (see Hawkins & Wixley, 1986)
        c (0.2, 1/3),
        df = df)

    return (ret$minimum)

}
