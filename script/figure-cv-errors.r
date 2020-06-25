### Analyse and plot results of cross-validation

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

### * Load local libraries
source ("paths.r")
source ("chisq-to-normal.r")
source ("compare-methods.r")

### * Load system packages
load.pkgs (c ("lme4", "lmerTest", "merTools", "emmeans", "Cairo"))

## *** Open the PDF file
cairo_pdf (file = file.path (figures.dir,"cv-errors.pdf"),
           width = 5, height = 7.5)

layout (matrix (seq (1, 4), ncol = 2, byrow = TRUE),
        heights = c (0.84, 1), widths = c (1, 0.1))

output <- c ("phy", "psy")

title <- list (phy = "Φ", psy = "Ψ")

### * Loop over output types
for (out in output) {

    ## ** Load results of cross-validation
    load (file.path (results.dir,
                     sprintf ("cross-validation-%s.dat", out)))

    ## ** Statistical analysis

    ## *** Add column for normal transformation of errors
    cv.df$sse.test.norm <- rep (NA, nrow (cv.df))

    ## *** Get lambda for chisq to normal transformation
    lambda <- c ()
    for (n in cv.nb.folds) {
        lambda [n] <- chisq.to.normal (5 * n)
        idx <- which (cv.df$nb.folds == n)
        cv.df$sse.test.norm [idx] <- cv.df$sse.test [idx] ^ lambda [n]
    }

    ## *** Run statistical analysis
    fm <- lmer (sse.test.norm ~ method * nb.folds + (1 | subject), cv.df)
    show (anova (fm))
    show (ranova (fm))
    show (pairs (emmeans (fm, "method")))

    ## ** Plot the results

    ## *** Confidence intervals

    alpha <- 0.95

    itv <- aggregate (sse.test.norm ~ method * nb.folds, cv.df, function (x) NA)
    lbd <- lambda [itv$nb.folds]

    n <- itv$nb.folds * 5
    CI <- qt (1 - (1 - alpha) / 2, n - 1) / sqrt (n)

    m <- aggregate (sse.test.norm ~ method * nb.folds, cv.df, mean)$sse.test.norm
    s <- aggregate (sse.test.norm ~ method * nb.folds, cv.df, sd)$sse.test.norm

    itv$mean <- m ^ (1 / lbd)
    itv$inf <- (m - s * CI) ^ (1 / lbd)
    itv$sup <- (m + s * CI) ^ (1 / lbd)

    min.v <- 0.02
    max.v <- 0.14

    ## *** Plotting parameters
    pchs <- seq (21, 24)
    cols <- c ("red", "blue", "gold3", "green4")

    par (mar = c (ifelse (out == "phy", 0, 4), 4, 1, 0.1))

    y.lab <- "mean squared error"

    ## *** Start plot without plotting
    plot (itv$mean, ylim = c (min.v, max.v), las = 1, log = "y",
          ylab = y.lab, bty = "n", xaxt = "n",
          xlab = ifelse (out == "phy", "", "number of folds"),
          type = "n")

    ## *** Plot regions for folds
    n <- length (methods)
    for (i in seq (1, length (cv.nb.folds)))
        polygon (i * n + c (-n + 0.5, -n + 0.5, 0.5, 0.5),
                 c (min.v, max.v, max.v, min.v), border = NA,
                 col = ifelse (i %% 2 == 1, "white", "#00000020"))
    if (out == "psy")
        axis (1, at = seq (0, length (cv.nb.folds) - 1) * n + (n + 1) / 2,
              labels = cv.nb.folds)

    ## *** Plot data for each fold
    for (i in seq (1, nrow (itv)))
        lines (c (i, i), c (itv$inf [i], itv$sup [i]), col = cols [itv$method [i]],
               lwd = 4)
    points (itv$mean, ylim = c(min (itv$inf), max (itv$sup)), cex = 2,
            pch = pchs [itv$method], bg = cols [itv$method])

    ## *** Legend
    if (out == "psy")
        legend ("bottom", ins = 0.05, pch = pchs, pt.bg = cols, pt.cex = 1.5,
                bg = "white", legend = names (methods))

    par (mar = c (ifelse (out == "phy", 0, 4), 0, 0, 0))
    plot (0, 0, type = "n", bty = "n", xlab = "", ylab = "",
          xaxt = "n", yaxt = "n")
    text (0, 0, title [[out]], adj = c (0.5, 0.5), cex = 2.5)


} # out

## ** Close PDF file
dummy <- dev.off ()
