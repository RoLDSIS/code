### Function for plotting the DWT scalogram

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
source ("paths.r")

### * Load system packages
load.pkgs ("wavelets")

### * Auxiliary function for generating the color scale
scale.color <- function (nb.cols, val, min.val, range)
    min (c (nb.cols, ceiling (nb.cols * (val - min.val) / range)))

### * Palettes for scalogram
palette.wb <- c ("white", "black")
palette.bwr <- c ("blue", "white", "red")

### * Main function for plotting the scalogram
plot.scalogram <- function (dwt.obj,
                            hi.level = dwt.start.level,
                            palette = palette.wb,
                            nb.cols = 51,
                            samp.freq = eeg.sampfreq,
                            main = NA,
                            x.axis = TRUE,
                            y.axis = TRUE) {

    if (! dwt.obj@aligned )
        dwt.obj <- align (dwt.obj)

    lo.level <- length (dwt.obj@W)
    if (hi.level > lo.level)
        stop (sprintf ("hi.level must not exceed %d", lo.level))

    max.val <- -Inf
    min.val <- Inf
    for (i in seq (hi.level, lo.level)) {
        wi <- dwt.obj@W [[i]]
        max.val <- max (max.val, wi)
        min.val <- min (min.val, wi)
    }
    vn <- dwt.obj@V [[lo.level]]
    max.val <- max (max.val, vn)
    min.val <- min (min.val, vn)
    if (length (palette) == 3)
        if (max.val > abs (min.val))
            min.val <- -max.val
        else
            max.val <- -min.val
    range <- max.val - min.val

    col.pal <- colorRampPalette (palette) (nb.cols)
    width <- length (dwt.obj@W [[hi.level]])
    height <- lo.level - hi.level + 2

    margins <- c (4.6, 3.6, 2.8, 0.4)
    if (is.na (main)) {
        margins [3] <- 0
        main = ""
    }
    if (! x.axis)
        margins [1] <- 0
    if (! y.axis)
        margins [2] <- 0
    par (mar = margins)
    plot (NA, NA, xlim = c (0, width), ylim = c (0, height), bty = "n",
          xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = main)

    for (i in seq (lo.level, hi.level)) {

        wi <- dwt.obj@W [[i]]
        n <- length (wi)
        wd <- width / n

        for (j in seq (1, n)) {
            y <- lo.level - i + 1 + c (0, 0, 1, 1)
            x <- (j - 1) * wd + wd * c (0, 1, 1, 0)
            idx.col <- scale.color (nb.cols, wi [j], min.val, range)
            polygon (x, y, col = col.pal [idx.col], border = "grey")
        }

        if (i == lo.level) {
            v <- dwt.obj@V [[lo.level]]
            for (j in seq (1, n)) {
                y <- lo.level - i + c (0, 0, 1, 1)
                x <- (j - 1) * wd + wd * c (0, 1, 1, 0)
                idx.col <- scale.color (nb.cols, vn [j], min.val, range)
                polygon (x, y, col = col.pal [idx.col], border = "grey")
            }
        }

    }

    if (y.axis)
        axis (2, at = 0.5 + seq (0, lo.level - hi.level + 1),
              labels = c (sprintf ("V%d", lo.level),
                          sapply (seq (lo.level, hi.level),
                                  function (x) sprintf ("W%d", x))),
              tick = FALSE, las = 1)

    ## ** Draw time axis
    if (x.axis) {
        par (new = TRUE, mar = margins)
        plot (NA, NA, xlim = c (0, 2 * length (dwt.obj@W [[1]]) / samp.freq),
              ylim = c (0, 1), bty = "n", xaxt = "n", yaxt = "n",
              xlab = "time (s)", ylab = "")
        axis (1)
    }

}
