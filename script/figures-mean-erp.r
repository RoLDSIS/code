### Plot per subject figures with averaged ERP for each stimulus
###
### The plotting is done only for the the experiment and the electrode
### selected in the cross-validation procedure.

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
source ("remove-spikes.r")

### * Load system package
load.pkgs ("wavelets")

### * Plotting parameters
vdist <- 10
vscale.x <- 0.23
vscale.x.width <- 0.01
vscale.y <- -42
vscale.y.range <- 10

### * File stem for the plotting
file.stem <- sprintf ("erp-%s-%s-", cv.exp.feature, cv.exp.type)

### * Time basis for plotting
t <- seq (0, by = 1 / eeg.sampfreq, length.out = dwt.length)

### * Loop over the cohort
for (subj in cohort) {

    ## ** Load the subject data
    load (cv.filename (cv.exp.feature, cv.exp.type, subj))

    ## ** Open PDF file
    pdf (file = file.path (figures.dir,
                           sprintf ("%sS%02d.pdf", file.stem, subj)),
         width = 7, height = 4)
    par (mar = c (5, 0.1, 0.1, 0.1))

    ## ** Get the averaged signals for the current subject
    x <- dwt.coefs.cv$mean

    ## ** Low pass filtering for detection of P2 peak
    x.lpf <- x
    for (i in seq (1, 5)) {
        d <- dwt (x.lpf [i, ])
        ## Eliminate bands W1 to W5
        for (j in seq (1, 5))
            d@W [[j]] <- 0 * d@W [[j]]
        x.lpf [i, ] <- idwt (d)
    }

    ## ** Get extension of the plot
    max.v <- max (x [1, ])
    min.v <- min (min (x [5, ] - vdist * 4), vscale.y - vscale.y.range)

    ## ** Loop over the stimuli
    for (i in seq (1, 5)) {

        ## *** Get vertical shift
        v.shift <- -vdist * (i - 1)

        ## *** Plot the signal. Stimuli #1 starts the plot
        if (i == 1)
            plot (t, remove.spikes (x [1, ]), type = "l", bty = "n",
                  col = stim.cols [1], ylim = c (min.v, max.v),
                  xlim = c (min (t) - 0.05, max (t)),
                  xlab = "time (s)", yaxt = "n", ylab = "")
        else
            lines (t, remove.spikes (x [i, ])
                   + v.shift, col = stim.cols [i])

        ## *** Plot markers for P1 peak
        ## This is experimental code and has been commented out.
        ## mx <- max (x.lpf [i, ])
        ## idx <- which.max (x.lpf [i, ])
        ## lines (rep (t [idx], 2),
        ##        v.shift + mx + c (1.5, 3.5) * vscale.y.range / 10,
        ##        lwd = 2)

    }

    ## ** Draw stimuli number
    points (rep (-0.025, 5), seq (0, -4) * vdist, cex = 4.5)
    text (rep (-0.025, 5), seq (0, -4) * vdist, adj = c (0.5, 0.5),
          labels = seq (1, 5), cex = 1.5)

    ## ** Draw vertical scale
    lines (vscale.x + c (0, -rep (vscale.x.width, 2), 0),
           vscale.y + c (0, 0, -rep (vscale.y.range, 2)))
    text (vscale.x, vscale.y - vscale.y.range / 2, adj = c (0, 0.5),
          bquote (.(vscale.y.range) * mu * V))

    ## ** Close the PDF file
    dummy <- dev.off ()

    ## ** Progress meter
    cat (sprintf ("\rSubject %2d", subj))
    flush (stdout ())

}

### * Clean the progress meter
cat ("\n")
flush (stdout ())

### * Compose the combo file with all subjects
system (sprintf ("pdftk %s/%sS*.pdf cat output %s/%sall.pdf",
                 figures.dir, file.stem, figures.dir, file.stem))
