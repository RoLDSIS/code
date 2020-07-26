### Plot the projection of the RoLDSIS procedure for all subjects

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
source ("roldsis.r")
source ("dwt-lib.r")
source ("cross-validation.r")
source ("scalogram.r")

### * Load the system packages
load.pkgs (c ("shape", "Cairo"))

### * Load slopes of the psychometric curves
load (file.path (results.dir, "id-slope.dat"))

### * Colors for the stimuli responses (add alpha level)
cols <- col2rgb (stim.cols)
cols [, 3] <- 0.4 * cols [, 3] # increase saturation of stimulus #3
cols <- rgb (t (cols / 255), alpha = 0.5)

### * Output types
outputs <- c ("phy", "psy")
title <- list (phy = "Φ", psy = "Ψ")

### * Specify the DWT coefficients on which the cross-validation will be done
nb.wavelets <- 2 * length (dwt (rep (0, dwt.length))@W [[dwt.start.level]])
idx.wavelets <- seq (dwt.length - nb.wavelets + 1, dwt.length)

### * Settings for the position of the stimulus labels
t <- seq (0, by = 1 / eeg.sampfreq, length.out = dwt.length)
ang <- c (-0.7, -1, 0, 0.2, 0.5) * pi
label.x <- cos (ang) * 0.05
label.y <- sin (ang) * 0.8
search.lim <- round (c (0.25, 0.28) * eeg.sampfreq)
search.idx <- seq (search.lim [1], search.lim [2])

direction <- signals <- list ()

### * Loop over the cohort
for (subj in cohort) {

    ## *** Load data for that subject
    load (cv.filename (cv.exp.feature, cv.exp.type, subj))

    ## *** Generate the scalogram figure
    cairo_pdf (file.path (figures.dir, sprintf ("cv-direction-S%02d.pdf", subj)),
               width = 4, height = 5.84)

    layout (matrix (seq (1, 8), ncol = 2, byrow = TRUE),
            heights = c (1.5, 5.4, 1.5, 6.6),
            widths = c (1, 0.1))

    ### ** Loop over output types
    for (out in outputs) {

        if (out == "phy")
            Y <- phy.out [subj, ] / 200
        else
            Y <- psy.out

        ## *** Run RolSIS on a single fold and get the results
        folds <- k.folds (dwt.coefs.cv$response, dwt.coefs.cv$stimulus, 1)
        sol <- roldsis (folds$x, Y)
        dir <- sol$direction
        proj <- sol$projection
        signals [[out]] <- apply (proj, 1, function (x) vec.to.signal (x, dwt.length))

        direction [[out]] <- rbind (direction[[out]], t(dir))

        par (mar = c (0, 4, 0.75, 0) + 0.1)

        x <- vec.to.signal (dir, dwt.length)
        t <- seq (0, by = 1 / eeg.sampfreq, length.out = length (x))
        time.shift <- 0.024 * max (t)
        plot (t, x, type = "l", bty = "n", lwd = 2, las = 1,
              xaxt = "n", xlab = "", yaxt = "n", ylab = "",
              xlim = c (0, max (t)) + time.shift)
        abline (h = 0, col = "#00000080")

        par (mar = c (0, 0, 0, 0))
        plot (0, 0, type = "n", bty = "n", xlab = "", ylab = "",
              xaxt = "n", yaxt = "n")

        plot.scalogram (dwt (x), palette = palette.bwr,
                        x.axis = ifelse (out == "phy", FALSE, TRUE))

        par (mar = c (ifelse (out == "phy", 0, 4), 0, 0, 0))
        plot (0, 0, type = "n", bty = "n", xlab = "", ylab = "",
              xaxt = "n", yaxt = "n")
        text (0, 0, title [[out]], adj = c (0.5, 0.5), cex = 2)


    } ## out

    dummy <- dev.off ()

    ## *** Generate the time-domain projections figure
    cairo_pdf (file.path (figures.dir,sprintf ("cv-projections-S%02d.pdf", subj)),
               width = 4, height = 4)

    layout (matrix (seq (1, 4), ncol = 2, byrow = TRUE),
            heights = c (0.74, 1), widths = c (1, 0.1))

    ### ** Loop over output types
    for (out in outputs) {

        sig <- signals [[out]]

        y.lim <- c (min (sig), max (sig) + 5)

        par (mar = c (ifelse (out == "phy", 0, 5), 4, 0, 0) + 0.1)
        plot (0, 0, xlim = c (min (t), max (t)), ylim = y.lim, type = "n",
              las = 1, col = cols [1], bty = "n", lwd = 2,
              xlab = ifelse (out == "phy", "", "time (s)"),
              ylab = "amplitude μV",
              xaxt = ifelse (out == "phy", "n", "s"))

        for (i in seq (1, 5))
            lines (t, sig [, i], col = cols [i], lwd = 2)

        label.idx <- (search.idx [1] - 1
            + which.max (sig [search.idx, 5] - sig [search.idx, 1]))

        x.end <- rep (t [label.idx], 5)
        y.end <- sig [label.idx, ]
        x.start <- x.end + label.x
        y.start <- y.end + label.y * diff (y.lim) / 3

        Arrows (x.start, y.start, x.end, y.end, arr.adj = 1, arr.length = 0.25)
        points (x.start, y.start, cex = 3, pch = 21, bg = "white")
        text (x.start, y.start, labels = seq (1, 5), cex = 1)

        par (mar = c (ifelse (out == "phy", 0, 4), 0, 0, 0))
        plot (0, 0, type = "n", bty = "n", xlab = "", ylab = "",
              xaxt = "n", yaxt = "n")
        text (0, 0, title [[out]], adj = c (0.5, 0.5), cex = 2)

    } # out

    dummy <- dev.off ()

    ## *** Progress meter
    cat (sprintf ("\rSubject %2d", subj))
    flush (stdout ())

} # subj

### * Clean the progress meter
cat ("\n")
flush (stdout ())

### * Compose the PDF file with the results for all subjects
system (sprintf ("pdftk %s/cv-projections-S*.pdf cat output %s/cv-projections-all.pdf",
                 figures.dir, figures.dir))
system (sprintf ("pdftk %s/cv-direction-S*.pdf cat output %s/cv-direction-all.pdf",
                 figures.dir, figures.dir))

### * Compute angle between physical and psychophysical direction vectors
ang <- rep (NA, length (cohort))
for (i in cohort)
    ang  [i] <- acos (sum (direction$phy [i, ] * direction$psy [i, ])) * 180 / pi

### * Plot angles × slopes of subjects' psychometric curves
cairo_pdf (file = file.path (figures.dir,
                             sprintf ("%s-%s-slope-angle.pdf",
                                      cv.exp.feature, cv.exp.type)),
           width = 4, height = 4)
par (mar = c (5, 4, 0, 0) + 0.1)
plot (id.slope, ang, bty = "n", pch = 19, las = 1, xlim = c (2, 11),
      ylim = c (20, 70), ylab = "Φ/Ψ angle (degrees)", xlab = "slope (%/ms)")
### ** Plot
pca <- prcomp (cbind (id.slope, ang))
### Gets slope of loading/eigenvector PC1
b <- pca$rotation [2, 1] / pca$rotation [1, 1]
a <- as.numeric (pca$center [2] - b * pca$center [1])
abline (a, b)
dummy <- dev.off ()
