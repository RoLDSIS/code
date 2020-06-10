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
load.pkgs ("shape")

### * Load slopes of the psychometric curves
load (file.path (results.dir, "id-slope.dat"))

### * Colors for the stimuli responses (add alpha level)
cols <- col2rgb (stim.cols)
cols [, 3] <- 0.4 * cols [, 3] # increase saturation of stimulus #3
cols <- rgb (t (cols / 255), alpha = 0.5)

### Responses
responses <- c ("phy", "psy")

### * Specify the DWT coefficients on which the cross-validation will be done
nb.wavelets <- 2 * length (dwt (rep (0, dwt.length))@W [[dwt.start.level]])
idx.wavelets <- seq (dwt.length - nb.wavelets + 1, dwt.length)

### * Settigns for the position of the stimulus labels
t <- seq (0, by = 1 / eeg.sampfreq, length.out = dwt.length)
ang <- c (-1, 0.83, 0.3, 0.35, 0.55) * pi
label.x <- cos (ang) * 0.05
label.y <- sin (ang)
search.lim <- round (c (0.1, 0.15) * eeg.sampfreq)
search.idx <- seq (search.lim [1], search.lim [2])

direction <- list ()

### * Loop over the cohort
for (subj in cohort) {

    ### ** Loop over responses
    for (resp in responses) {

        if (resp == "phy")
            Y <- phy.resp [subj, ] / 200
        else
            Y <- psycho.resp

        ## *** Load data for that subject
        load (cv.filename (cv.exp.feature, cv.exp.type, subj))

        ## *** Run RolSIS on a single fold and get the results
        folds <- k.folds (dwt.coefs.cv$response, dwt.coefs.cv$stimulus, 1)
        sol <- roldsis (folds$x, Y)
        dir <- sol$direction
        proj <- sol$projection
        signals <- apply (proj, 1, function (x) vec.to.signal (x, dwt.length))

        direction [[resp]] <- rbind(direction[[resp]], t(dir))

        ## *** Generate the scalogram figure
        pdf (file.path (figures.dir, sprintf ("cv-direction-%s-S%02d.pdf",
                                              resp, subj)),
             width = 5 * 1.2, height = 3.5 * 1.2)

        layout (matrix (c (1, 2), ncol = 1), heights = c (1.5, 6))
        par (mar = c (0, 4, 0.75, 0) + 0.1)

        x <- vec.to.signal (dir, dwt.length)
        t <- seq (0, by = 1 / eeg.sampfreq, length.out = length (x))
        time.shift <- 0.024 * max (t)
        plot (t, x, type = "l", bty = "n", lwd = 2, las = 1, xaxt = "n",
              xlab = "", yaxt = "n", ylab = "", main = resp,
              xlim = c (0, max (t)) + time.shift)
        abline (h = 0, col = "#00000080")

        plot.scalogram (dwt (x), palette = palette.bwr)

        dummy <- dev.off ()


        ## *** Generate the time-domain projections figure
        pdf (file.path (figures.dir,
                        sprintf ("cv-projections-%s-S%02d.pdf", resp, subj)),
             width = 6, height = 3)

        par (mar = c (5, 4, 0, 0) + 0.1)

        y.lim <- c (min (signals), max (signals)+0.5)

        plot (0, 0, xlim = c (min (t), max (t)), ylim = y.lim, type = "n",
              las = 1, col = cols [1], bty = "n", lwd = 2, main = resp,
              xlab = "time (s)", ylab = bquote ("amplitude (" * mu * "V)"))
        for (i in seq (1, 5))
            lines (t, signals [, i], col = cols [i], lwd = 2)

        label.idx <- (search.idx [1] - 1
            + which.max (signals [search.idx, 5] - signals [search.idx, 1]))

        x.end <- rep (t [label.idx], 5)
        y.end <- signals [label.idx, ]
        x.start <- x.end + label.x
        y.start <- y.end + label.y * diff (y.lim) / 3

        Arrows (x.start, y.start, x.end, y.end, arr.adj = 1)
        points (x.start, y.start, cex = 4, pch = 21, bg = "white")
        text (x.start, y.start, labels = seq (1, 5), cex = 1.5)

        dummy <- dev.off ()

        ## *** Progress meter
        cat (sprintf ("\rResponse: %s    Subject %2d", resp, subj))
        flush (stdout ())

    } # resp

} # subj

### * Clean the progress meter
cat ("\n")
flush (stdout ())

### * Compose the PDF file with the results for all subjects
system (sprintf ("pdftk %s/cv-projections-*-S*.pdf cat output %s/cv-projections-all.pdf",
                 figures.dir, figures.dir))
system (sprintf ("pdftk %s/cv-direction-*-S*.pdf cat output %s/cv-direction-all.pdf",
                 figures.dir, figures.dir))

### * Compute angle between physical and psychophysical direction vectors
ang <- rep (NA, length (cohort))
for (i in cohort)
    ang  [i] <- acos (sum (direction$phy [i, ] * direction$psy [i, ])) * 180 / pi

### * Plot angles × slopes of subjects' psychometric curves
pdf (file = file.path (figures.dir,
                       sprintf ("%s-%s-slope-angle.pdf",
                                cv.exp.feature, cv.exp.type)))
plot (id.slope, ang, pch = 19, las = 1, ylab = "psy/phy angle", xlab = "slope (%/ms)",
      main = sprintf ("%s %s", cv.exp.feature, cv.exp.type),)
abline (lm (ang ~ id.slope))
dummy <- dev.off ()
