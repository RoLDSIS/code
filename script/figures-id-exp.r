### Script to generate the results of the phonemic identification experiemnt

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
#source ("experiments.r")

### * Load the system libraries
load.pkgs (c ("robustbase", "signal", "shape", "Cairo"))

### * Set see of the random number generator (for reproducible graphics)
set.seed (0)

### * Values of the physical feature in the stimuli continua
### VOT: initial = -52ms, final = 17ms
### Formantes: initial = 533 Hz, final = 1387 Hz
continuum <- list (VOT = function (x) -52 + 68 * (x - 1) / 199,
                   Formantes = function (x) 533 + 854 * (x - 1) / 199)

### * Parameters for plotting
x.label <- list (VOT = "VOT (ms)", Formantes = "F2 – F1 (Hz)")
y.label <- list (VOT = "/ta/", Formantes = "/pɛ/")
shift.label <- list (VOT = 5, Formantes = 63)

### * Scale for the graphics in the PDF file
scale <- 0.6

### * Initialize vector for slope of psychometric identification curve
id.slope <- c ()

### * Loop over the experimental features
### for (feature in experiment.features) {
for (feature in "VOT") {

    ## ** Loop over the subjects
    for (subj in cohort) {

        ## *** Load the subject data
        load (file.path (id.dir, feature, sprintf ("S%02d.dat", subj)))

        ## *** Get the experimental values
        stim <- continuum [[feature]] (id.result$stimulus)
        resp <- id.result$response

        ## *** Open the output PDF file and set margins for the plot
        cairo_pdf (file = file.path (figures.dir, sprintf ("%s-S%02d.pdf",
                                                     feature, subj)),
                   width = scale * 12, height = scale * 7)
        par (mar = c (5, 4, 1, 0.1))

        ## *** Plot the experimental points, with jitter in the y axis
        plot (stim, resp + rnorm (length (resp), sd = 0.03), pch = 16,
              col = "#00000050", bty = "n", las = 1, xaxt = "n", yaxt = "n",
              xlab = x.label [[feature]], yaxt = "n",
              ylab = sprintf ("percentage of %s response", y.label [[feature]]))
        axis (2, at = c (0, 0.5, 1), labels = c (0, 50, 100), las = 1)

        ## *** Fit the psychometric model
        df <- data.frame (stim = stim, resp = resp)
        fit <- glmrob (resp ~ stim, df, family = "binomial")

        ## *** Store slope of psychometric curve (in %/ms)
        id.slope <- c (id.slope, 25 * coefficients (fit) [2])

        ## *** Plot the predict psychometric curve
        pred.stim <- seq (min (stim), max (stim), length.out = 200)
        pred.resp <- predict (fit, data.frame (stim = pred.stim),
                              type = "response")
        lines (pred.stim, pred.resp, col = "#00000080", lwd = 3.5)

        ## *** Get the values for the horizontal positions of the five stimuli
        ## **** Cope with the case where the psychometric curve did not reach
        ## the value 0.05 for the first stimulus in the continuum.  The third
        ## stimulus in the continuum was taken as stim2
        if (pred.resp [1] > psy.out [2])
            interp.stim <- c (min (stim),
                              continuum [[feature]] (3),
                              interp1 (pred.resp, pred.stim, psy.out [3:4]),
                              max (stim))
        else
            interp.stim <- c (min (stim),
                              interp1 (pred.resp, pred.stim, psy.out [2:4]),
                              max (stim))
        axis (1, at = sapply (interp.stim, function (x) sprintf ("%.0f", x)))
        points (interp.stim, psy.out, cex = 3, pch = 21, bg = stim.cols)

        ## *** Cope with superposition of stim#2 on stim#1 and stim#4 on stim#5
        if ((interp.stim [2] - interp.stim [1])
            / (interp.stim [5] - interp.stim [1]) < 0.23)
            reflect.stim2 <- 1
        else
            reflect.stim2 <- -1
        if ((interp.stim [5] - interp.stim [4])
            / (interp.stim [5] - interp.stim [1]) < 0.2)
            reflect.stim4 <- -1
        else
            reflect.stim4 <- 1

        ## *** Add the labels for the stimuli
        pos.beg.x <- (interp.stim
            + shift.label [[feature]] * c (1, reflect.stim2, -2,
                                           reflect.stim4, -1))
        pos.end.x <- interp.stim
        pos.beg.y <- psy.out + c (1, 1, 0, -1, -1) * 0.2
        pos.end.y <- psy.out

        pos.end.x <- pos.end.x - 0.35 * (pos.end.x - pos.beg.x)
        pos.end.y <- pos.end.y - 0.35 * (pos.end.y - pos.beg.y)

        Arrows (pos.beg.x, pos.beg.y, pos.end.x, pos.end.y)

        points (pos.beg.x, pos.beg.y, cex = 6, pch = 21, bg = "white")
        text (pos.beg.x, pos.beg.y, cex = 2, adj = c (0.5, 0.5))

        ## *** Close the PDF file
        dummy <- dev.off ()

        ## ** Progress meter
        cat (sprintf ("\rFeature: %8s    Subject %2d", feature, subj))
        flush (stdout ())

    }

}

### * Clean the progress meter
cat ("\n")
flush (stdout ())

### * Save slopes of the psychometric curves
save (id.slope, file = file.path (results.dir, "id-slope.dat"))

### * Compose the PDF file with the results for all subjects
system (sprintf ("pdftk %s/*-S*.pdf cat output %s/all.pdf",
                 figures.dir, figures.dir))

