### General DWT parameters

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

### * Valid subjects
cohort <- seq (1, 11)

### * Mother wavelet
dwt.filter <- "la8"

### * Length of the analysis window
dwt.length <- 2048

### * Number of levels of the DWT
dwt.nb.levels <- 8

### * Start of region of interest
dwt.start.level <- 5

### * Number of eeg channels
eeg.channels <- 17

### * Length of baseline (in number of samples)
baseline.length <- 750

### * EEG acquisition rate
eeg.sampfreq <- 5000

### * Do shrinkage for denoising?
do.shrinkage <- TRUE

### * Cross-validation settings
cv.exp.feature <- "VOT"
cv.exp.type <- "Ativo"
cv.electrode <- "TP9"
cv.nb.folds <- seq (3, 6)

### * Levels of the psychometric curve for determining the five stimuli
psy.out <- c (0, 0.05, 0.5, 0.95, 1)

### Physical responses of each subject for vot continuum
vot <- read.csv (file = "../data/identification/VOT/stim-VOT.csv")
phy.out <- as.matrix (vot [, c ("stim1", "stim2", "stim3", "stim4", "stim5")])

### * Color of stimuli used in cros-validation plots
stim.cols <- c ("red", "orange", "gray", "green3", "blue")
