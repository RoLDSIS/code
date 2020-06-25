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

FIGURES = ../figures/VOT-S09.pdf		\
          ../figures/erp-VOT-Ativo-S09.pdf	\
          ../figures/cv-direction-S09.pdf	\
          ../figures/cv-projections-S09.pdf	\
          ../figures/cv-errors.pdf		\
          ../figures/cv-scalograms.pdf	\

.PHONY: all
all: $(FIGURES)

### Fig. 1 of the manuscript
../figures/VOT-S09.pdf ../results/id-slope.dat: figures-id-exp.r
	Rscript figures-id-exp.r

### Fig. 2 of the manuscript
../figures/erp-VOT-Ativo-S09.pdf: figures-mean-erp.r
	Rscript figures-mean-erp.r

### Figs. 3 and 4 of the manuscript
../figures/cv-direction-S09.pdf	\
../figures/cv-projection-S09.pdf: figures-roldsis-results.r ../results/id-slope.dat
	Rscript figures-roldsis-results.r

### Cross-validation results
../results/cross-validation-phy.dat ../results/cross-validation-psy.dat: run-cv.r
	Rscript run-cv.r

### Fig. 5 of the manuscript
../figures/cv-errors.pdf: figure-cv-errors.r ../results/cross-validation-psy.dat
	Rscript figure-cv-errors.r

### Fig. 6 of the manuscript
../figures/cv-scalograms.pdf: figure-cv-scalogram.r
	Rscript figure-cv-scalogram.r

.PHONY: clean
clean:
	rm -f $(FIGURES) ../results/*.dat
