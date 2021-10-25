# RoLDSIS – Regression on Low-Dimension Spanned Input Space

## Introduction

The present repository contains the data and the code for analysing the
data related to a publication involving the RoLDSIS method.  RoLDSIS stands
for “Regression on Low-Dimension Spanned Input Space”.  The data for an
electroencephalography (EEG) experiment involving phonemic identification
is included in the repository.  The scripts in the R language, for
processing and analysing the data, as well as producing the figures that
appear in the following article:

> Santana, A., Barbosa A., Yehia H., and Laboissière R. (2021) A
dimension reduction technique applied to regression on high dimension,
low sample size neurophysiological data sets. BMC Neurosci 22(1). DOI:
[10.1186/s12868-020-00605-0](http://dx.doi.org/10.1186/s12868-020-00605-0).

## Installation

The code can be obtained with the following command:

```
git clone https://github.com/RoLDSIS/code RoLDSIS
```

N.B.: Developers with write access right to the repository should do,
instead:

```
git clone git@github.com:RoLDSIS/code.git RoLDSIS
```

## Running the analysis and producing the figures

The system packages `R`, `make`, and `pdftk` are needed.  The necessary
R packages are installed as required by the scripts.

In GNU/Linux systems, the analysis of the data and the generation of the
figures are done with the following commands :

```
cd RoLDSIS/script
make
```

## Description of the contents of this repository

### Scripts containing supporting functions and general parameters definition

* [paths.r](script/paths.r): Definition pf paths for saving figures and results
* [dwt-parameters.r](script/dwt-parameters.r): Definition of DWT parameters
* [geom-lib.r](script/geom-lib.r): Geometry supporting function
* [dwt-lib.r](script/dwt-lib.r): DWT supporting functions
* [scalogram.r](script/scalogram.r): Plot the DWT scalogram
* [remove-spikes.r](script/remove-spikes.r): Remove spikes from signals, using DWT
* [chisq-to-normal.r](script/chisq-to-normal.r): Convert Chi-squre to normal distribution

### RoLDSIS functions

* [roldsis.r](script/roldsis.r)

### Cross-validation for RolDSIS, LASSO, Ridge Regression and SPLS

* [cross-validation.r](script/cross-validation.r): Basic cross-validation function
* [compare-methods.r](script/compare-methods.r): Regression methods
* [run-cv.r](script/run-cv.r): Main script for running the cross-validation

### Generation of figures

* [figure-cv-scalogram.r](script/figure-cv-scalogram.r)
* [figure-cv-errors.r](script/figure-cv-errors.r)
* [figures-mean-erp.r](script/figures-mean-erp.r)
* [figures-roldsis-results.r](script/figures-roldsis-results.r)
* [figures-id-exp.r](script/figures-id-exp.r)
* [figure-trials-observation.r](script/figure-trials-observation.r)

### Data

* [VOT/Snn.dat](data/identification/VOT/): Results from the phonemic identification task
* [Snn/cross-validation.dat](data/la8-2048-8/VOT/Ativo/): Raw EEG data

## Description of the generated figures:

* `psy-VOT-Snn.pdf`: Results of phonemic identification experiment for subject nn
* `erp-VOT-Ativo-Snn.pdf`: Average ERPs for the five stimuli for subject nn
* `cv-direction-Snn.pdf`: RoldSIS result represented in scalogram for subject nn
* `cv-projections-Snn.pdf`: RoLDSIS projected responses for subject nn
* `cv-errors.pdf`: Cross-validation errors for RoLDSIS, LASSO, Ridge Regression, and SPLS
* `cv-scalograms.pdf`: Histograms for the regression methods of regressed coefficients
* `trials-observation.pdf`: RMS prediction error for all subjects using different number of averaged points

## Documentation

* [doc/RoLDSIS.pdf](doc/RoLDSIS.pdf): Early version of the paper describing
  the method

## Licensing conditions

The files in this repository are made available under the conditions of the
[GNU Public License, version 3 or later](COPYING).  No warranties.  The
RolDSIS has been submitted for publication. If you use the software or the
data of this repository, please give credit.  Stay tuned for more
information on the eventual acceptation for publication of the manuscript .


## Authors

* Copyright © 2020 [Rafael Laboissière](https://github.com/rlaboiss)
* Copyright © 2020 [Adrielle de Carvalho Santana](https://github.com/Adrielle-Santana)
* Copyright © 2020 [Hani Camille Yehia](https://github.com/haniyehia)
