### Run de cross-validation on the whole set of subjects

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
source ("compare-methods.r")
source ("dwt-parameters.r")

outputs <- c ("phy", "psy")

### * Loop over output types
for (out in outputs) {

    ## ** Initialize output variables
    cv.df <- data.frame (method = c (), subject = c (), nb.folds = c (),
                         sse.train = c (), sse.test = c ())
    cv.results <- list ()

    ## ** Loop over the cohort
    for (subj in cohort) {

        ## *** Initialize subject slot
        cv.results [[subj]] <- list ()

        ## *** Load data for that subject
        load (cv.filename (cv.exp.feature, cv.exp.type, subj))

        ## *** Loop over the required numbers of folds
        for (i in seq (1, length (cv.nb.folds))) {

            ## *** Get number of folds
            n <- cv.nb.folds [i]
            cat (sprintf ("subject: %02d    nb.folds: %d\n", subj, n))

            ## *** Get the AER and the output values
            x <- dwt.coefs.cv$response

            if (out == "phy")
                Y <- phy.out [subj, ] / 200
            else
                Y <- psy.out

            y <- Y [dwt.coefs.cv$stimulus]

            ## ** Run cross-validations
            r <- compare.methods (x, y, n)
            cat ("\n")

            ## ** Store the results
            cv.results [[subj]] [[i]] <- r

            ## ** Compose the output data frame
            for (j in names (r))
                cv.df <- rbind (cv.df,
                                data.frame (method = j,
                                            subject = subj,
                                            nb.folds = n,
                                            sse.train = r [[j]]$sse.train,
                                            sse.test = mean (r [[j]]$sse.test)))

            ## ** Flush the progress meter
            flush (stdout ())

        }

    } # subj

    ## ** Save results
    save (file = file.path (results.dir,
                            sprintf ("cross-validation-%s.dat", out)),
          cv.results, cv.df)

} # resp
