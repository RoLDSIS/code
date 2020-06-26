### Basic functions for running bootstrap with RoLDSIS

### * Load the local library
source ("roldsis.r")

### Load system package
load.pkgs ("boot")

### * Statistic function for boot
###
### This function, which is used as the "statistic" argument of boot, must
### provide two arguments:
###    input: A list of matrix.  Each element must correspond to one of the
###           experimental stimuli.  Each line of the matrices correspond to
###           a trial and the columns correspond to a feature (wavelet
###           coefficients)
###   output: The response for the regression.  Must be a vector whose length is
###           equal to the number of stimuli.

boot.fun <- function (input, output = NULL) {
    ## ** Initialize the input matix for roldsis
    x <- c ()
    ## ** Compute the grand average for each stimulus
    for (i in seq (1, length (input)))
        x <- rbind (x, colMeans (input [[i]]))
    ## ** Do the RoLDSIS regression
    reg <- roldsis (x, output)
    ## ** Return the result
    return (reg$dir)
}

### * Random values generation function
###
### This function do the sampling with replacement inside the set of reponses
### for each stimulus.  It is used as the "ran.gen" argument of boot.

ran.gen.fun <- function (input, mle) {
    ## ** Initialize return value
    retval <- list ()
    ## ** Sample trials for each stimulus
    for (i in seq (1, length (input))) {
        n <- nrow (input [[i]])
        retval [[i]] <- input [[i]] [sample (seq (1, n), replace = TRUE), ]
    }
    ## ** Return value
    return (retval)
}

### * Bootstrap function for the RoLDSIS procedure

roldsis.boot <- function (input, output, boot.rep)
    return (boot (input, boot.fun, boot.rep, sim = "parametric",
                  ran.gen = ran.gen.fun, output = output))
