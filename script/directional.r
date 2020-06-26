crt.to.sph <- function (crt, adjust.last = TRUE) {
    sph <- t (apply (crt, 1, crt.to.sph.vector))
    if (adjust.last) {
        nc <- ncol (sph)
        nr <- nrow (sph)
        x <- sph [, nc]
        s <- sort (x, index.return = TRUE)
        ix <- s$ix
        y <- c ()
        for (i in seq (1, nr)) {
            idx <- seq (1, i)
            yi <- x
            yi [ix [i : 1]] <- yi [ix [i : 1]] + 2 * pi
            y <- cbind (y, yi)
        }
        idx <- which.min (apply (y, 2, function (x) max (x) - min (x)))
        sph [, nc] <- y [, idx]
    }
    return (sph)
}

crt.to.sph.vector <- function (crt) {
    n <- length (crt)
    sph <- rep (NA, n - 1)
    sph [1] <- acos (crt [1])
    s <- sin (sph [1])
    if (n > 3)
        for (i in seq (2, n - 2)) {
            if (s == 0)
                sph [i] <- 0
            else
                sph [i] <- acos (crt [i] / s)
            s <- s * sin (sph [i])
        }
    if (s == 0)
        sph [n - 1] <- 0
    else
        sph [n - 1] <- atan2 (crt [n] / s, crt [n - 1] / s)
    return (sph)
}

sph.to.crt <- function (sph)
    crt <- t (apply (sph, 1, sph.to.crt.vector))

sph.to.crt.vector <- function (sph) {
    n <- length (sph)
    crt <- rep (NA, n + 1)
    crt [1] <- cos (sph [1])
    s <- 1
    for (i in seq (2, n)) {
        s <- s * sin (sph [i - 1])
        crt [i] <- s * cos (sph [i])
    }
    crt [n + 1] <- s * sin (sph [n])
    return (crt)
}
