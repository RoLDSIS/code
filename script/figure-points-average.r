source ("paths.r")

load.pkgs ("limSolve")

fold <- function (x, n) {
    ret <- c ()
    i <- 1
    x <- x [sample (seq (1, nrow (x))), ]
    while (i <= nrow (x)) {
        idx <- seq (i, min (nrow (x), i + n - 1))
        if (length (idx) == 1)
            ret <- rbind (ret, x [idx, ])
        else
            ret <- rbind (ret, colMeans (x [idx,]))
        i <- i + n
    }
    return (ret)
}

results <- data.frame (subject = integer (),
                       fold.size = integer (),
                       nb.points = integer (),
                       rms = numeric ())


for (subj in cohort) {

    load (cv.filename ("VOT", "Ativo", subj))

    s <- dwt.coefs.cv$stimulus
    r <- dwt.coefs.cv$response

    for (n in seq (1, 10)) {

        s2 <- r2 <- c ()

        for (i in seq (1, 5)) {
            idx <- which (s == i)
            ri <- fold (r [idx,], n)
            s2 <- c (s2, rep (i, nrow (ri)))
            r2 <- rbind (r2, ri)
        }

        if (length (s2) < ncol (r2))
            break

        x <- cbind (r2, rep (1, nrow (r2)))
        y <- c (0, 0.05, 0.5, 0.95, 1) [s2]
        b <- Solve (x, y)

        results <- rbind (results,
                          data.frame (subject = subj,
                                      fold.size = n,
                                      nb.points = length (s2),
                                      rms = sqrt (mean ((y - x %*% b)^2))))
    }

}

m <- aggregate(rms ~ fold.size, results, mean)
std <- aggregate(rms ~ fold.size, results, sd)

pdf (file = file.path (figures.dir, "points-average.pdf"),
     width = 5, height = 5)

par (mar = c (5, 4, 0, 0) + 0.1)

plot (m$fold.size[1:6], m$rms[1:6], pch = 19, cex = 2, ylim = c(0,.4),
      xlim = c (0.5,6.5), bty = "n", las = 1,
      ylab = "RMS prediction error",
      xlab = "number of points per grouping")

for (i in seq (1,6))
    lines (rep (std$fold.size[i], 2), m$rms [i] + c(-1, 1) * std$rms [i], lwd = 3)

dummy <- dev.off ()
