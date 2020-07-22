
### * Load local libraries
source ("paths.r")
source ("dwt-lib.r")
source ("roldsis-bootstrap.r")
source ("directional.r")

set.seed (1234)

### * Load system package
load.pkgs (c ("wavelets", "MASS", "Cairo"))

output.type <- c ("phy", "psy")
title <- list (phy = "Φ", psy = "Ψ")

output.cols <- c ("red", "blue")

boot.rep <- 100

ef <- "VOT"
et <- "Ativo"

boot.result <- pca.all <- list ()

for (subj in cohort) {

    ## Defines which matrix of physical responses to use
    output <- list (psy = psy.out,
                    phy = phy.out [subj, ])

    boot.result [[subj]] <- list ()

    load (cv.filename (ef, et, subj))

    input <- list ()

    ## Start stim and electrodes loops to compose the big matrix
    ## to be used in the bootstrap procedure.
    for (stim in seq (1, 5)) {
        idx.stim <- which (dwt.coefs.cv$stimulus == stim)
        input [[stim]] <- dwt.coefs.cv$response [idx.stim, ]
    }

    ## Run bootstrap
    for (ot in output.type)
        boot.result [[ot]] <- roldsis.boot (input, output [[ot]],
                                            boot.rep)

    phy <- boot.result $ phy $ t
    psy <- boot.result $ psy $ t

    nb.items <- nrow (psy)

    ang <- crt.to.sph (rbind (phy, psy))
    pca <- prcomp (ang)
    pca.all [[subj]] <- pca

    cols <- rep (output.cols, rep (nb.items, 2))
    cols.alpha <- rgb (t (col2rgb (cols) / 255), alpha = 0.5)

    pca.var <- sum (pca$sdev ^2)

    cairo_pdf (file.path (figures.dir, sprintf ("lda-pca-S%02d.pdf", subj)),
               width = 5, height = 5)
    par (mar = c (5, 4, 0, 0) + 0.1)
    plot (pca$x [, 1], pca$x [, 2], pch = 19, col = cols.alpha,
          las = 1, asp = 1, bty = "n",
          xlab = sprintf ("PC1 (%.1f%%)", 100 * pca$sdev [1] ^2 / pca.var),
          ylab = sprintf ("PC2 (%.1f%%)", 100 * pca$sdev [2] ^2 / pca.var))
    class <- rep (c ("psy", "phy"), rep (nb.items, 2))
    df <- data.frame (pca$x [, 1 : 2], class = class)
    z <- lda (class ~ ., df)
    cf <- coefficients (z)
    abline (0, -cf [1] / cf [2], lwd = 2, col = "#00000060")
    legend ("bottomleft", ins = 0.05, pch = 19, col = output.cols,
            legend = c (title$phy, title$psy))
    dummy <- dev.off ()

    cat (sprintf ("Subject %02d: %3d misclassifications\n",
                  subj, length (which (predict (z) $ class != class))))

}

var12 <- c ()
for (i in seq (1, length (pca.all))) {
    sdev <- pca.all [[i]]$sdev
    var12 [i] <- (sdev [1] ^ 2 + sdev [2] ^ 2) / sum (sdev ^ 2)
}

cat ("\n")
flush (stdout ())
