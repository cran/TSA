`plot1.acf` <-
function (x, ci = 0.95, type = "h", xlab = "Lag", ylab = NULL, 
    ylim = NULL, main = NULL, ci.col = "blue", ci.type = c("white", 
        "ma"), max.mfrow = 6, ask = Npgs > 1 && dev.interactive(), 
    mar = if (nser > 2) c(3, 2, 2, 0.8) else par("mar"), oma = if (nser > 
        2) c(1, 1.2, 1, 1) else par("oma"), mgp = if (nser > 
        2) c(1.5, 0.6, 0) else par("mgp"), xpd = par("xpd"), 
    cex.main = if (nser > 2) 1 else par("cex.main"), verbose = getOption("verbose"), 
    ...) 
{
    ci.type <- match.arg(ci.type)
    if ((nser <- ncol(x$lag)) < 1) 
        stop("x$lag must have at least 1 column")
    if (is.null(ylab)) 
        ylab <- switch(x$type, correlation = "ACF", covariance = "ACF (cov)", 
            partial = "Partial ACF")
    if (is.null(snames <- x$snames)) 
        snames <- paste("Series ", if (nser == 1) 
            x$series
        else 1:nser)
    with.ci <- ci > 0 && x$type != "covariance"
    with.ci.ma <- with.ci && ci.type == "ma" && x$type == "correlation"
    clim0 <- if (with.ci) 
        qnorm((1 + ci)/2)/sqrt(x$n.used)
    else c(0, 0)
    Npgs <- 1
    nr <- nser
    if (nser > 1) {
        sn.abbr <- if (nser > 2) 
            abbreviate(snames)
        else snames
        if (nser > max.mfrow) {
            Npgs <- ceiling(nser/max.mfrow)
            nr <- ceiling(nser/Npgs)
        }
        opar <- par(mfrow = rep(nr, 2), mar = mar, oma = oma, 
            mgp = mgp, ask = ask, xpd = xpd, cex.main = cex.main)
        on.exit(par(opar))
        if (verbose) {
            cat("par(*) : ")
            str(par("mfrow", "cex", "cex.main", "cex.axis", "cex.lab", 
                "cex.sub"))
        }
    }
    for (I in 1:Npgs) for (J in 1:Npgs) {
        iind <- (I - 1) * nr + 1:nr
        jind <- (J - 1) * nr + 1:nr
        if (verbose) 
            cat("Page [", I, ",", J, "]: i =", paste(iind, collapse = ","), 
                "; j =", paste(jind, collapse = ","), "\n")
        for (i in iind) for (j in jind) if (max(i, j) > nser) {
            frame()
            box(col = "light gray")
        }
        else {
            clim <- if (with.ci.ma && i == j) 
                clim0 * sqrt(cumsum(c(1, 2 * x$acf[-1, i, j]^2)))
            else clim0
            if (is.null(ylim)) {
                ymin <- min(c(x$acf[, i, j], -clim), na.rm = TRUE)
                ymax <- max(c(x$acf[, i, j], clim), na.rm = TRUE)
                ylim <- c(ymin, ymax)
            }
            plot(x$lag[, i, j], x$acf[, i, j], type = type, xlab = xlab, 
                ylab = if (j == 1) 
                  ylab
                else "", ylim = ylim, ...)
            abline(h = 0)
            if (with.ci && ci.type == "white") 
                abline(h = c(clim, -clim), col = ci.col, lty = 2)
            else if (with.ci.ma && i == j) {
                lines(x$lag[, i, j], clim, col = ci.col, lty = 2)
                lines(x$lag[, i, j], -clim, col = ci.col, lty = 2)
            }
            title(if (!is.null(main)) 
                main
            else if (i == j) 
                snames[i]
            else paste(sn.abbr[i], "&", sn.abbr[j]), line = if (nser > 
                2) 
                1
            else 2)
        }
        if (Npgs > 1) {
            mtext(paste("[", I, ",", J, "]"), side = 1, line = -0.2, 
                adj = 1, col = "dark gray", cex = 1, outer = TRUE)
        }
    }
}
