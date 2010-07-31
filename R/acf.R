`acf` <-
function (x, lag.max = NULL, type = c("correlation", "covariance", 
    "partial"), plot = TRUE, na.action = na.fail, demean = TRUE, 
    drop.lag.0 = TRUE, ...) 
{
    type <- match.arg(type)
    if (type == "partial") {
        m <- match.call()
        m[[1]] <- as.name("pacf")
        m$type <- NULL
        return(eval(m, parent.frame()))
    }
    series <- deparse(substitute(x))
    x <- na.action(as.ts(x))
    x.freq <- frequency(x)
    x <- as.matrix(x)
    if (!is.numeric(x)) 
        stop("'x' must be numeric")
    sampleT <- nrow(x)
    nser <- ncol(x)
    if (is.null(lag.max)) 
        lag.max <- floor(10 * (log10(sampleT) - log10(nser)))
    lag.max <- min(lag.max, sampleT - 1)
    if (lag.max < 1) 
        stop("'lag.max' must be at least 1")
    if (demean) 
        x <- sweep(x, 2, colMeans(x, na.rm = TRUE))
    lag <- matrix(1, nser, nser)
    lag[lower.tri(lag)] <- -1
    acf <- array(.C("acf", as.double(x), as.integer(sampleT), 
        as.integer(nser), as.integer(lag.max), as.integer(type == 
            "correlation"), acf = double((lag.max + 1) * nser * 
            nser), NAOK = TRUE, PACKAGE = "stats")$acf, c(lag.max + 
        1, nser, nser))
    lag <- outer(0:lag.max, lag/x.freq)
    if (drop.lag.0) {
        if (type == "correlation") {
            acf = acf[-1, , , drop = FALSE]
            lag = lag[-1, , , drop = FALSE]
        }
    }
    acf.out <- structure(.Data = list(acf = acf, type = type, 
        n.used = sampleT, lag = lag, series = series, snames = colnames(x)), 
        class = "acf")
    if (plot) {
        plot1.acf(acf.out, ...)
        return(invisible(acf.out))
    }
    else return(acf.out)
}
