`arimax` <-
function (x, order = c(0, 0, 0), seasonal = list(order = c(0, 
    0, 0), period = NA), xreg = NULL, include.mean = TRUE, transform.pars = TRUE, 
    fixed = NULL, init = NULL, method = c("CSS-ML", "ML", "CSS"), 
    n.cond, optim.control = list(), kappa = 1e+06, io = NULL, 
    xtransf, transfer = NULL) 
{
    if (!missing(transfer)) {
        for (i in seq(length(transfer))) transfer[[i]][2] = transfer[[i]][2] + 
            1
    }
    checkIO = function(x, phi = NULL, theta = NULL, Delta = NULL) {
        name1 = colnames(x)
        d = max(length(phi) + length(Delta), length(theta)) + 
            1
        ind = 0
        for (i in name1) {
            ind = ind + 1
            if (length(grep("IO-", i)) >= 1) 
                x[, ind] = armafilter(c(rep(0, d), x[, ind]), 
                  phi = phi, theta = c(1, theta), Delta = Delta)[-(1:d)]
        }
        x
    }
    armafilter = function(x, phi, theta, Delta) {
        if (length(theta) > 0) 
            x = filter(x, filter = theta, side = 1, method = "convolution")
        if (length(phi) > 0 && any(phi != 0)) 
            x = filter(x, filter = phi, side = 1, method = "recursive")
        if (length(Delta) > 0) 
            x = filter(x, filter = Delta, side = 1, method = "recursive")
        x
    }
    makePulse = function(point, xtsp) {
        if (!is.ts(x)) {
            pulse = 0 * x
            pulse[point] = 1
            return(pulse)
        }
        pulse = as.numeric(0 * x)
        if (length(point) == 1) {
            x[point] = 1
            return(x)
        }
        realp = as.integer(round((point[1] + (point[2] - 1)/xtsp[3] - 
            xtsp[1]) * xtsp[3] + 1))
        pulse[realp] = 1
        pulse
    }
    upx = function(x, par) {
        ncoef = narma + ncxreg
        ind = 0
        for (modorder in transfer) {
            ind = ind + 1
            p = modorder[1]
            q = modorder[2]
            if (p > 0) 
                phi = par[ncoef + (1:p)]
            else phi = NULL
            if (q > 0) 
                theta = par[ncoef + p + (1:q)]
            else theta = NULL
            ncoef = ncoef + p + q
            x = x - armafilter(xtransf[, ind], phi = phi, theta = theta, 
                Delta = NULL)
        }
        x
    }
    `%+%` <- function(a, b) .Call("TSconv", a, b, PACKAGE = "stats")
    upARIMA <- function(mod, phi, theta) {
        p <- length(phi)
        q <- length(theta)
        mod$phi <- phi
        mod$theta <- theta
        r <- max(p, q + 1)
        if (p > 0) 
            mod$T[1:p, 1] <- phi
        if (r > 1) 
            mod$Pn[1:r, 1:r] <- .Call("getQ0", phi, theta, PACKAGE = "stats")
        else if (p > 0) 
            mod$Pn[1, 1] <- 1/(1 - phi^2)
        else mod$Pn[1, 1] <- 1
        mod$a[] <- 0
        mod
    }
    arimaSS <- function(y, mod) {
        .Call("ARIMA_Like", y, mod$phi, mod$theta, mod$Delta, 
            mod$a, mod$P, mod$Pn, as.integer(0), TRUE, PACKAGE = "stats")
    }
    armafn <- function(p, trans) {
        par <- coef
        par[mask] <- p
        trarma <- .Call("ARIMA_transPars", par, arma, trans, 
            PACKAGE = "stats")
        Z <- upARIMA(mod, trarma[[1]], trarma[[2]])
        if (ncxreg > 0) {
            nxreg = checkIO(xreg, phi = Z$phi, theta = Z$theta, 
                Delta = Z$Delta)
            x <- x - nxreg %*% par[narma + (1:ncxreg)]
        }
        if (ntransf.term) 
            x = upx(x, par)
        res <- .Call("ARIMA_Like", x, Z$phi, Z$theta, Z$Delta, 
            Z$a, Z$P, Z$Pn, as.integer(0), FALSE, PACKAGE = "stats")
        s2 <- res[1]/res[3]
        0.5 * (log(s2) + res[2]/res[3])
    }
    armaCSS <- function(p) {
        par <- as.double(fixed)
        par[mask] <- p
        trarma <- .Call("ARIMA_transPars", par, arma, FALSE, 
            PACKAGE = "stats")
        if (ncxreg > 0) {
            nxreg = checkIO(xreg, phi = trarma[[1]], theta = trarma[[2]], 
                Delta = Delta)
            x <- x - nxreg %*% par[narma + (1:ncxreg)]
        }
        if (ntransf.term) 
            x = upx(x, par)
        res <- .Call("ARIMA_CSS", x, arma, trarma[[1]], trarma[[2]], 
            as.integer(ncond), FALSE, PACKAGE = "stats")
        0.5 * log(res)
    }
    arCheck <- function(ar) {
        p <- max(which(c(1, -ar) != 0)) - 1
        if (!p) 
            return(TRUE)
        all(Mod(polyroot(c(1, -ar[1:p]))) > 1)
    }
    maInvert <- function(ma) {
        q <- length(ma)
        q0 <- max(which(c(1, ma) != 0)) - 1
        if (!q0) 
            return(ma)
        roots <- polyroot(c(1, ma[1:q0]))
        ind <- Mod(roots) < 1
        if (all(!ind)) 
            return(ma)
        if (q0 == 1) 
            return(c(1/ma[1], rep(0, q - q0)))
        roots[ind] <- 1/roots[ind]
        x <- 1
        for (r in roots) x <- c(x, 0) - c(0, x)/r
        c(Re(x[-1]), rep(0, q - q0))
    }
    series <- deparse(substitute(x))
    if (NCOL(x) > 1) 
        stop("only implemented for univariate time series")
    method <- match.arg(method)
    x <- as.ts(x)
    if (!is.numeric(x)) 
        stop("'x' must be numeric")
    storage.mode(x) <- "double"
    dim(x) <- NULL
    n <- length(x)
    if (!missing(order)) 
        if (!is.numeric(order) || length(order) != 3 || any(order < 
            0)) 
            stop("'order' must be a non-negative numeric vector of length 3")
    if (!missing(seasonal)) 
        if (is.list(seasonal)) {
            if (is.null(seasonal$order)) 
                stop("'seasonal' must be a list with component 'order'")
            if (!is.numeric(seasonal$order) || length(seasonal$order) != 
                3 || any(seasonal$order < 0)) 
                stop("'seasonal$order' must be a non-negative numeric vector of length 3")
        }
        else if (is.numeric(order)) {
            if (length(order) == 3) 
                seasonal <- list(order = seasonal)
            else ("'seasonal' is of the wrong length")
        }
        else stop("'seasonal' must be a list with component 'order'")
    if (is.null(seasonal$period) || is.na(seasonal$period) || 
        seasonal$period == 0) 
        seasonal$period <- frequency(x)
    arma <- as.integer(c(order[-2], seasonal$order[-2], seasonal$period, 
        order[2], seasonal$order[2]))
    narma <- sum(arma[1:4])
    xtsp <- tsp(x)
    tsp(x) <- NULL
    Pm = NULL
    if (!missing(io)) {
        Pm = sapply(io, makePulse, xtsp)
        name1 = substitute(io)
        namev = NULL
        for (i in 2:length(name1)) namev = c(namev, paste("IO", 
            deparse(name1[[i]]), sep = "-"))
        colnames(Pm) = namev
    }
    Delta <- 1
    for (i in seq(len = order[2])) Delta <- Delta %+% c(1, -1)
    for (i in seq(len = seasonal$order[2])) Delta <- Delta %+% 
        c(1, rep(0, seasonal$period - 1), -1)
    Delta <- -Delta[-1]
    nd <- order[2] + seasonal$order[2]
    n.used <- sum(!is.na(x)) - length(Delta)
    if(!is.null(xreg)) xreg=data.frame(xreg)
    if (!is.null(xreg)) 
        xreg = as.matrix(xreg)
    xreg = cbind(xreg, Pm)
    if (is.null(xreg)) {
        ncxreg <- 0
    }
    else {
        nmxreg <- deparse(substitute(xreg))
        if (NROW(xreg) != n) 
            stop("lengths of 'x' and 'xreg' do not match")
        ncxreg <- NCOL(xreg)
        xreg <- as.matrix(xreg)
        storage.mode(xreg) <- "double"
    }
    class(xreg) <- NULL
    if (ncxreg > 0 && is.null(colnames(xreg))) 
        colnames(xreg) <- if (ncxreg == 1) 
            nmxreg
        else paste(nmxreg, 1:ncxreg, sep = "")
    if (include.mean && (nd == 0)) {
        xreg <- cbind(intercept = rep(1, n), xreg = xreg)
        ncxreg <- ncxreg + 1
    }
    ntransf = 0
    ntransf.term = 0
    ct = NULL
    if (!is.null(transfer)) {
        xtransf = as.matrix(xtransf)
        storage.mode(xtransf) = "double"
        ntransf.term = length(transfer)
        ntransf = sum(sapply(transfer, sum))
        if (missing(xtransf) && length(transfer) != ncol(xtransf)) 
            stop("no of columns of xtransf need to equal the number of transfer sub-models")
        ncoef = narma + ncxreg
        ind = 0
        if (!is.null(colnames(xtransf))) 
            tnames = colnames(xtransf)
        else tnames = paste("T", 1:ntransf.term, sep = "")
        for (modorder in transfer) {
            ind = ind + 1
            p = modorder[1]
            q = modorder[2]
            if (p > 0) 
                ct = c(ct, paste(paste(tnames[ind], "AR", sep = "-"), 
                  1:p, sep = ""))
            if (q > 0) 
                ct = c(ct, paste(paste(tnames[ind], "MA", sep = "-"), 
                  1:q - 1, sep = ""))
        }
    }
    if (method == "CSS-ML") {
        anyna <- any(is.na(x))
        if (ncxreg) 
            anyna <- anyna || any(is.na(xreg))
        if (anyna) 
            method <- "ML"
    }
    if (method == "CSS" || method == "CSS-ML") {
        ncond <- order[2] + seasonal$order[2] * seasonal$period
        ncond1 <- order[1] + seasonal$period * seasonal$order[1]
        ncond <- if (!missing(n.cond)) 
            ncond + max(n.cond, ncond1)
        else ncond + ncond1
    }
    else ncond <- 0
    if (is.null(fixed)) 
        fixed <- rep(as.numeric(NA), narma + ncxreg + ntransf)
    else if (length(fixed) != narma + ncxreg + ntransf) 
        stop("wrong length for 'fixed'")
    mask <- is.na(fixed)
    no.optim <- !any(mask)
    if (no.optim) 
        transform.pars <- FALSE
    if (transform.pars) {
        ind <- arma[1] + arma[2] + seq(length = arma[3])
        if (any(!mask[seq(length = arma[1])]) || any(!mask[ind])) {
            warning("some AR parameters were fixed: setting transform.pars = FALSE")
            transform.pars <- FALSE
        }
    }
    init0 <- rep(0, narma)
    parscale <- rep(1, narma)
    if (ncxreg) {
        cn <- colnames(xreg)
        orig.xreg <- (ncxreg == 1) || any(!mask[narma + 1:ncxreg])
        if (!orig.xreg) {
            S <- svd(na.omit(xreg))
            xreg <- xreg %*% S$v
        }
        if (include.mean && (nd == 0)) {
            fit <- lm(x ~ xreg - 1, na.action = na.omit)
            n.used <- sum(!is.na(resid(fit))) - length(Delta)
            init0 <- c(init0, coef(fit))
            ses <- summary(fit)$coefficients[, 2]
        }
        else {
            fit <- lm(x ~ xreg, na.action = na.omit)
            n.used <- sum(!is.na(resid(fit))) - length(Delta)
            init0 <- c(init0, coef(fit)[-1])
            ses <- summary(fit)$coef[-1, 2]
        }
        parscale <- c(parscale, 10 * ses)
    }
    if (ntransf) {
        init0 = c(init0, rep(0, ntransf))
    }
    if (ntransf.term) {
        for (i in 1:ntransf.term) {
            fit = lm(x ~ xtransf[, i], na.action = na.omit)
            ses = summary(fit)$coef[2, 2]
            parscale = c(parscale, rep(10 * ses, sum(transfer[[i]])))
        }
    }
    if (n.used <= 0) 
        stop("too few non-missing observations")
    if (!is.null(init)) {
        if (length(init) != length(init0)) 
            stop("'init' is of the wrong length")
        if (any(ind <- is.na(init))) 
            init[ind] <- init0[ind]
        if (method == "ML") {
            if (arma[1] > 0) 
                if (!arCheck(init[1:arma[1]])) 
                  stop("non-stationary AR part")
            if (arma[3] > 0) 
                if (!arCheck(init[sum(arma[1:2]) + 1:arma[3]])) 
                  stop("non-stationary seasonal AR part")
            if (transform.pars) 
                init <- .Call("ARIMA_Invtrans", as.double(init), 
                  arma, PACKAGE = "stats")
        }
    }
    else init <- init0
    coef <- as.double(fixed)
    if (!("parscale" %in% names(optim.control))) 
        optim.control$parscale <- parscale[mask]
    if (method == "CSS") {
        res <- if (no.optim) 
            list(convergence = 0, par = numeric(0), value = armaCSS(numeric(0)))
        else optim(init[mask], armaCSS, method = "BFGS", hessian = TRUE, 
            control = optim.control)
        if (res$convergence > 0) 
            warning("possible convergence problem: optim gave code=", 
                res$convergence)
        coef[mask] <- res$par
        trarma <- .Call("ARIMA_transPars", coef, arma, FALSE, 
            PACKAGE = "stats")
        mod <- makeARIMA(trarma[[1]], trarma[[2]], Delta, kappa)
        if (ncxreg > 0) {
            nxreg = checkIO(xreg, phi = trarma[[1]], theta = trarma[[2]], 
                Delta = Delta)
            x <- x - nxreg %*% coef[narma + (1:ncxreg)]
        }
        if (ntransf.term) 
            x = upx(x, coef)
        arimaSS(x, mod)
        val <- .Call("ARIMA_CSS", x, arma, trarma[[1]], trarma[[2]], 
            as.integer(ncond), TRUE, PACKAGE = "stats")
        sigma2 <- val[[1]]
        var <- if (no.optim) 
            numeric(0)
        else solve(res$hessian * n.used)
    }
    else {
        if (method == "CSS-ML") {
            res <- if (no.optim) 
                list(convergence = 0, par = numeric(0), value = armaCSS(numeric(0)))
            else optim(init[mask], armaCSS, method = "BFGS", 
                hessian = FALSE, control = optim.control)
            if (res$convergence == 0) 
                init[mask] <- res$par
            if (arma[1] > 0) 
                if (!arCheck(init[1:arma[1]])) 
                  stop("non-stationary AR part from CSS")
            if (arma[3] > 0) 
                if (!arCheck(init[sum(arma[1:2]) + 1:arma[3]])) 
                  stop("non-stationary seasonal AR part from CSS")
            ncond <- 0
        }
        if (transform.pars) {
            init <- .Call("ARIMA_Invtrans", init, arma, PACKAGE = "stats")
            if (arma[2] > 0) {
                ind <- arma[1] + 1:arma[2]
                init[ind] <- maInvert(init[ind])
            }
            if (arma[4] > 0) {
                ind <- sum(arma[1:3]) + 1:arma[4]
                init[ind] <- maInvert(init[ind])
            }
        }
        trarma <- .Call("ARIMA_transPars", init, arma, transform.pars, 
            PACKAGE = "stats")
        mod <- makeARIMA(trarma[[1]], trarma[[2]], Delta, kappa)
        res <- if (no.optim) 
            list(convergence = 0, par = numeric(0), value = armafn(numeric(0), 
                as.logical(transform.pars)))
        else optim(init[mask], armafn, method = "BFGS", hessian = TRUE, 
            control = optim.control, trans = as.logical(transform.pars))
        if (res$convergence > 0) 
            warning("possible convergence problem: optim gave code=", 
                res$convergence)
        coef[mask] <- res$par
        if (transform.pars) {
            if (arma[2] > 0) {
                ind <- arma[1] + 1:arma[2]
                if (all(mask[ind])) 
                  coef[ind] <- maInvert(coef[ind])
            }
            if (arma[4] > 0) {
                ind <- sum(arma[1:3]) + 1:arma[4]
                if (all(mask[ind])) 
                  coef[ind] <- maInvert(coef[ind])
            }
            if (any(coef[mask] != res$par)) {
                oldcode <- res$convergence
                res <- optim(coef[mask], armafn, method = "BFGS", 
                  hessian = TRUE, control = list(maxit = 0, parscale = optim.control$parscale), 
                  trans = TRUE)
                res$convergence <- oldcode
                coef[mask] <- res$par
            }
            A <- .Call("ARIMA_Gradtrans", as.double(coef), arma, 
                PACKAGE = "stats")
            A <- A[mask, mask]
            var <- t(A) %*% solve(res$hessian * n.used) %*% A
            coef <- .Call("ARIMA_undoPars", coef, arma, PACKAGE = "stats")
        }
        else var <- if (no.optim) 
            numeric(0)
        else solve(res$hessian * n.used)
        trarma <- .Call("ARIMA_transPars", coef, arma, FALSE, 
            PACKAGE = "stats")
        mod <- makeARIMA(trarma[[1]], trarma[[2]], Delta, kappa)
        if (ncxreg > 0) {
            nxreg = checkIO(xreg, phi = trarma[[1]], theta = trarma[[2]], 
                Delta = Delta)
            x = x - nxreg %*% coef[narma + (1:ncxreg)]
        }
        if (ntransf.term) 
            x = upx(x, coef)
        val = arimaSS(x, mod)
        sigma2 <- val[[1]][1]/n.used
    }
    value <- 2 * n.used * res$value + n.used + n.used * log(2 * 
        pi)
    aic <- if (method != "CSS") 
        value + 2 * sum(mask)
    else NA
    nm <- NULL
    if (arma[1] > 0) 
        nm <- c(nm, paste("ar", 1:arma[1], sep = ""))
    if (arma[2] > 0) 
        nm <- c(nm, paste("ma", 1:arma[2], sep = ""))
    if (arma[3] > 0) 
        nm <- c(nm, paste("sar", 1:arma[3], sep = ""))
    if (arma[4] > 0) 
        nm <- c(nm, paste("sma", 1:arma[4], sep = ""))
    if (ncxreg > 0) {
        nm <- c(nm, cn)
        if (!orig.xreg) {
            ind <- narma + 1:ncxreg
            coef[ind] <- S$v %*% coef[ind]
            A <- diag(length(coef))
            A[ind, ind] <- S$v
            A <- A[mask, mask]
            var <- A %*% var %*% t(A)
        }
    }
    nm = c(nm, ct)
    names(coef) <- nm
    if (!no.optim) 
        dimnames(var) <- list(nm[mask], nm[mask])
    resid <- val[[2]]
    tsp(resid) <- xtsp
    class(resid) <- "ts"
    res <- list(coef = coef, sigma2 = sigma2, var.coef = var, 
        mask = mask, loglik = -0.5 * value, aic = aic, arma = arma, 
        residuals = resid, call = match.call(), series = series, 
        code = res$convergence, n.cond = ncond, model = mod)
    class(res) <- c("Arimax", "Arima")
    res
}

