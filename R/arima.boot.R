`arima.boot` <-
function(arima.fit,cond.boot=FALSE,is.normal=TRUE,B=1000,init,ntrans=100)
{
#
# Programmed by Kung-Sik Chan
# Date: April 19, 2006
#
# This function bootstraps time series according to the 
# fitted ARMA(p,d,q) model supplied by the fitted object arima.fit.
#
# Input: arima.fit=a fitted object from arima function. (Currently, it excludes
#          SARIMA)
#        cond.boot=does bootstrap conditional on the (p+d) initial values), if
#          it is set true. If false (default), stationary bootstrap is used. 
#        is.normal=if true (default), errors are normally distributed, otherwise
#          errors are drawn randomly and with replacement from the residuals of the
#          fitted model.
#        B=no of bootstrap replicates (1000, default)
#        init=initial values for the bootstrap; needed if cond.boot=True
#          default values are the initial values of the time series of the fitted
#          model.
#        ntrans=number of tranient values for the stationary bootstrap. 
#          Default=100
#                            
# Output: a matrix each row of which consists of the coefficient estimates of 
#         a bootstrap time-series.
#
# Examples: 
#
# arima.hare=arima(sqrt(hare2),order=c(3,0,0))
# boot.hare=arima.boot(arima.hare,B=50,init=sqrt(hare2)[1:3],ntrans=100)
# apply(boot.hare,2,quantile, c(.025,.975))
#
# bootstrap the quasi-period estimates.
#period.boot=apply(boot.hare,1,function(x){
#roots=polyroot(c(1,-x[1:3]))
#min1=1.e+9
#rootc=NA
#for (root in roots) {
#if( abs(Im(root))<1e-10) next
#if (Mod(root)< min1) {min1=Mod(root); rootc=root}
#}
#if(is.na(rootc)) period=NA else period=2*pi/abs(Arg(rootc))
#period
#})
#hist(period.boot)
#quantile(period.boot,c(0.025,.975))
#

`arimab.sim` <-
function(model=list(ar=NA,ma=NA,d=0),n,ntrans=100,init=rep(0,100),cond.boot=FALSE,
residuals,include.mean=TRUE,mean,is.normal=TRUE,sd,...){
#
# workhorse of the arima.boot function
#
if(!is.null(model$ar)) ar=model$ar else ar=NA
if(!is.null(model$ma)) ma=model$ma else ma=NA
if(!is.null(model$d)) d=model$d else d=0
if(include.mean) intercept=mean else intercept=0
if(!any(is.na(ar))) 
{d1=d
if(d1>0) {
bar=c(1,-ar)
while(d1>0) {bar=c(bar,0)-c(0,bar);d1=d1-1}
ar=-bar[-1]}
p=length(ar)} else {
p=0
d1=d
if(d1>0){ar=c(1,0)
while(d1>0) {bar=c(bar,0)-c(0,bar);d1=d1-1}
ar=-bar[-1]}
if(!any(is.na(ar))) p=length(ar)
}
if(!is.na(ma)) q=length(ma) else q=0
initial=rev(init[1:p])-intercept
if(cond.boot) {ntrans=0}
if(is.normal){
if(q>0) noise=filter(rnorm(n=n+ntrans,mean=0,sd=sd), 
init=rnorm(n=q,mean=0,sd=sd),filter=ma,method='convolution',side=1) else
noise=rnorm(n=n+ntrans,mean=0,sd=sd)
} else {
if(q>0) 
noise=filter(sample(residuals,replace=TRUE,size=n+ntrans),
init=sample(residuals,size=q,replace=TRUE),
filter=ma,method='convolution',side=1) else
noise=sample(residuals,size=n+ntrans,replace=TRUE)
}
boot=filter(noise,filter=ar,method='recursive',
init=initial,side=1)+intercept
boot=boot[(ntrans+1):(n+ntrans)]
if(cond.boot) boot=c(rev(initial)+intercept,rev(rev(boot)[-seq(initial)]))
#browser()
boot
}



`arimab` <-
function (x, order = c(0, 0, 0), seasonal = list(order = c(0, 
    0, 0), period = NA), xreg = NULL, include.mean = TRUE, transform.pars = TRUE, 
    fixed = NULL, init = NULL, method = c("CSS-ML", "ML", "CSS"), 
    n.cond, optim.control = list(), kappa = 1e+06) 
{
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
        if (ncxreg > 0) 
            x <- x - xreg %*% par[narma + (1:ncxreg)]
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
        if (ncxreg > 0) 
            x <- x - xreg %*% par[narma + (1:ncxreg)]
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
    Delta <- 1
    for (i in seq(len = order[2])) Delta <- Delta %+% c(1, -1)
    for (i in seq(len = seasonal$order[2])) Delta <- Delta %+% 
        c(1, rep(0, seasonal$period - 1), -1)
    Delta <- -Delta[-1]
    nd <- order[2] + seasonal$order[2]
    n.used <- sum(!is.na(x)) - length(Delta)
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
        fixed <- rep(as.numeric(NA), narma + ncxreg)
    else if (length(fixed) != narma + ncxreg) 
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
        fit <- lm(x ~ xreg - 1, na.action = na.omit)
        n.used <- sum(!is.na(resid(fit))) - length(Delta)
        init0 <- c(init0, coef(fit))
        ses <- summary(fit)$coef[, 2]
        parscale <- c(parscale, 10 * ses)
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
                  return(NA)
            if (arma[3] > 0) 
                if (!arCheck(init[sum(arma[1:2]) + 1:arma[3]])) 
                  return(NA)
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
        if (ncxreg > 0) 
            x <- x - xreg %*% coef[narma + (1:ncxreg)]
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
                  return(NA)
            if (arma[3] > 0) 
                if (!arCheck(init[sum(arma[1:2]) + 1:arma[3]])) 
                  return(NA)
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
        val <- if (ncxreg > 0) 
            arimaSS(x - xreg %*% coef[narma + (1:ncxreg)], mod)
        else arimaSS(x, mod)
        sigma2 <- val[[1]][1]/n.used
    }
    value <- 2 * n.used * res$value + n.used + n.used * log(2 * 
        pi)
    aic <- if (method != "CSS") 
        value + 2 * sum(mask) + 2
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
            A <- diag(narma + ncxreg)
            A[ind, ind] <- S$v
            A <- A[mask, mask]
            var <- A %*% var %*% t(A)
        }
    }
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
    class(res) <- "Arima"
    res
}


        
if(!missing(arima.fit)){
order=arima.fit$call$order
p=order[[2]]
d=order[[3]]
q=order[[4]]
order=c(p,d,q)
n=length(eval(arima.fit$call[[2]]))
if(p>0) ar=coef(arima.fit)[1:p] else ar=NA
if(q>0) ma=coef(arima.fit)[(p+1):(p+q)] else ma=NA
model=list(ar=ar,ma=ma,d=d)
if(any(is.null(arima.fit$call$include.mean))) include.mean=TRUE else include.mean=(arima.fit$call$include.mean)=="T"
if(include.mean) mean=coef(arima.fit)[(p+q+1)] else mean=0
sd=arima.fit$sigma2^.5
residuals=residuals(arima.fit)
if(missing(init)) {
init=eval(arima.fit$call[[2]])[1:(p+d)]
}
} else stop("Need to input the fitted model")
coefm=NULL
for (i in 1:B){
boot=arimab.sim(model=model,n=n,ntrans=ntrans,init=init,cond.boot=cond.boot,
residuals=residuals,include.mean=include.mean,mean=mean,is.normal=is.normal,
sd=sd)
res=arimab(boot,order=order,include.mean=include.mean)
if(!any(is.na(res))) coefm=rbind(coefm,c(coef(res),res$sigma2))
}
colnames(coefm)=c(names(coef(res)),"noise var")
invisible(coefm)
}

