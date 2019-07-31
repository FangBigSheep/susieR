ridge <- function (X, y, sigma.b2, sigma2) {
    n <- nrow(X)
    p <- ncol(X)

    S.inv <- solve(crossprod(X) + (sigma2 / sigma.b2) * diag(p))
    beta.hat <- S.inv %*% t(X) %*% y
    beta.var <- sigma2 * crossprod(X %*% S.inv)

    result <- list()
    result$beta.hat <- beta.hat
    result$var      <- beta.var

    return (result)
}

eq <- function (log.sigma.b, XtX, Xty, sigma2) {
    sigma.b2 <- exp(2 * log.sigma.b)

    p <- ncol(XtX)
    S.inv <- solve(sigma.b2 * XtX + sigma2 * diag(p))

    u <- S.inv %*% Xty

    return (sum(u ^ 2) - p)
}

minus.llh <- function (log.sigma.b, XtX, Xty, sigma2) {
    # all irrelevent constants removed
    # and scaled by 2

    sigma.b2 <- exp(log.sigma.b * 2)
    p <- nrow(XtX)

    S <- XtX + (sigma2 / sigma.b2) * diag(p)

    term1 <- t(Xty) %*% solve(S) %*% Xty / sigma2
    term2 <- log(det(S)) - p * log(sigma.b2)

    return (- term1 - term2)
}

compute.llh <- function (sigma.b2, X, y, sigma2) {
    XtX <- crossprod(X)
    Xty <- crossprod(X, y)

    p <- nrow(XtX)

    S <- XtX + (sigma2 / sigma.b2) * diag(p)

    term1 <- t(Xty) %*% solve(S) %*% Xty / sigma2
    term2 <- log(det(S)) - p * log(sigma.b2)
    term3 <- - n * log (sigma2) - sum(y ^ 2) / sigma2

    return ((term1 + term2 + term3) / 2)
}

MLE <- function (X, y, sigma2, method = 'uniroot') {
    XtX <- crossprod(X)
    Xty <- crossprod(X, y)

    if (method == "uniroot")
        log.sigma.b <- uniroot(eq, interval = c(-6, 6),
                               XtX = XtX, Xty = Xty, sigma2 = sigma2)$root
    if (method == 'optim')
        log.sigma.b <- optim(par = 0, fn = minus.llh, XtX = XtX, Xty = Xty, sigma2 = sigma2,
                             method = 'Brent', upper = 4, lower = -4)$par

    return (exp(log.sigma.b * 2))
}

