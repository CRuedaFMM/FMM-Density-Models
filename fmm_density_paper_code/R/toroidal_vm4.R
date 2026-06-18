###############################################################
## VM4: mixture of bivariate product von Mises components
###############################################################

vm_norm <- function(kappa) 1 / (2*pi*besselI(kappa, 0, expon.scaled = FALSE))

dvm_torus <- function(x, mu, kappa) vm_norm(kappa) * exp(kappa * cos(x - mu))

kappa_from_Rbar <- function(Rbar) {
  Rbar <- min(max(Rbar, 1e-8), 0.999999)
  if (Rbar < 0.53) return(2*Rbar + Rbar^3 + (5*Rbar^5)/6)
  if (Rbar < 0.85) return(-0.4 + 1.39*Rbar + 0.43/(1 - Rbar))
  1/(Rbar^3 - 4*Rbar^2 + 3*Rbar)
}

circ_mean <- function(x, w) atan2(sum(w*sin(x)), sum(w*cos(x))) %% (2*pi)

circ_Rbar <- function(x, w) sqrt(sum(w*cos(x))^2 + sum(w*sin(x))^2) / sum(w)

fit_vm_prod_mix <- function(X, Mi = 4, maxit = 200, tol = 1e-6, seed = 1) {
  set.seed(seed)
  X <- as.matrix(X)
  n <- nrow(X); th <- X[,1]; ph <- X[,2]
  z <- sample.int(Mi, n, replace = TRUE)
  pi <- tabulate(z, Mi) / n
  mu1 <- as.numeric(tapply(th, z, mean)); mu2 <- as.numeric(tapply(ph, z, mean))
  mu1[is.na(mu1)] <- runif(sum(is.na(mu1)), 0, 2*pi)
  mu2[is.na(mu2)] <- runif(sum(is.na(mu2)), 0, 2*pi)
  k1 <- rep(5, Mi); k2 <- rep(5, Mi)
  ll_old <- -Inf
  for (it in seq_len(maxit)) {
    logdens <- matrix(NA_real_, n, Mi)
    for (j in seq_len(Mi)) {
      logdens[,j] <- log(pi[j] + 1e-300) +
        log(vm_norm(k1[j]) + 1e-300) + k1[j]*cos(th - mu1[j]) +
        log(vm_norm(k2[j]) + 1e-300) + k2[j]*cos(ph - mu2[j])
    }
    mmax <- apply(logdens, 1, max)
    w <- exp(logdens - mmax); w <- w / rowSums(w)
    Nj <- colSums(w); pi <- Nj / n
    for (j in seq_len(Mi)) {
      mu1[j] <- circ_mean(th, w[,j]); mu2[j] <- circ_mean(ph, w[,j])
      k1[j] <- kappa_from_Rbar(circ_Rbar(th, w[,j]))
      k2[j] <- kappa_from_Rbar(circ_Rbar(ph, w[,j]))
    }
    ll <- sum(mmax + log(rowSums(exp(logdens - mmax))))
    if (abs(ll - ll_old) < tol) break
    ll_old <- ll
  }
  list(pmix = pi, mu1 = mu1, mu2 = mu2, kappa1 = k1, kappa2 = k2, ll = ll_old)
}

eval_vm_prod_mix_grid <- function(par, theta_grid, phi_grid) {
  Q <- matrix(0, length(theta_grid), length(phi_grid))
  for (j in seq_along(par$pmix)) {
    a <- dvm_torus(theta_grid, par$mu1[j], par$kappa1[j])
    b <- dvm_torus(phi_grid, par$mu2[j], par$kappa2[j])
    Q <- Q + par$pmix[j] * (a %o% b)
  }
  normalize_grid(Q)
}

fit_mix_L2_from_FF <- function(FF, Mi = 4, N = 3000, seed = 1) {
  FF <- as.matrix(FF); storage.mode(FF) <- "double"
  grids <- make_grid(nrow(FF), ncol(FF))
  theta_grid <- grids$theta; phi_grid <- grids$phi
  P <- normalize_grid(FF)
  X <- sample_from_grid(P, theta_grid, phi_grid, N = N, seed = seed)
  par <- fit_vm_prod_mix(X, Mi = Mi, seed = seed)
  Q <- eval_vm_prod_mix_grid(par, theta_grid, phi_grid)
  L2sq <- L2_sq_grid(P, Q)
  list(L2 = sqrt(L2sq), L2sq = L2sq, par = par, P = P, Q = Q)
}
