############################################################
## Product von Mises mixtures on the torus (VM4 benchmark)
############################################################

vm_norm <- function(kappa) 1/(2*pi*besselI(kappa, 0, expon.scaled = FALSE))

dvm <- function(x, mu, kappa){
  vm_norm(kappa) * exp(kappa*cos(x - mu))
}

kappa_from_Rbar <- function(Rbar){
  if(Rbar < 1e-8) return(0)
  if(Rbar < 0.53) return(2*Rbar + Rbar^3 + (5*Rbar^5)/6)
  if(Rbar < 0.85) return(-0.4 + 1.39*Rbar + 0.43/(1 - Rbar))
  1/(Rbar^3 - 4*Rbar^2 + 3*Rbar)
}

circ_mean <- function(x, w){
  atan2(sum(w*sin(x)), sum(w*cos(x))) %% (2*pi)
}

circ_Rbar <- function(x, w){
  C <- sum(w*cos(x)); S <- sum(w*sin(x))
  sqrt(C^2 + S^2)/sum(w)
}

sample_from_grid <- function(P, theta_grid, phi_grid, N = 3000, seed = 1){
  set.seed(seed)
  P <- as.matrix(P); nr <- nrow(P); nc <- ncol(P)
  P[!is.finite(P)] <- 0
  s <- sum(P)
  if(!is.finite(s) || s <= 0) stop("P is not a valid probability grid.")
  P <- P/s
  idx <- sample.int(nr*nc, size = N, replace = TRUE, prob = as.vector(P))
  i <- ((idx - 1) %% nr) + 1
  j <- ((idx - 1) %/% nr) + 1
  X <- cbind(theta_grid[i], phi_grid[j])
  colnames(X) <- c("theta", "phi")
  X
}

fit_vm_prod_mix <- function(X, M = 4, maxit = 200, tol = 1e-6, seed = 1){
  set.seed(seed)
  X <- as.matrix(X)
  n <- nrow(X)
  th <- wrap_2pi(X[,1]); ph <- wrap_2pi(X[,2])

  z <- sample.int(M, n, replace = TRUE)
  pmix <- tabulate(z, M)/n

  mu1 <- as.numeric(tapply(th, z, mean))
  mu2 <- as.numeric(tapply(ph, z, mean))
  mu1[is.na(mu1)] <- runif(sum(is.na(mu1)), 0, 2*pi)
  mu2[is.na(mu2)] <- runif(sum(is.na(mu2)), 0, 2*pi)

  k1 <- rep(5, M); k2 <- rep(5, M)
  ll_old <- -Inf

  for(it in seq_len(maxit)){
    logdens <- matrix(NA_real_, n, M)
    for(j in 1:M){
      logdens[,j] <- log(pmix[j] + 1e-300) +
        log(vm_norm(k1[j]) + 1e-300) + k1[j]*cos(th - mu1[j]) +
        log(vm_norm(k2[j]) + 1e-300) + k2[j]*cos(ph - mu2[j])
    }

    mmax <- apply(logdens, 1, max)
    w <- exp(logdens - mmax)
    w <- w/rowSums(w)

    Nj <- colSums(w)
    pmix <- Nj/n

    for(j in 1:M){
      if(Nj[j] <= .Machine$double.eps){
        mu1[j] <- runif(1, 0, 2*pi); mu2[j] <- runif(1, 0, 2*pi)
        k1[j] <- 0.5; k2[j] <- 0.5
      } else {
        mu1[j] <- circ_mean(th, w[,j])
        mu2[j] <- circ_mean(ph, w[,j])
        k1[j] <- kappa_from_Rbar(pmin(pmax(circ_Rbar(th, w[,j]), 1e-8), 0.999999))
        k2[j] <- kappa_from_Rbar(pmin(pmax(circ_Rbar(ph, w[,j]), 1e-8), 0.999999))
      }
    }

    ll <- sum(mmax + log(rowSums(exp(logdens - mmax))))
    if(abs(ll - ll_old) < tol) break
    ll_old <- ll
  }

  list(pmix = pmix, mu1 = mu1, mu2 = mu2, kappa1 = k1, kappa2 = k2, ll = ll_old)
}

eval_vm_prod_mix_grid <- function(par, theta_grid, phi_grid){
  n_theta <- length(theta_grid); n_phi <- length(phi_grid)
  Q <- matrix(0, n_theta, n_phi)
  for(j in seq_along(par$pmix)){
    a <- dvm(theta_grid, par$mu1[j], par$kappa1[j])
    b <- dvm(phi_grid, par$mu2[j], par$kappa2[j])
    Q <- Q + par$pmix[j]*(a %o% b)
  }
  normalize_grid(Q)
}

loglik_vm_prod_mix <- function(X, par){
  X <- as.matrix(X)
  th <- wrap_2pi(X[,1]); ph <- wrap_2pi(X[,2])
  M <- length(par$pmix)
  ll <- 0
  for(i in seq_len(nrow(X))){
    v <- numeric(M)
    for(j in 1:M){
      v[j] <- log(par$pmix[j] + 1e-300) +
        log(vm_norm(par$kappa1[j]) + 1e-300) + par$kappa1[j]*cos(th[i] - par$mu1[j]) +
        log(vm_norm(par$kappa2[j]) + 1e-300) + par$kappa2[j]*cos(ph[i] - par$mu2[j])
    }
    ll <- ll + logsumexp(v)
  }
  ll
}

fit_mix_L2_from_FF <- function(FF, Mi = 4, N = 3000, seed = 1){
  FF <- as.matrix(FF)
  grids <- make_grid(nrow(FF), ncol(FF))
  theta_grid <- grids$theta
  phi_grid <- grids$phi
  P <- normalize_grid(FF)
  X <- sample_from_grid(P, theta_grid, phi_grid, N = N, seed = seed)
  par <- fit_vm_prod_mix(X, M = Mi, seed = seed)
  Q <- eval_vm_prod_mix_grid(par, theta_grid, phi_grid)
  L2sq <- L2_sq_grid(P, Q)
  list(L2 = sqrt(L2sq), L2sq = L2sq, par = par, P = P, Q = Q)
}
