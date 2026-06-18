###############################################################
## Circular density models for aFMM methodology
###############################################################

wrap2pi <- function(x) x %% (2*pi)

safe_density <- function(x, eps = 1e-12) pmax(x, eps)

Phi_aFMM <- function(theta, alpha, omega) {
  u <- theta - alpha
  2 * atan2(omega * sin(u / 2), cos(u / 2))
}

d_aFMM <- function(theta, alpha, omega, beta, r, log = FALSE) {
  theta <- wrap2pi(theta)
  Phi <- Phi_aFMM(theta, alpha, omega)
  c1 <- (1 - omega) / (1 + omega)
  Z <- 2*pi * (1 + r * c1 * cos(beta))
  dens <- (1 + r * cos(Phi + beta)) / Z
  dens <- safe_density(dens)
  if (log) return(log(dens))
  dens
}

dvm <- function(theta, mu, kappa) {
  exp(kappa * cos(theta - mu)) / (2*pi*besselI(kappa, 0))
}

dwc <- function(theta, mu, rho) {
  (1 - rho^2) / (2*pi * (1 + rho^2 - 2*rho*cos(theta - mu)))
}

dmix_vm2 <- function(theta, mu1, k1, mu2, k2, pi1) {
  pi1 * dvm(theta, mu1, k1) + (1 - pi1) * dvm(theta, mu2, k2)
}

dmix_wc2 <- function(theta, mu1, rho1, mu2, rho2, pi1) {
  pi1 * dwc(theta, mu1, rho1) + (1 - pi1) * dwc(theta, mu2, rho2)
}

dmix_afmm2 <- function(theta,
                       alpha1, omega1, r1, beta1,
                       alpha2, omega2, r2, beta2,
                       pi1) {
  pi1 * d_aFMM(theta, alpha1, omega1, beta1, r1) +
    (1 - pi1) * d_aFMM(theta, alpha2, omega2, beta2, r2)
}

dfourier2 <- function(theta, a1, b1, a2, b2) {
  d <- (1 + a1*cos(theta) + b1*sin(theta) +
          a2*cos(2*theta) + b2*sin(2*theta)) / (2*pi)
  safe_density(d)
}

d_cardioid_model <- function(theta, mu, rho) {
  safe_density((1 + 2*rho*cos(theta - mu)) / (2*pi))
}

d_vonmises_model <- function(theta, mu, kappa) dvm(theta, mu, kappa)

d_wrapped_cauchy_model <- function(theta, mu, gamma) dwc(theta, mu, gamma)

d_sineskew_vm_model <- function(theta, mu, kappa, lambda) {
  safe_density(d_vonmises_model(theta, mu, kappa) * (1 + lambda * sin(theta - mu)))
}

d_extended_wc_model <- function(theta, mu, gamma, lambda) {
  safe_density(d_wrapped_cauchy_model(theta, mu, gamma) * (1 + lambda * sin(theta - mu)))
}

d_gc2_model <- function(theta, mu, a1, b1, a2, b2) {
  z <- 1 + a1*cos(theta - mu) + b1*sin(theta - mu) +
    a2*cos(2*(theta - mu)) + b2*sin(2*(theta - mu))
  z <- safe_density(z)
  Z <- integrate(function(t) {
    safe_density(1 + a1*cos(t - mu) + b1*sin(t - mu) +
                   a2*cos(2*(t - mu)) + b2*sin(2*(t - mu)))
  }, lower = 0, upper = 2*pi)$value
  z / Z
}

nll_aFMM <- function(par, theta) {
  alpha <- wrap2pi(par[1]); omega <- plogis(par[2])
  beta  <- wrap2pi(par[3]); r <- plogis(par[4])
  -sum(log(safe_density(d_aFMM(theta, alpha, omega, beta, r))))
}

nll_vm <- function(par, theta) {
  mu <- wrap2pi(par[1]); kappa <- exp(par[2])
  -sum(log(safe_density(dvm(theta, mu, kappa))))
}

nll_wc <- function(par, theta) {
  mu <- wrap2pi(par[1]); rho <- plogis(par[2])
  -sum(log(safe_density(dwc(theta, mu, rho))))
}

nll_mix_vm2 <- function(par, theta) {
  d <- dmix_vm2(theta, wrap2pi(par[1]), exp(par[3]),
                wrap2pi(par[2]), exp(par[4]), plogis(par[5]))
  if (any(!is.finite(d)) || any(d <= 0)) return(1e10)
  -sum(log(d))
}

nll_mix_wc2 <- function(par, theta) {
  d <- dmix_wc2(theta, wrap2pi(par[1]), plogis(par[3]),
                wrap2pi(par[2]), plogis(par[4]), plogis(par[5]))
  if (any(!is.finite(d)) || any(d <= 0)) return(1e10)
  -sum(log(d))
}

nll_mixture_aFMM2 <- function(par, theta) {
  d <- dmix_afmm2(
    theta,
    wrap2pi(par[1]), plogis(par[3]), plogis(par[5]), wrap2pi(par[7]),
    wrap2pi(par[2]), plogis(par[4]), plogis(par[6]), wrap2pi(par[8]),
    plogis(par[9])
  )
  if (any(!is.finite(d)) || any(d <= 0)) return(1e10)
  -sum(log(d))
}

nll_fourier2 <- function(par, theta) {
  d <- dfourier2(theta, par[1], par[2], par[3], par[4])
  if (any(!is.finite(d)) || any(d <= 0)) return(1e10)
  -sum(log(d))
}

density_from_fit <- function(fit, grid) {
  p <- fit$raw_par
  if (fit$model == "mix_vM2") {
    return(dmix_vm2(grid, wrap2pi(p[1]), exp(p[3]), wrap2pi(p[2]), exp(p[4]), plogis(p[5])))
  }
  if (fit$model == "mix_WC2") {
    return(dmix_wc2(grid, wrap2pi(p[1]), plogis(p[3]), wrap2pi(p[2]), plogis(p[4]), plogis(p[5])))
  }
  if (fit$model == "mix_aFMM2") {
    return(dmix_afmm2(grid,
      wrap2pi(p[1]), plogis(p[3]), plogis(p[5]), wrap2pi(p[7]),
      wrap2pi(p[2]), plogis(p[4]), plogis(p[6]), wrap2pi(p[8]),
      plogis(p[9])
    ))
  }
  if (fit$model == "Fourier_K2") return(dfourier2(grid, p[1], p[2], p[3], p[4]))
  stop("Unknown model: ", fit$model)
}
