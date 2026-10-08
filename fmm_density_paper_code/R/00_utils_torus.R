############################################################
## General utilities for circular/toroidal analyses
############################################################

wrap_2pi <- function(x) x %% (2*pi)
wrap2pi <- wrap_2pi
vm_wrap2pi <- wrap_2pi

safe_exp <- function(z) exp(pmin(z, 700))
clip_pos <- function(x, eps = 1e-8) pmax(x, eps)

sigmoid <- function(u) 1/(1 + exp(-u))
logit <- function(p) log(p/(1-p))

omega_from_u01 <- function(u, eps = 1e-6){
  eps + (1 - 2*eps) * sigmoid(u)
}

u_from_omega01 <- function(omega, eps = 1e-6){
  z <- (omega - eps)/(1 - 2*eps)
  z <- pmin(pmax(z, 1e-12), 1 - 1e-12)
  logit(z)
}

pack_theta01 <- function(alpha_P, omega_P, alpha_T, omega_T, eps = 1e-6){
  c(alpha_P,
    u_from_omega01(omega_P, eps = eps),
    alpha_T,
    u_from_omega01(omega_T, eps = eps))
}

unpack_theta01 <- function(theta, K, L, eps = 1e-6){
  alpha_P <- theta[1:K]
  uP <- theta[(K+1):(2*K)]
  alpha_T <- theta[(2*K+1):(2*K+L)]
  uT <- theta[(2*K+L+1):(2*K+2*L)]
  list(
    alpha_P = wrap_2pi(alpha_P),
    omega_P = omega_from_u01(uP, eps = eps),
    alpha_T = wrap_2pi(alpha_T),
    omega_T = omega_from_u01(uT, eps = eps)
  )
}

circ_diff <- function(a, b){
  ((a - b + pi) %% (2*pi)) - pi
}

ang_diff <- function(a, b){
  d <- abs(a - b)
  pmin(d, 2*pi - d)
}

circ_sqdist <- function(a, b) sum(circ_diff(a, b)^2)
sqdist <- function(a, b) sum((a - b)^2)

make_grid <- function(n_theta, n_phi){
  theta_grid <- seq(0, 2*pi, length.out = n_theta + 1)[-(n_theta + 1)]
  phi_grid <- seq(0, 2*pi, length.out = n_phi + 1)[-(n_phi + 1)]
  list(theta = theta_grid, phi = phi_grid)
}

grid_torus <- function(ng = 100){
  phi <- seq(0, 2*pi, length.out = ng + 1)[-(ng + 1)]
  psi <- seq(0, 2*pi, length.out = ng + 1)[-(ng + 1)]
  as.matrix(expand.grid(phi = phi, psi = psi))
}

grid.torus <- grid_torus

normalize_mat <- function(Z){
  Z <- as.matrix(Z)
  Z[!is.finite(Z)] <- 0
  Z[Z < 0] <- 0
  s <- sum(Z)
  if (!is.finite(s) || s <= 0) stop("Matrix cannot be normalized: non-positive total mass.")
  Z/s
}
############################################################
## Sampling pairs from a toroidal probability surface
############################################################

sample_pairs_from_Z <- function(Z, n) {

  Z <- normalize_mat(Z)

  nphi <- nrow(Z)
  npsi <- ncol(Z)

  idx <- sample(
    length(Z),
    size = n,
    replace = TRUE,
    prob = as.vector(Z)
  )

  i <- ((idx - 1) %% nphi) + 1
  j <- ((idx - 1) %/% nphi) + 1

  phi <- seq(0, 2*pi, length.out = nphi)[i]
  psi <- seq(0, 2*pi, length.out = npsi)[j]

  data.frame(
    phi = phi,
    psi = psi
  )
}
normalize_grid <- function(M){
  M <- as.matrix(M)
  M[!is.finite(M)] <- 0
  M[M < 0] <- 0
  dA <- (2*pi/nrow(M)) * (2*pi/ncol(M))
  s <- sum(M * dA)
  if (!is.finite(s) || s <= 0) stop("Grid cannot be normalized: non-positive integral.")
  M/s
}

R2_2D <- function(Z_true, Z_fit){
  1 - sum((Z_true - Z_fit)^2) / sum((Z_true - mean(Z_true))^2)
}

L2_sq_grid <- function(P, Q){
  dA <- (2*pi/nrow(P)) * (2*pi/ncol(P))
  sum((P - Q)^2) * dA
}

############################################################
## Integrated squared error against known generating surface
############################################################

grid_integral <- function(Z){
  Z <- as.matrix(Z)
  dA <- (2*pi/nrow(Z)) * (2*pi/ncol(Z))
  sum(Z) * dA
}

as_density_grid <- function(Z){
  Z <- as.matrix(Z)

  mass <- grid_integral(Z)

  if(!is.finite(mass) || abs(mass) < 1e-12)
    stop("Surface has zero or non-finite integral.")

  Z / mass
}

ISE_density_grid <- function(Z_true, Z_fit){

  if(!all(dim(Z_true) == dim(Z_fit)))
    stop("Z_true and Z_fit must have the same dimensions.")

  P <- as_density_grid(Z_true)
  Q <- as_density_grid(Z_fit)

  dA <- (2*pi/nrow(P)) * (2*pi/ncol(P))

  sum((P - Q)^2) * dA
}


cv_circular <- function(theta, na.rm = TRUE){
  if (na.rm) theta <- theta[is.finite(theta)]
  if (length(theta) < 2) return(NA_real_)
  C <- mean(cos(theta)); S <- mean(sin(theta))
  R <- sqrt(C^2 + S^2)
  if (R <= 0) return(Inf)
  sqrt(-2 * log(R))
}

mean_circular_2pi <- function(x){
  m <- atan2(mean(sin(x)), mean(cos(x)))
  if (m < 0) m <- m + 2*pi
  m
}

order_pi <- function(x){
  x_shift <- (x - pi) %% (2*pi)
  x[order(x_shift)]
}

logsumexp <- function(v){
  m <- max(v)
  m + log(sum(exp(v - m)))
}
