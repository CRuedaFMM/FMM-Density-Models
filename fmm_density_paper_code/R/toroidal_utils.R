###############################################################
## Toroidal utilities
###############################################################

wrap2pi <- function(x) x %% (2*pi)
wrap_angle <- wrap2pi

make_grid <- function(n_theta, n_phi) {
  theta_grid <- seq(0, 2*pi, length.out = n_theta + 1)[-(n_theta + 1)]
  phi_grid   <- seq(0, 2*pi, length.out = n_phi + 1)[-(n_phi + 1)]
  list(theta = theta_grid, phi = phi_grid)
}

grid_torus <- function(ng = 100) {
  phi <- seq(0, 2*pi, length.out = ng + 1)[-(ng + 1)]
  psi <- seq(0, 2*pi, length.out = ng + 1)[-(ng + 1)]
  as.matrix(expand.grid(phi = phi, psi = psi))
}

normalize_mat <- function(M) {
  M <- as.matrix(M)
  M[!is.finite(M)] <- 0
  M[M < 0] <- 0
  s <- sum(M)
  if (!is.finite(s) || s <= 0) stop("Matrix has non-positive total mass.")
  M / s
}

normalize_grid <- function(M) {
  M <- as.matrix(M)
  dA <- (2*pi/nrow(M)) * (2*pi/ncol(M))
  M[!is.finite(M)] <- 0
  M <- pmax(M, 0)
  M / sum(M * dA)
}

R2_2D <- function(Z_true, Z_fit) {
  1 - sum((Z_true - Z_fit)^2) / sum((Z_true - mean(Z_true))^2)
}

ang_diff <- function(a, b) {
  d <- abs(a - b)
  pmin(d, 2*pi - d)
}

cv_circular <- function(theta, na.rm = TRUE) {
  if (na.rm) theta <- theta[is.finite(theta)]
  if (length(theta) < 2) return(NA_real_)
  C <- mean(cos(theta)); S <- mean(sin(theta))
  R <- sqrt(C^2 + S^2)
  if (R <= 0) return(Inf)
  sqrt(-2 * log(R))
}

mean_circular_2pi <- function(x) {
  m <- atan2(mean(sin(x)), mean(cos(x)))
  if (m < 0) m <- m + 2*pi
  m
}

order_pi <- function(x) {
  x_shift <- (x - pi) %% (2*pi)
  x[order(x_shift)]
}

sample_pairs_from_Z <- function(Z, n) {
  Z <- normalize_mat(Z)
  nphi <- nrow(Z); npsi <- ncol(Z)
  idx <- sample(length(Z), size = n, replace = TRUE, prob = as.vector(Z))
  i <- ((idx - 1) %% nphi) + 1
  j <- ((idx - 1) %/% nphi) + 1
  phi <- seq(0, 2*pi, length.out = nphi + 1)[-(nphi + 1)][i]
  psi <- seq(0, 2*pi, length.out = npsi + 1)[-(npsi + 1)][j]
  data.frame(phi = phi, psi = psi)
}

sample_from_grid <- function(P, theta_grid, phi_grid, N = 3000, seed = 1) {
  set.seed(seed)
  P <- normalize_mat(P)
  nr <- nrow(P); nc <- ncol(P)
  idx <- sample.int(nr*nc, size = N, replace = TRUE, prob = as.vector(P))
  i <- ((idx - 1) %% nr) + 1
  j <- ((idx - 1) %/% nr) + 1
  cbind(theta_grid[i], phi_grid[j])
}

L2_sq_grid <- function(P, Q) {
  dA <- (2*pi/nrow(P)) * (2*pi/ncol(P))
  sum((P - Q)^2) * dA
}
