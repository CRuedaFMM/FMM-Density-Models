############################################################
## Synthetic toroidal density generators
############################################################

generate_2modes1 <- function(n = 100){
  phi <- seq(-pi, pi, length.out = n)
  psi <- seq(-pi, pi, length.out = n)
  grid <- expand.grid(phi = phi, psi = psi)
  wrap_dist <- function(x, mu) ((x - mu + pi) %% (2*pi)) - pi
  gauss2d_torus <- function(phi, psi, mu, Sigma){
    D <- cbind(wrap_dist(phi, mu[1]), wrap_dist(psi, mu[2]))
    exp(-0.5*rowSums((D %*% solve(Sigma))*D))
  }
  mu1 <- c(-2.2, -2.0)
  Sigma1 <- matrix(c(0.18, 0.14, 0.14, 0.22), 2, 2)
  mu2 <- c(1.0, -1.1)
  Sigma2 <- matrix(c(0.22, -0.15, -0.15, 0.18), 2, 2)
  f <- 0.52*gauss2d_torus(grid$phi, grid$psi, mu1, Sigma1) +
    0.48*gauss2d_torus(grid$phi, grid$psi, mu2, Sigma2)
  Z <- matrix(f/sum(f), nrow = n, ncol = n)
  list(phi = phi, psi = psi, Z = Z, grid = grid,
       mu = rbind(mu1, mu2), Sigma = list(Sigma1, Sigma2), weights = c(0.52, 0.48))
}

generate_2modes2 <- function(n = 100){
  phi <- seq(-pi, pi, length.out = n)
  psi <- seq(-pi, pi, length.out = n)
  grid <- expand.grid(phi = phi, psi = psi)
  wrap_dist <- function(x, mu) ((x - mu + pi) %% (2*pi)) - pi
  gauss2d_torus <- function(phi, psi, mu, Sigma){
    D <- cbind(wrap_dist(phi, mu[1]), wrap_dist(psi, mu[2]))
    exp(-0.5*rowSums((D %*% solve(Sigma))*D))
  }
  mu1 <- c(-1.6, -1.9)
  Sigma1 <- matrix(c(0.22, 0.14, 0.14, 0.26), 2, 2)
  mu2 <- c(0.2, -1.1)
  Sigma2 <- matrix(c(0.26, -0.12, -0.12, 0.22), 2, 2)
  f <- 0.55*gauss2d_torus(grid$phi, grid$psi, mu1, Sigma1) +
    0.45*gauss2d_torus(grid$phi, grid$psi, mu2, Sigma2)
  Z <- matrix(f/sum(f), nrow = n, ncol = n)
  list(phi = phi, psi = psi, Z = Z, grid = grid,
       mu = rbind(mu1, mu2), Sigma = list(Sigma1, Sigma2), weights = c(0.55, 0.45))
}
