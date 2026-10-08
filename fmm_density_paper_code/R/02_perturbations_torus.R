############################################################
## Toroidal perturbation schemes
############################################################

shift_torus_mat <- function(Z, dphi = 0, dpsi = 0){
  nr <- nrow(Z); nc <- ncol(Z)
  i_shift <- ((seq_len(nr) - 1 - dphi) %% nr) + 1
  j_shift <- ((seq_len(nc) - 1 - dpsi) %% nc) + 1
  Z[i_shift, j_shift]
}

smooth_torus_mat <- function(Z, kernel = NULL){
  nr <- nrow(Z); nc <- ncol(Z)
  if(is.null(kernel)){
    kernel <- matrix(c(1,2,1,2,4,2,1,2,1), nrow = 3, byrow = TRUE)
    kernel <- kernel/sum(kernel)
  }
  out <- matrix(0, nr, nc)
  for(a in -1:1){
    for(b in -1:1){
      out <- out + kernel[a+2, b+2] * shift_torus_mat(Z, dphi = a, dpsi = b)
    }
  }
  normalize_mat(out)
}

perturb_Z_torus <- function(Z,
                            shift_sd = 3,
                            alpha_mix = 0.12,
                            smooth_iter = 1,
                            noise_level = 0.03){
  nr <- nrow(Z); nc <- ncol(Z)
  dphi <- round(rnorm(1, 0, shift_sd))
  dpsi <- round(rnorm(1, 0, shift_sd))
  Z_shift <- shift_torus_mat(Z, dphi = dphi, dpsi = dpsi)
  Z_new <- (1 - alpha_mix)*Z + alpha_mix*Z_shift
  eps <- matrix(rlnorm(nr*nc, meanlog = 0, sdlog = noise_level), nr, nc)
  Z_new <- Z_new * eps
  for(s in seq_len(smooth_iter)) Z_new <- smooth_torus_mat(Z_new)
  normalize_mat(Z_new)
}

wrap_angle <- function(x) x %% (2*pi)

torus_bump <- function(G, mu1, mu2, sd1 = 0.25, sd2 = 0.25, rho = 0){
  d1 <- ang_diff(G[,1], mu1)
  d2 <- ang_diff(G[,2], mu2)
  if(abs(rho) < 1e-10){
    z <- exp(-0.5*((d1/sd1)^2 + (d2/sd2)^2))
  } else {
    q <- (d1/sd1)^2 - 2*rho*(d1/sd1)*(d2/sd2) + (d2/sd2)^2
    z <- exp(-q/(2*(1-rho^2)))
  }
  z
}

perturb_Z_torus2 <- function(Z,
                             shift_sd = 0.15,
                             alpha_mix = 0.10,
                             noise_level = 0.03,
                             smooth_iter = 1,
                             shear_sd = 0.10,
                             n_bumps = sample(1:3, 1),
                             bump_weight = 0.15,
                             bump_sd = 0.30){
  ng <- nrow(Z)
  G <- grid_torus(ng)
  delta1 <- rnorm(1, 0, shift_sd)
  delta2 <- rnorm(1, 0, shift_sd)
  shear <- rnorm(1, 0, shear_sd)

  Gp <- cbind(
    wrap_angle(G[,1] + delta1 + shear*sin(G[,2])),
    wrap_angle(G[,2] + delta2 + shear*sin(G[,1]))
  )

  idx1 <- pmax(1, pmin(ng, round(Gp[,1]/(2*pi)*(ng - 1)) + 1))
  idx2 <- pmax(1, pmin(ng, round(Gp[,2]/(2*pi)*(ng - 1)) + 1))

  Z_shift <- matrix(0, ng, ng)
  for(k in seq_len(length(idx1))) Z_shift[k] <- Z[idx1[k], idx2[k]]
  Z_shift <- normalize_mat(matrix(Z_shift, ng, ng))

  base_part <- normalize_mat((1 - alpha_mix)*Z + alpha_mix*Z_shift)

  if(n_bumps > 0 && bump_weight > 0){
    bump_field <- rep(0, nrow(G))
    for(j in seq_len(n_bumps)){
      bump_field <- bump_field + torus_bump(
        G,
        mu1 = runif(1, 0, 2*pi),
        mu2 = runif(1, 0, 2*pi),
        sd1 = bump_sd,
        sd2 = bump_sd,
        rho = runif(1, -0.4, 0.4)
      )
    }
    bump_field <- normalize_mat(matrix(bump_field, ng, ng))
    Z_mix <- normalize_mat((1 - bump_weight)*base_part + bump_weight*bump_field)
  } else {
    Z_mix <- base_part
  }

  E <- normalize_mat(matrix(rexp(ng*ng, rate = 1), ng, ng))
  Z_out <- normalize_mat((1 - noise_level)*Z_mix + noise_level*E)

  if(smooth_iter > 0){
    for(it in seq_len(smooth_iter)) Z_out <- smooth_torus_mat(Z_out)
  }
  Z_out
}

generate_replicates <- function(Z_true,
                                B = 100,
                                n = 300,
                                kappa = 25,
                                ng = nrow(Z_true),
                                shift_sd = 3,
                                alpha_mix = 0.12,
                                smooth_iter = 1,
                                noise_level = 0.03){
  requireNamespace("ClusTorus")
  G <- grid_torus(ng)
  out <- vector("list", B)
  for(b in seq_len(B)){
    Z_b <- perturb_Z_torus(Z_true, shift_sd, alpha_mix, smooth_iter, noise_level)
    dat <- sample_pairs_from_Z(Z_b, n)
    fhat <- ClusTorus::kde.torus(as.matrix(dat), eval.point = G, concentration = kappa)
    M <- normalize_mat(matrix(fhat, nrow = ng, ncol = ng))
    out[[b]] <- list(dat = dat, Z_b = Z_b, M = M)
  }
  out
}

generate_replicates2 <- function(Z,
                                 B = 100,
                                 n = 300,
                                 kappa = 25,
                                 ng = nrow(Z),
                                 shift_sd = 0.15,
                                 alpha_mix = 0.12,
                                 smooth_iter = 1,
                                 noise_level = 0.03,
                                 shear_sd = 0.10,
                                 n_bumps = NULL,
                                 bump_weight = 0.15,
                                 bump_sd = 0.30){
  requireNamespace("ClusTorus")
  G <- grid_torus(ng)
  out <- vector("list", B)
  for(b in seq_len(B)){
    nb <- if(is.null(n_bumps)) sample(1:3, 1) else n_bumps
    Z_b <- perturb_Z_torus2(Z, shift_sd, alpha_mix, noise_level, smooth_iter,
                            shear_sd, nb, bump_weight, bump_sd)
    dat <- sample_pairs_from_Z(Z_b, n)
    fhat <- ClusTorus::kde.torus(as.matrix(dat), eval.point = G, concentration = kappa)
    M <- normalize_mat(matrix(fhat, nrow = ng, ncol = ng))
    out[[b]] <- list(dat = dat, Z_b = Z_b, M = M,
                     settings = list(shift_sd = shift_sd, alpha_mix = alpha_mix,
                                     noise_level = noise_level, shear_sd = shear_sd,
                                     n_bumps = nb, bump_sd = bump_sd))
  }
  out
}
