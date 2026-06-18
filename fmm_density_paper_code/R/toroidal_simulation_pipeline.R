###############################################################
## Toroidal simulation pipeline: Fourier33, S-FMM2, VM4
###############################################################

check_sfmm_available <- function() {
  if (!exists("S_fmm", mode = "function")) {
    stop("The function S_fmm() is not available. Source the core S-FMM implementation before running this script.")
  }
}

fit_one_toroidal_replicate <- function(M, theta_grid, phi_grid, N_vm = 3000, seed = 1) {
  check_sfmm_available()
  fourier_fit <- fit_fourier_surface_2d(M, theta_grid, phi_grid, K = 3, L = 3)
  Z_fourier <- fourier_fit$fitted_matrix
  sfmm_fit <- S_fmm(M, 30, 2, 2)
  vm_fit <- fit_mix_L2_from_FF(M, Mi = 4, N = N_vm, seed = seed)
  list(
    R2_Fourier33 = R2_2D(M, Z_fourier),
    R2_SFMM2 = sfmm_fit$r2,
    R2_VM4 = R2_2D(normalize_grid(M), vm_fit$Q),
    sfmm_fit = sfmm_fit,
    vm_fit = vm_fit
  )
}

summarize_toroidal_replicates <- function(Z_ref, replicas, N_vm = 3000, seed = 1) {
  check_sfmm_available()
  grids <- make_grid(nrow(Z_ref), ncol(Z_ref))
  theta_grid <- grids$theta; phi_grid <- grids$phi
  B <- length(replicas)
  R2_mat <- matrix(NA_real_, B, 3)
  colnames(R2_mat) <- c("Fourier33", "S-FMM2", "VM4")
  SFMM_par_mat <- matrix(NA_real_, B, 4)
  VM_centers_mat <- matrix(NA_real_, B, 8)

  ref_fit <- S_fmm(Z_ref, 30, 2, 2)
  SFMM_ref <- ref_fit$par_current[c(1:2, 5:6)]
  SFMM_ref <- c(order_pi(SFMM_ref[1:2]), order_pi(SFMM_ref[3:4]))

  for (b in seq_len(B)) {
    M <- replicas[[b]]$M
    res <- fit_one_toroidal_replicate(M, theta_grid, phi_grid, N_vm = N_vm, seed = seed)
    R2_mat[b,] <- c(res$R2_Fourier33, res$R2_SFMM2, res$R2_VM4)
    SFMM_par_mat[b,] <- res$sfmm_fit$par_current[c(1:2, 5:6)]
    VM_centers_mat[b,] <- c(res$vm_fit$par$mu1, res$vm_fit$par$mu2)
  }

  SFMM_par_mat <- cbind(t(apply(SFMM_par_mat[,1:2, drop=FALSE], 1, order_pi)),
                        t(apply(SFMM_par_mat[,3:4, drop=FALSE], 1, order_pi)))
  VM_centers_mat <- cbind(t(apply(VM_centers_mat[,1:4, drop=FALSE], 1, order_pi)),
                          t(apply(VM_centers_mat[,5:8, drop=FALSE], 1, order_pi)))

  CD_SFMM <- apply(SFMM_par_mat, 2, cv_circular)
  CD_VM <- apply(VM_centers_mat, 2, cv_circular)
  mean_SFMM <- apply(SFMM_par_mat, 2, mean_circular_2pi)
  AD_SFMM <- mean(mapply(ang_diff, mean_SFMM, SFMM_ref), na.rm = TRUE)

  list(
    mean_R2 = colMeans(R2_mat, na.rm = TRUE),
    CD_SFMM = mean(CD_SFMM, na.rm = TRUE),
    CD_VM4 = mean(CD_VM, na.rm = TRUE),
    AD_SFMM = AD_SFMM,
    R2_mat = R2_mat,
    SFMM_par_mat = SFMM_par_mat,
    VM_centers_mat = VM_centers_mat,
    SFMM_ref = SFMM_ref
  )
}

one_row_summary <- function(scenario, res) {
  data.frame(
    Scenario = scenario,
    R2_Fourier33 = unname(res$mean_R2["Fourier33"]),
    R2_SFMM2 = unname(res$mean_R2["S-FMM2"]),
    R2_VM4 = unname(res$mean_R2["VM4"]),
    CD_SFMM2 = res$CD_SFMM,
    CD_VM4 = res$CD_VM4,
    AD_SFMM2 = res$AD_SFMM
  )
}
