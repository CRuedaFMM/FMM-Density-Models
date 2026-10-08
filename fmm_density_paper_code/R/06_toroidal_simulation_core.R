############################################################
## Core functions to run toroidal simulation experiments
############################################################

fit_one_toroidal_replicate <- function(M,
                                       theta_grid,
                                       phi_grid,
                                       fourier_K = 3,
                                       fourier_L = 3,
                                       sfmm_fun = S_fmm,
                                       sfmm_K = 2,
                                       sfmm_L = 2,
                                       vm_M = 4,
                                       vm_N = 3000,
                                       vm_seed = 1){
  fourier_fit <- fit_fourier_surface_2d(M, theta_grid, phi_grid, K = fourier_K, L = fourier_L)
  Z_fourier <- fourier_fit$fitted_matrix
  R2_fourier <- R2_2D(M, Z_fourier)

  sfmm_fit <- sfmm_fun(M, 30, sfmm_K, sfmm_L)
  R2_sfmm <- sfmm_fit$r2
  sfmm_par <- sfmm_fit$par_current[c(1:2, 5:6)]
 
 Z_sfmm <- sfmm_fit$y_hat
vm_fit <- fit_mix_L2_from_FF(M, Mi = vm_M, N = vm_N, seed = vm_seed)

Q <- vm_fit$Q

dA <- (2*pi/nrow(M)) * (2*pi/ncol(M))
mass_M <- sum(M) * dA

Z_vm <- Q * mass_M

R2_vm <- R2_2D(M, Z_vm)

list(
    R2 = c(
      Fourier33 = R2_fourier,
      S_FMM2 = R2_sfmm,
      VM4 = R2_vm
    ),

    sfmm_par = sfmm_par,

    vm_centers = c(
      vm_fit$par$mu1,
      vm_fit$par$mu2
    ),

    fits = list(
      fourier = fourier_fit,
      sfmm = sfmm_fit,
      vm4 = vm_fit
    ),

    surfaces = list(
      fourier = Z_fourier,
      sfmm = Z_sfmm,
      vm4 = Z_vm
    )
)
}

summarize_toroidal_replicates <- function(results, S_FMM_reference){
  B <- length(results)
  R2_mat <- do.call(rbind, lapply(results, `[[`, "R2"))
  
ISE_mat <- do.call(
  rbind,
  lapply(results, `[[`, "ISE")
)

min_SFMM_vec <- sapply(
  results,
  `[[`,
  "min_SFMM"
)

prop_negative_SFMM_vec <- sapply(
  results,
  `[[`,
  "prop_negative_SFMM"
)

 SFMM_par_mat <- do.call(rbind, lapply(results, `[[`, "sfmm_par"))
  VM_centers_mat <- do.call(rbind, lapply(results, `[[`, "vm_centers"))

  S_ref <- c(order_pi(S_FMM_reference[1:2]), order_pi(S_FMM_reference[3:4]))
  SFMM_par_mat <- cbind(t(apply(SFMM_par_mat[,1:2, drop = FALSE], 1, order_pi)),
                        t(apply(SFMM_par_mat[,3:4, drop = FALSE], 1, order_pi)))
  VM_centers_mat <- cbind(t(apply(VM_centers_mat[,1:4, drop = FALSE], 1, order_pi)),
                          t(apply(VM_centers_mat[,5:8, drop = FALSE], 1, order_pi)))

  CD_SFMM <- apply(SFMM_par_mat, 2, cv_circular)
  CD_VM <- apply(VM_centers_mat, 2, cv_circular)
  mean_SFMM <- apply(SFMM_par_mat, 2, mean_circular_2pi)
  AD_SFMM <- mean(mapply(ang_diff, mean_SFMM, S_ref))

  list(

  ## R2 against KDE
  mean_R2 = colMeans(R2_mat, na.rm = TRUE),
  sd_R2 = apply(R2_mat, 2, sd, na.rm = TRUE),

  ## ISE against true perturbed surface Z_b
  mean_ISE = colMeans(ISE_mat, na.rm = TRUE),
  sd_ISE = apply(ISE_mat, 2, sd, na.rm = TRUE),

  ## Existing structural summaries
  CD_SFMM = mean(CD_SFMM, na.rm = TRUE),
  CD_VM4 = mean(CD_VM, na.rm = TRUE),
  AD_SFMM = AD_SFMM,

  ## Negativity diagnostics
  min_SFMM = min(min_SFMM_vec, na.rm = TRUE),

  mean_prop_negative_SFMM =
    mean(prop_negative_SFMM_vec, na.rm = TRUE),

  prop_replicates_negative =
    mean(min_SFMM_vec < 0, na.rm = TRUE),

  ## Raw replicate results
  R2_mat = R2_mat,
  ISE_mat = ISE_mat,
  min_SFMM_vec = min_SFMM_vec,
  prop_negative_SFMM_vec = prop_negative_SFMM_vec,
  SFMM_par_mat = SFMM_par_mat,
  VM_centers_mat = VM_centers_mat
)
}

run_toroidal_experiment <- function(Z,
                                    replicas,
                                    sfmm_fun = S_fmm,
                                    fourier_K = 3,
                                    fourier_L = 3,
                                    sfmm_K = 2,
                                    sfmm_L = 2,
                                    vm_M = 4,
                                    vm_N = 3000,
                                    vm_seed = 1){
  grids <- make_grid(nrow(Z), ncol(Z))
  theta_grid <- grids$theta
  phi_grid <- grids$phi

  ref_fit <- sfmm_fun(Z, 30, sfmm_K, sfmm_L)
  S_ref <- ref_fit$par_current[c(1:2, 5:6)]
results <- vector("list", length(replicas))

for(b in seq_along(replicas)){

  message("Replicate ", b, "/", length(replicas))

  one <- fit_one_toroidal_replicate(
    M = replicas[[b]]$M,
    theta_grid = theta_grid,
    phi_grid = phi_grid,
    fourier_K = fourier_K,
    fourier_L = fourier_L,
    sfmm_fun = sfmm_fun,
    sfmm_K = sfmm_K,
    sfmm_L = sfmm_L,
    vm_M = vm_M,
    vm_N = vm_N,
    vm_seed = vm_seed
  )

  Z_b <- replicas[[b]]$Z_b

  one$ISE <- c(
    Fourier33 = ISE_density_grid(
      Z_b,
      one$surfaces$fourier
    ),
    S_FMM2 = ISE_density_grid(
      Z_b,
      one$surfaces$sfmm
    ),
    VM4 = ISE_density_grid(
      Z_b,
      one$surfaces$vm4
    )
  )

  one$min_SFMM <- min(
    one$surfaces$sfmm,
    na.rm = TRUE
  )

  one$prop_negative_SFMM <- mean(
    one$surfaces$sfmm < 0,
    na.rm = TRUE
  )

  results[[b]] <- one
}

  summary <- summarize_toroidal_replicates(results, S_ref)
  list(reference_fit = ref_fit, replicate_results = results, summary = summary)
}
