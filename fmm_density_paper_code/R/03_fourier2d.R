############################################################
## Two-dimensional Fourier surface approximation
############################################################

build_fourier_design_2d <- function(theta_vec, phi_vec, K = 3, L = 3){
  X_list <- list("(Intercept)" = rep(1, length(theta_vec)))

  for(k in 1:K){
    X_list[[paste0("cos_t_", k)]] <- cos(k*theta_vec)
    X_list[[paste0("sin_t_", k)]] <- sin(k*theta_vec)
  }
  for(l in 1:L){
    X_list[[paste0("cos_p_", l)]] <- cos(l*phi_vec)
    X_list[[paste0("sin_p_", l)]] <- sin(l*phi_vec)
  }
  for(k in 1:K){
    for(l in 1:L){
      X_list[[paste0("cc_", k, "_", l)]] <- cos(k*theta_vec)*cos(l*phi_vec)
      X_list[[paste0("sc_", k, "_", l)]] <- sin(k*theta_vec)*cos(l*phi_vec)
      X_list[[paste0("cs_", k, "_", l)]] <- cos(k*theta_vec)*sin(l*phi_vec)
      X_list[[paste0("ss_", k, "_", l)]] <- sin(k*theta_vec)*sin(l*phi_vec)
    }
  }
  as.matrix(as.data.frame(X_list))
}

fit_fourier_surface_2d <- function(Z, theta_grid, phi_grid, K = 3, L = 3){
  Gt <- length(theta_grid); Gp <- length(phi_grid)
  grid_df <- expand.grid(theta = theta_grid, phi = phi_grid)
  y <- as.vector(Z)
  X <- build_fourier_design_2d(grid_df$theta, grid_df$phi, K, L)
  fit <- lm.fit(x = X, y = y)
  yhat <- as.vector(X %*% fit$coefficients)
  list(coefficients = fit$coefficients,
       fitted_matrix = matrix(yhat, nrow = Gt, ncol = Gp),
       K = K, L = L)
}
