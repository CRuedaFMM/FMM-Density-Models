###############################################################
## Circular benchmark: probabilistic mixtures vs functional FMM
## Paper: FMM Density Models for Circular and Toroidal Data
###############################################################


## Required functions from source files:
## d_aFMM
## nll_mix_vm2
## nll_mix_wc2
## nll_mixture_aFMM2
## nll_fourier2
## density_from_fit

## Required before running:
## source("R/circular_density_models.R")
## source("R/circular_likelihoods.R")
##
wrap2pi <- function(x) {
  x %% (2*pi)
}

AIC_from_nll <- function(nll, npar) {
  2*nll + 2*npar
}

###############################################################
## 1. Generic likelihood fitting
###############################################################

fit_one_model <- function(theta,
                          nll_fun,
                          init_fun,
                          npar,
                          model_name,
                          nstart = 100,
                          maxit = 2000) {
  
  theta <- wrap2pi(theta)
  best <- NULL
  
  for (s in seq_len(nstart)) {
    
    opt <- try(
      optim(
        par = init_fun(),
        fn = nll_fun,
        theta = theta,
        method = "Nelder-Mead",
        control = list(maxit = maxit, reltol = 1e-8)
      ),
      silent = TRUE
    )
    
    if (!inherits(opt, "try-error") && is.finite(opt$value)) {
      if (is.null(best) || opt$value < best$value) {
        best <- opt
      }
    }
  }
  
  if (is.null(best)) {
    stop("Model fitting failed for: ", model_name)
  }
  
  list(
    model = model_name,
    raw_par = best$par,
    nll = best$value,
    npar = npar,
    AIC = AIC_from_nll(best$value, npar),
    optim = best
  )
}

###############################################################
## 2. Probabilistic circular models
###############################################################

fit_mixture_vm2 <- function(theta, nstart = 100) {
  
  fit_one_model(
    theta = theta,
    nll_fun = nll_mix_vm2,
    init_fun = function() c(
      runif(1, 0, 2*pi),                 # mu 1
      runif(1, 0, 2*pi),                 # mu 2
      log(runif(1, 0.5, 30)),            # kappa 1
      log(runif(1, 0.5, 30)),            # kappa 2
      qlogis(runif(1, 0.15, 0.85))       # mixture weight
    ),
    npar = 5,
    model_name = "mix_vM2",
    nstart = nstart
  )
}

fit_mixture_wc2 <- function(theta, nstart = 100) {
  
  fit_one_model(
    theta = theta,
    nll_fun = nll_mix_wc2,
    init_fun = function() c(
      runif(1, 0, 2*pi),                 # mu 1
      runif(1, 0, 2*pi),                 # mu 2
      qlogis(runif(1, 0.05, 0.95)),      # gamma 1
      qlogis(runif(1, 0.05, 0.95)),      # gamma 2
      qlogis(runif(1, 0.15, 0.85))       # mixture weight
    ),
    npar = 5,
    model_name = "mix_WC2",
    nstart = nstart
  )
}

fit_mixture_aFMM2 <- function(theta,
                              nstart = 100,
                              omega_min_init = 0.15,
                              omega_max_init = 0.95) {
  
  fit_one_model(
    theta = theta,
    nll_fun = nll_mixture_aFMM2,
    init_fun = function() c(
      runif(1, 0, 2*pi),                              # alpha 1
      runif(1, 0, 2*pi),                              # alpha 2
      qlogis(runif(1, omega_min_init, omega_max_init)), # omega 1
      qlogis(runif(1, omega_min_init, omega_max_init)), # omega 2
      qlogis(runif(1, 0.10, 0.98)),                   # r 1
      qlogis(runif(1, 0.10, 0.98)),                   # r 2
      runif(1, 0, 2*pi),                              # beta 1
      runif(1, 0, 2*pi),                              # beta 2
      qlogis(runif(1, 0.15, 0.85))                    # mixture weight
    ),
    npar = 9,
    model_name = "mix_aFMM2",
    nstart = nstart
  )
}

fit_fourier_density_k2 <- function(theta, nstart = 100) {
  
  fit_one_model(
    theta = theta,
    nll_fun = nll_fourier2,
    init_fun = function() rnorm(4, 0, 0.2),
    npar = 4,
    model_name = "Fourier_K2",
    nstart = nstart
  )
}

fit_all_probabilistic_models <- function(theta, nstart = 100) {
  
  fits <- list(
    fit_mixture_vm2(theta, nstart = nstart),
    fit_mixture_wc2(theta, nstart = nstart),
    fit_mixture_aFMM2(theta, nstart = nstart),
    fit_fourier_density_k2(theta, nstart = nstart)
  )
  
  table <- do.call(
    rbind,
    lapply(fits, function(fit) {
      data.frame(
        model = fit$model,
        npar = fit$npar,
        nll = fit$nll,
        AIC = fit$AIC
      )
    })
  )
  
  table <- table[order(table$AIC), ]
  rownames(table) <- NULL
  
  list(
    fits = fits,
    table = table
  )
}

###############################################################
## 3. Circular KDE target
###############################################################

make_circular_kde <- function(theta,
                              grid,
                              bw = 0.20) {
  
  theta <- wrap2pi(as.numeric(theta))
  theta <- theta[is.finite(theta)]
  
  D <- outer(
    grid,
    theta,
    function(g, x) atan2(sin(g - x), cos(g - x))
  )
  
  K <- exp(-(D^2) / (2*bw^2))
  kde <- rowMeans(K)
  
  delta <- 2*pi / length(grid)
  kde <- pmax(kde, 0)
  kde <- kde / (sum(kde) * delta)
  
  kde
}

###############################################################
## 4. FMM functional approximation to circular KDE
###############################################################

Phi_FMM <- function(theta, alpha, omega) {
  2 * atan(omega * tan((theta - alpha) / 2))
}

design_FMM_functional_1d <- function(grid,
                                     alpha,
                                     omega) {
  
  J <- length(alpha)
  X <- matrix(1, nrow = length(grid), ncol = 1 + 2*J)
  
  col <- 2
  
  for (j in seq_len(J)) {
    Phi_j <- Phi_FMM(grid, alpha[j], omega[j])
    X[, col] <- cos(Phi_j)
    X[, col + 1] <- sin(Phi_j)
    col <- col + 2
  }
  
  X
}

fit_FMM_functional_kde <- function(theta,
                                   J = 2,
                                   n_grid = 512,
                                   nstart = 80,
                                   omega_min = 0.03,
                                   omega_max = 0.95,
                                   bw = 0.20,
                                   seed = 1) {
  
  set.seed(seed)
  
  theta <- wrap2pi(as.numeric(theta))
  theta <- theta[is.finite(theta)]
  
  grid <- seq(0, 2*pi, length.out = n_grid + 1)[-(n_grid + 1)]
  delta <- 2*pi / n_grid
  
  fhat <- make_circular_kde(theta, grid, bw = bw)
  
  objective_phase <- function(par_phase) {
    
    alpha <- wrap2pi(par_phase[seq_len(J)])
    eta <- par_phase[(J + 1):(2*J)]
    
    omega <- omega_min +
      (omega_max - omega_min) * plogis(eta)
    
    X <- design_FMM_functional_1d(grid, alpha, omega)
    beta <- lm.fit(X, fhat)$coefficients
    fitted <- as.numeric(X %*% beta)
    
    sum((fitted - fhat)^2) * delta
  }
  
  best <- NULL
  
  for (s in seq_len(nstart)) {
    
    init_alpha <- runif(J, 0, 2*pi)
    init_omega <- runif(J, 0.10, 0.80)
    
    init_eta <- qlogis(
      (init_omega - omega_min) /
        (omega_max - omega_min)
    )
    
    init <- c(init_alpha, init_eta)
    
    opt <- try(
      optim(
        par = init,
        fn = objective_phase,
        method = "Nelder-Mead",
        control = list(maxit = 1500, reltol = 1e-9)
      ),
      silent = TRUE
    )
    
    if (!inherits(opt, "try-error") && is.finite(opt$value)) {
      if (is.null(best) || opt$value < best$value) {
        best <- opt
      }
    }
  }
  
  if (is.null(best)) {
    stop("Functional FMM fit failed")
  }
  
  alpha <- wrap2pi(best$par[seq_len(J)])
  eta <- best$par[(J + 1):(2*J)]
  
  omega <- omega_min +
    (omega_max - omega_min) * plogis(eta)
  
  X <- design_FMM_functional_1d(grid, alpha, omega)
  beta <- lm.fit(X, fhat)$coefficients
  
  fitted_raw <- as.numeric(X %*% beta)
  
  fitted_pos <- pmax(fitted_raw, 0)
  fitted_pos <- fitted_pos / (sum(fitted_pos) * delta)
  
  list(
    grid = grid,
    fhat = fhat,
    fitted_raw = fitted_raw,
    fitted_pos = fitted_pos,
    alpha = alpha,
    omega = omega,
    beta = beta,
    ISE_raw = sum((fitted_raw - fhat)^2) * delta,
    ISE_pos = sum((fitted_pos - fhat)^2) * delta,
    opt = best
  )
}

###############################################################
## 5. Fourier functional approximation to circular KDE
###############################################################

fit_fourier_kde <- function(fhat,
                            grid,
                            K = 4) {
  
  X <- matrix(1, nrow = length(grid), ncol = 1 + 2*K)
  
  col <- 2
  
  for (k in seq_len(K)) {
    X[, col] <- cos(k * grid)
    X[, col + 1] <- sin(k * grid)
    col <- col + 2
  }
  
  beta <- lm.fit(X, fhat)$coefficients
  fitted_raw <- as.numeric(X %*% beta)
  
  delta <- 2*pi / length(grid)
  
  fitted_pos <- pmax(fitted_raw, 0)
  fitted_pos <- fitted_pos / (sum(fitted_pos) * delta)
  
  list(
    K = K,
    beta = beta,
    fitted_raw = fitted_raw,
    fitted_pos = fitted_pos,
    ISE = sum((fitted_pos - fhat)^2) * delta
  )
}

###############################################################
## 6. Evaluation metrics
###############################################################

ISE_grid <- function(fitted,
                     target,
                     grid) {
  
  delta <- 2*pi / length(grid)
  
  sum((fitted - target)^2) * delta
}

R2_grid <- function(fitted,
                    target) {
  
  1 - sum((target - fitted)^2) /
    sum((target - mean(target))^2)
}

###############################################################
## 7. Main comparison function
###############################################################

compare_circular_functional_models <- function(theta,
                                               dataset_name = "dataset",
                                               bw = 0.20,
                                               n_grid = 512,
                                               nstart_prob = 100,
                                               nstart_FMM = 80,
                                               seed = 1) {
  
  theta <- wrap2pi(as.numeric(theta))
  theta <- theta[is.finite(theta)]
  
  ## Probabilistic models
  prob <- fit_all_probabilistic_models(
    theta,
    nstart = nstart_prob
  )
  
  ## Common functional target
  grid <- seq(0, 2*pi, length.out = n_grid + 1)[-(n_grid + 1)]
  fhat <- make_circular_kde(theta, grid, bw = bw)
  
  ## FMM functional approximation
  FMM_functional <- fit_FMM_functional_kde(
    theta = theta,
    J = 2,
    n_grid = n_grid,
    nstart = nstart_FMM,
    bw = bw,
    seed = seed
  )
  
  ## Fourier functional approximations
  fourier_k2 <- fit_fourier_kde(fhat, grid, K = 2)
  fourier_k4 <- fit_fourier_kde(fhat, grid, K = 4)
  
  ## Evaluation table
  table <- data.frame(
    dataset = dataset_name,
    model = c(
      "FMM_functional",
      "Fourier_K4",
      "mix_aFMM2",
      "mix_WC2",
      "mix_vM2",
      "Fourier_K2"
    ),
    ISE = c(
      FMM_functional$ISE_pos,
      fourier_k4$ISE,
      NA,
      NA,
      NA,
      fourier_k2$ISE
    ),
    R2 = c(
      R2_grid(FMM_functional$fitted_pos, fhat),
      R2_grid(fourier_k4$fitted_pos, fhat),
      NA,
      NA,
      NA,
      R2_grid(fourier_k2$fitted_pos, fhat)
    ),
    AIC = NA_real_
  )
  
  ## Add AIC values for probabilistic models
  for (m in c("mix_aFMM2", "mix_WC2", "mix_vM2", "Fourier_K2")) {
    
    idx <- table$model == m
    aic_value <- prob$table$AIC[prob$table$model == m]
    
    if (length(aic_value) == 1) {
      table$AIC[idx] <- aic_value
    }
  }
  
  table <- table[order(table$ISE, na.last = TRUE), ]
  rownames(table) <- NULL
  
  list(
    table = table,
    probabilistic = prob,
    FMM_functional = FMM_functional,
    fourier_k2 = fourier_k2,
    fourier_k4 = fourier_k4
  )
}