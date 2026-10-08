
phi_mobius <- function(tau, alpha, omega) {
  2 * atan(omega * tan((tau - alpha) / 2))
}

fit_fmm_linear5 <- function(theta, phi, f_data,
                            alpha_P, omega_P,
                            alpha_T, omega_T,
                            M = 0) {

  K <- length(alpha_P)
  L <- length(alpha_T)

  n_theta <- length(theta)
  n_phi <- length(phi)

  Phi_P <- sapply(1:K, function(k)
    phi_mobius(theta, alpha_P[k], omega_P[k]))

  Phi_T <- sapply(1:L, function(l)
    phi_mobius(phi, alpha_T[l], omega_T[l]))

  ncols <- 4 * K * L + 2 * L + 2 * K

  X <- matrix(0, nrow = n_theta * n_phi, ncol = ncols)
  col_idx <- 1

  for (l in 1:L) {
    phi_l <- Phi_T[, l]

    X[, col_idx] <-
      as.vector(outer(rep(1, n_theta), cos(phi_l)))

    X[, col_idx + 1] <-
      as.vector(outer(rep(1, n_theta), sin(phi_l)))

    col_idx <- col_idx + 2
  }

  for (k in 1:K) {
    theta_k <- Phi_P[, k]

    X[, col_idx] <-
      as.vector(outer(cos(theta_k), rep(1, n_phi)))

    X[, col_idx + 1] <-
      as.vector(outer(sin(theta_k), rep(1, n_phi)))

    col_idx <- col_idx + 2
  }

  for (k in 1:K) {
    theta_k <- Phi_P[, k]

    for (l in 1:L) {
      phi_l <- Phi_T[, l]

      base1 <- outer(cos(theta_k), cos(phi_l))
      base2 <- outer(sin(theta_k), cos(phi_l))
      base3 <- outer(cos(theta_k), sin(phi_l))
      base4 <- outer(sin(theta_k), sin(phi_l))

      X[, col_idx]     <- as.vector(base1)
      X[, col_idx + 1] <- as.vector(base2)
      X[, col_idx + 2] <- as.vector(base3)
      X[, col_idx + 3] <- as.vector(base4)

      col_idx <- col_idx + 4
    }
  }

  y <- as.vector(f_data - M)

  X_int <- cbind(X, intercept = 1)

  fit <- lm.fit(X_int, y)

  beta_hat <- fit$coefficients
  y_hat <- X_int %*% beta_hat + M

  ss_res <- sum((y - (y_hat - M))^2)
  ss_tot <- sum((y - mean(y))^2)

  r2 <- 1 - ss_res / ss_tot

  f_hat <- matrix(
    y_hat,
    nrow = n_theta,
    ncol = n_phi
  )

  list(
    fitted = f_hat,
    r2 = r2,
    coefficients = beta_hat
  )
}

Ref <- function(
    alpha2, omega2, alpha1, omega1,
    theta, phi, f_data,
    tole = 0.05, it = 5,
    fit_fmm_linear5, phi_mobius, M = 0
) {

  par_init <- c(alpha2, omega2, alpha1, omega1)
  par_current <- par_init
  r2_current <- -Inf
  K <- length(alpha2)
  L <- length(alpha1)

  generate_pairs <- function(K, L) {
    pairs <- list()
    idx <- 1

    for (i in 1:K) {
      pairs[[idx]] <- c(i, i + K)
      idx <- idx + 1
    }

    for (j in 1:L) {
      start_index <- 1 + 2 * K
      pairs[[idx]] <- c(
        start_index + j - 1,
        start_index + j - 1 + L
      )
      idx <- idx + 1
    }

    return(pairs)
  }

  pairs_to_update <- generate_pairs(K, L)

  for (iter in 1:it) {

    for (pair in pairs_to_update) {

      mask <- rep(FALSE, length(par_init))
      mask[pair] <- TRUE

      make_obj_fmm <- function(
          theta, phi, f_data,
          par_init, mask, M = 0
      ) {

        function(par_free) {

          full_par <- par_init
          full_par[mask] <- par_free

          alpha_P <- full_par[1:K]
          omega_P <- full_par[(K + 1):(2 * K)]

          alpha_T <- full_par[
            (2 * K + 1):(2 * K + L)
          ]

          omega_T <- full_par[
            (2 * K + L + 1):(2 * K + 2 * L)
          ]

          resultado <- tryCatch({
  res <- fit_fmm_linear5(theta, phi, f_data,
                         alpha_P, omega_P,
                         alpha_T, omega_T, M)

  val <- -res$r2
  if (!is.finite(val)) 1e10 else val

}, error = function(e) 1e10)

          return(resultado)
        }
      }

      obj <- make_obj_fmm(
        theta,
        phi,
        f_data,
        par_current,
        mask
      )

      par_free_init <- par_current[mask]

      lower <- par_current * (1 - tole)
      upper <- par_current * (1 + tole)

      idx_omega_P <- (K + 1):(2 * K)

      idx_omega_T <- (
        2 * K + L + 1
      ):(
        2 * K + 2 * L
      )

      upper[idx_omega_P] <-
        pmin(1, upper[idx_omega_P])

      upper[idx_omega_T] <-
        pmin(1, upper[idx_omega_T])

      optim_res <- optim(
        par = par_free_init,
        fn = obj,
        method = "L-BFGS-B",
        lower = lower[mask],
        upper = upper[mask]
      )

      par_new <- par_current
      par_new[mask] <- optim_res$par

      par_new[1:K] <-
        par_new[1:K] %% (2 * pi)

      par_new[
        (2 * K + 1):(2 * K + L)
      ] <-
        par_new[
          (2 * K + 1):(2 * K + L)
        ] %% (2 * pi)

      res_new <- fit_fmm_linear5(
        theta,
        phi,
        f_data,
        par_new[1:K],
        par_new[(K + 1):(2 * K)],
        par_new[(2 * K + 1):(2 * K + L)],
        par_new[
          (2 * K + L + 1):(2 * K + 2 * L)
        ],
        M
      )

      r2_new <- res_new$r2

      if (r2_new > r2_current) {
        par_current <- par_new
        r2_current <- r2_new
      }
    }
  }

  n_theta <- length(theta)
  n_phi <- length(phi)

  alpha_P <- par_current[1:K]
  omega_P <- par_current[(K + 1):(2 * K)]

  alpha_T <- par_current[
    (2 * K + 1):(2 * K + L)
  ]

  omega_T <- par_current[
    (2 * K + L + 1):(2 * K + 2 * L)
  ]

  Phi_P <- sapply(
    1:K,
    function(k)
      phi_mobius(
        theta,
        alpha_P[k],
        omega_P[k]
      )
  )

  Phi_T <- sapply(
    1:L,
    function(l)
      phi_mobius(
        phi,
        alpha_T[l],
        omega_T[l]
      )
  )

  X <- matrix(
    0,
    nrow = n_theta * n_phi,
    ncol = 4 * K * L + 2 * L + 2 * K
  )

  col_idx <- 1

  for (l in 1:L) {

    phi_l <- Phi_T[, l]

    X[, col_idx] <-
      as.vector(
        outer(
          rep(1, n_theta),
          cos(phi_l)
        )
      )

    X[, col_idx + 1] <-
      as.vector(
        outer(
          rep(1, n_theta),
          sin(phi_l)
        )
      )

    col_idx <- col_idx + 2
  }

  for (k in 1:K) {

    theta_k <- Phi_P[, k]

    X[, col_idx] <-
      as.vector(
        outer(
          cos(theta_k),
          rep(1, n_phi)
        )
      )

    X[, col_idx + 1] <-
      as.vector(
        outer(
          sin(theta_k),
          rep(1, n_phi)
        )
      )

    col_idx <- col_idx + 2
  }

  for (k in 1:K) {

    for (l in 1:L) {

      phi_k <- Phi_P[, k]
      phi_l <- Phi_T[, l]

      X[, col_idx] <-
        as.vector(
          outer(
            cos(phi_k),
            cos(phi_l)
          )
        )

      X[, col_idx + 1] <-
        as.vector(
          outer(
            sin(phi_k),
            cos(phi_l)
          )
        )

      X[, col_idx + 2] <-
        as.vector(
          outer(
            cos(phi_k),
            sin(phi_l)
          )
        )

      X[, col_idx + 3] <-
        as.vector(
          outer(
            sin(phi_k),
            sin(phi_l)
          )
        )

      col_idx <- col_idx + 4
    }
  }

  X <- cbind(X, intercept = 1)

  y <- as.vector(f_data - M)

  fit <- lm.fit(X, y)

  beta_hat <- fit$coefficients

  y_hat <- X %*% beta_hat + M

  return(
    list(
      par_current = par_current,
      r2 = r2_current,
      beta_hat = beta_hat,
      y_hat = matrix(
        y_hat,
        nrow = n_theta,
        ncol = n_phi
      )
    )
  )
}


S_FMM <- function(K, L, vDataMatrix, path, iter = 30) {
  # Load required libraries
  library(ggplot2)
  library(FMM)
  library(tidyr)
  library(dplyr)
  library(readxl)
  library(akima)
  library(RcppCNPy)
  library(reticulate)
  library(fields)
  library(gplots)
  library(MASS)
  
  # Load auxiliary functions
source(paste0(path, "/requiredFunctionsPreprocessing_v4.1.R"), chdir = TRUE)
source(paste0(path, "/auxMultiFMM.R"), chdir = TRUE)
source(paste0(path, "/FMM_internal.R"), chdir = TRUE)
source(paste0(path, "/auxMultiFMM_inference.R"), chdir = TRUE)
  # Set up data
  ncom <- K
  ncom1 <- L
  Zmat <- vDataMatrix
  dim1 <- nrow(Zmat)
  dim2 <- ncol(Zmat)
  n_theta <- dim1
  n_phi <- dim2
  
  # SVD decomposition
  svd_result <- svd(Zmat)
  U <- svd_result$u
  V <- svd_result$v
  
  # Fit FMM to first singular vectors
  S1 <- fitFMM(vData = U[,1], nback = 3, omegaMin = 0.01, maxiter = 100)
  S2 <- fitFMM(vData = V[,1], nback = 3, omegaMin = 0.01, maxiter = 100)
  
  # Estimation
  ncan <- min(dim1, 40)
  M_time2 <- t(Zmat)
  ind <- seq(1, dim1, length.out = ncan)
  M_time2 <- M_time2[, round(ind)]
  
  fitml <- fitMultiFMM(vDataMatrix = M_time2[,1:ncan], nBack = ncom, maxIter = iter,
                       plotToFile = FALSE, filename = "", omegaMin = 0.01,
                       showPredeterminedPlot = FALSE)
  
  alpha1 <- fitml[[1]]$Alpha
  omega1 <- fitml[[1]]$Omega
  
  # Second stage
  ncoe <- ncom * 2 + 1
  coe <- matrix(0, nrow = dim1, ncol = ncoe)
  theta <- seq(0, 2 * pi, length.out = dim2 + 1)[-(dim2+1)]
  
  Xv <- array(0, c(ncom, 2, dim2))
  for (j in 1:ncom) {
    Xv[j,1,] <- cos(2 * atan2(omega1[j] * sin((theta - alpha1[j]) / 2), cos((theta - alpha1[j]) / 2)))
    Xv[j,2,] <- sin(2 * atan2(omega1[j] * sin((theta - alpha1[j]) / 2), cos((theta - alpha1[j]) / 2)))
  }
  
  Xmat <- matrix(0, nrow = dim2, ncol = ncom * 2)
  for (j in 1:ncom) {
    Xmat[, (j - 1) * 2 + 1] <- Xv[j,1,]
    Xmat[, (j - 1) * 2 + 2] <- Xv[j,2,]
  }
  
  rl2 <- numeric(dim1)
  for (i in 1:dim1) {
    y <- Zmat[i, ]
    ajus <- lm(y ~ Xmat)
    coe[i, ] <- coef(ajus)
    pred <- ajus$fitted.values
    ssres <- sum((Zmat[i,] - pred)^2)
    sstot <- sum((Zmat[i,] - mean(Zmat[i,]))^2)
    rl2[i] <- 1 - ssres / sstot
  }
  
  fitml1 <- fitMultiFMM(vDataMatrix = coe[,1:ncoe], nBack = ncom1, maxIter = iter,
                        plotToFile = FALSE, filename = "", showPredeterminedPlot = FALSE,
                        omegaMin = 0.01)
  
  alpha2 <- fitml1[[1]]$Alpha
  omega2 <- fitml1[[1]]$Omega
  
  # Final fit
  theta <- seq(0, 2 * pi, length.out = dim1)
  phi <- seq(0, 2 * pi, length.out = dim2)
  f_data <- Zmat
  
  final_fit <- fit_fmm_linear5(theta, phi, f_data, alpha2, omega2, alpha1, omega1, M = 0)
  
  # Refinement
  res <- Ref(alpha2, omega2, alpha1, omega1, theta, phi, f_data, tole = 0.9, it = 2,
             fit_fmm_linear5 = fit_fmm_linear5, phi_mobius = phi_mobius)
  par_current <- res$par_current
  K <- length(alpha2)
  L <- length(alpha1)
  alpha2 <- par_current[1:K]
  omega2 <- par_current[(K+1):(2*K)]
  alpha1 <- par_current[(2*K+1):(2*K+L)]
  omega1 <- par_current[(2*K+L+1):(2*K+2*L)]
  
  res <- Ref(alpha2, omega2, alpha1, omega1, theta, phi, f_data, tole = 0.05, it = 2,
             fit_fmm_linear5 = fit_fmm_linear5, phi_mobius = phi_mobius)
  par_current <- res$par_current
  alpha2 <- par_current[1:K]
omega2 <- par_current[(K + 1):(2 * K)]
alpha1 <- par_current[(2 * K + 1):(2 * K + L)]
omega1 <- par_current[(2 * K + L + 1):(2 * K + 2 * L)]
  res <- Ref(alpha2, omega2, alpha1, omega1, theta, phi, f_data, tole = 0.9, it = 2,
             fit_fmm_linear5 = fit_fmm_linear5, phi_mobius = phi_mobius)
  
  len <- ncol(res$X) - 1
  r2 <- res$r2
  beta_hat <- res$beta_hat
  PredFMM <- res$y_hat

  return(list(
    PredFMM = PredFMM,
    r2 = r2,
    alpha1 = alpha1,
    omega1 = omega1,
    alpha2 = alpha2,
    omega2 = omega2,
    beta_hat = beta_hat
  ))
}

S_fmm <- function(M, n_iter, K, L){

  out <- S_FMM(
    K = K,
    L = L,
    vDataMatrix = M,
    path = "R",
  iter = n_iter
  )

  par_current <- c(
    out$alpha2,
    out$omega2,
    out$alpha1,
    out$omega1
  )

  list(
    r2 = out$r2,
    par_current = par_current,
    y_hat = out$PredFMM
  )
}