## ============================================================
## Comparison of Möbius-warped trigonometric circular density
## with classical circular densities
## ============================================================

setwd("C:/Users/cristina/Desktop/FMM/MObiusDENsity/")

## ============================================================
## Möbius-warped trigonometric circular density
## Simulation, fitting, AIC comparison and plots
## ============================================================

rm(list = ls())

## ============================================================
## 1. DENSITIES
## ============================================================

Phi_mobius <- function(theta, alpha, omega) {
  u <- theta - alpha
  2 * atan2(omega * sin(u / 2), cos(u / 2))
}

d_mobius_trig <- function(theta, alpha, omega, beta, rho, log = FALSE) {
  
  num <- 1 + rho * cos(Phi_mobius(theta, alpha, omega) + beta)
  
  Z <- 2*pi * (1 + rho * cos(beta) * (1 - omega) / (1 + omega))
  
  dens <- num / Z
  
  if (any(dens <= 0) || Z <= 0 || any(!is.finite(dens))) {
    if (log) return(rep(-Inf, length(theta)))
    return(rep(0, length(theta)))
  }
  
  if (log) log(dens) else dens
}

d_cardioid <- function(theta, mu = 0, rho = 0.3) {
  (1 + 2*rho*cos(theta - mu)) / (2*pi)
}

d_vonmises <- function(theta, mu = 0, kappa = 2) {
  exp(kappa*cos(theta - mu)) / (2*pi*besselI(kappa, 0))
}

d_wrapped_cauchy <- function(theta, mu = 0, gamma = 0.5) {
  (1 - gamma^2) /
    (2*pi * (1 + gamma^2 - 2*gamma*cos(theta - mu)))
}

d_sineskew_vm <- function(theta, mu = 0, kappa = 2, lambda = 0.3) {
  base <- d_vonmises(theta, mu, kappa)
  base * (1 + lambda * sin(theta - mu))
}

d_ext_wrapped_cauchy <- function(theta, mu = 0, gamma = 0.5, lambda = 0.3) {
  base <- d_wrapped_cauchy(theta, mu, gamma)
  base * (1 + lambda * sin(theta - mu))
}

d_gen_cardioid2 <- function(theta,
                            mu = 0,
                            a1 = 0.3,
                            b1 = 0,
                            a2 = 0.1,
                            b2 = 0) {
  
  x <- theta - mu
  
  (1 +
      a1*cos(x) + b1*sin(x) +
      a2*cos(2*x) + b2*sin(2*x)) / (2*pi)
}

## ============================================================
## 2. NEGATIVE LOG-LIKELIHOODS
## ============================================================

neglog_mobius <- function(par, theta) {
  
  alpha <- par[1] %% (2*pi)
  omega <- plogis(par[2])
  beta  <- par[3] %% (2*pi)
  rho   <- tanh(par[4])
  
  ll <- d_mobius_trig(theta, alpha, omega, beta, rho, log = TRUE)
  
  if (any(!is.finite(ll))) return(1e10)
  
  -sum(ll)
}

neglog_cardioid <- function(par, theta) {
  
  mu <- par[1] %% (2*pi)
  rho <- 0.5 * tanh(par[2])
  
  f <- d_cardioid(theta, mu, rho)
  
  if (any(f <= 0) || any(!is.finite(f))) return(1e10)
  
  -sum(log(f))
}

neglog_vm <- function(par, theta) {
  
  mu <- par[1] %% (2*pi)
  kappa <- exp(par[2])
  
  f <- d_vonmises(theta, mu, kappa)
  
  if (any(f <= 0) || any(!is.finite(f))) return(1e10)
  
  -sum(log(f))
}

neglog_wc <- function(par, theta) {
  
  mu <- par[1] %% (2*pi)
  gamma <- plogis(par[2])
  
  f <- d_wrapped_cauchy(theta, mu, gamma)
  
  if (any(f <= 0) || any(!is.finite(f))) return(1e10)
  
  -sum(log(f))
}

neglog_sineskew_vm <- function(par, theta) {
  
  mu <- par[1] %% (2*pi)
  kappa <- exp(par[2])
  lambda <- tanh(par[3])
  
  f <- d_sineskew_vm(theta, mu, kappa, lambda)
  
  if (any(f <= 0) || any(!is.finite(f))) return(1e10)
  
  -sum(log(f))
}

neglog_ext_wc <- function(par, theta) {
  
  mu <- par[1] %% (2*pi)
  gamma <- plogis(par[2])
  lambda <- tanh(par[3])
  
  f <- d_ext_wrapped_cauchy(theta, mu, gamma, lambda)
  
  if (any(f <= 0) || any(!is.finite(f))) return(1e10)
  
  -sum(log(f))
}

neglog_gen_cardioid2 <- function(par, theta) {
  
  mu <- par[1] %% (2*pi)
  
  a1 <- tanh(par[2])
  b1 <- tanh(par[3])
  a2 <- tanh(par[4])
  b2 <- tanh(par[5])
  
  f <- d_gen_cardioid2(theta, mu, a1, b1, a2, b2)
  
  if (any(f <= 0) || any(!is.finite(f))) return(1e10)
  
  -sum(log(f))
}

## ============================================================
## 3. FITTING FUNCTIONS
## ============================================================

fit_generic <- function(theta, fn, npar, nstart = 20) {
  
  theta <- theta %% (2*pi)
  
  starts <- replicate(
    nstart,
    c(runif(1, 0, 2*pi), rnorm(npar - 1)),
    simplify = FALSE
  )
  
  fits <- lapply(starts, function(st) {
    
    optim(
      par = st,
      fn = fn,
      theta = theta,
      method = "BFGS",
      control = list(maxit = 2000)
    )
    
  })
  
  fits[[which.min(sapply(fits, function(z) z$value))]]
}

fit_mobius_trig <- function(theta, nstart = 30) {
  
  theta <- theta %% (2*pi)
  
  starts <- replicate(nstart, {
    c(
      runif(1, 0, 2*pi),
      rnorm(1, 0, 1),
      runif(1, 0, 2*pi),
      rnorm(1, 0, 1)
    )
  }, simplify = FALSE)
  
  fits <- lapply(starts, function(st) {
    
    optim(
      par = st,
      fn = neglog_mobius,
      theta = theta,
      method = "BFGS",
      control = list(maxit = 2000)
    )
    
  })
  
  best <- fits[[which.min(sapply(fits, function(z) z$value))]]
  
  est <- c(
    alpha = best$par[1] %% (2*pi),
    omega = plogis(best$par[2]),
    beta  = best$par[3] %% (2*pi),
    rho   = tanh(best$par[4])
  )
  
  list(
    par = est,
    logLik = -best$value,
    optim = best
  )
}

## ============================================================
## 4. RANDOM GENERATION
## ============================================================

r_from_density_reject <- function(n, dfun, ...) {
  
  grid <- seq(0, 2*pi, length.out = 5000)
  
  fmax <- max(dfun(grid, ...))
  
  out <- numeric(0)
  
  while (length(out) < n) {
    
    cand <- runif(n, 0, 2*pi)
    u <- runif(n)
    
    f <- dfun(cand, ...)
    
    acc <- u < f/fmax
    
    out <- c(out, cand[acc])
  }
  
  out[1:n]
}

r_vonmises_reject <- function(n, mu = 0, kappa = 2) {
  r_from_density_reject(n, d_vonmises, mu = mu, kappa = kappa)
}

## ============================================================
## 5. SIMULATE ONE SAMPLE FROM EACH TRUTH
## ============================================================

set.seed(123)

n <- 800

theta_vm <- r_vonmises_reject(
  n = n,
  mu = 2.4,
  kappa = 4
)

theta_wc <- r_from_density_reject(
  n = n,
  dfun = d_wrapped_cauchy,
  mu = 2.4,
  gamma = 0.65
)

theta_ssvm <- r_from_density_reject(
  n = n,
  dfun = d_sineskew_vm,
  mu = 2.4,
  kappa = 2.2,
  lambda = 0.90
)

theta_ewc <- r_from_density_reject(
  n = n,
  dfun = d_ext_wrapped_cauchy,
  mu = 2.4,
  gamma = 0.75,
  lambda = 0.90
)

theta_gc2 <- r_from_density_reject(
  n = n,
  dfun = d_gen_cardioid2,
  mu = 2.4,
  a1 = 0.85,
  b1 = 0.55,
  a2 = 0.25,
  b2 = -0.18
)
## ============================================================
## 6. MODEL COMPARISON
## ============================================================

simulate_and_compare <- function(theta, label = "truth",
                                 nstart_main = 30,
                                 nstart_gc2 = 50) {
  
  fit_mob <- fit_mobius_trig(theta, nstart = nstart_main)
  
  fit_car <- fit_generic(theta, neglog_cardioid, npar = 2, nstart = nstart_main)
  fit_vm  <- fit_generic(theta, neglog_vm,       npar = 2, nstart = nstart_main)
  fit_wc  <- fit_generic(theta, neglog_wc,       npar = 2, nstart = nstart_main)
  
  fit_ssvm <- fit_generic(theta, neglog_sineskew_vm, npar = 3, nstart = nstart_main)
  fit_ewc  <- fit_generic(theta, neglog_ext_wc,      npar = 3, nstart = nstart_main)
  fit_gc2  <- fit_generic(theta, neglog_gen_cardioid2,
                          npar = 5,
                          nstart = nstart_gc2)
  
  out <- data.frame(
    Truth = label,
    Model = c(
      "Mobius-trig",
      "Cardioid",
      "von Mises",
      "Wrapped Cauchy",
      "Sine-skewed von Mises",
      "Extended wrapped Cauchy",
      "Generalized cardioid 2"
    ),
    Parameters = c(4, 2, 2, 2, 3, 3, 5),
    logLik = c(
      fit_mob$logLik,
      -fit_car$value,
      -fit_vm$value,
      -fit_wc$value,
      -fit_ssvm$value,
      -fit_ewc$value,
      -fit_gc2$value
    )
  )
  
  out$AIC <- -2*out$logLik + 2*out$Parameters
  out$deltaAIC <- out$AIC - min(out$AIC)
  
  out[order(out$AIC), ]
}

res_all <- rbind(
  simulate_and_compare(theta_vm,   "von Mises"),
  simulate_and_compare(theta_wc,   "Wrapped Cauchy"),
  simulate_and_compare(theta_ssvm, "Sine-skewed von Mises"),
  simulate_and_compare(theta_ewc,  "Extended wrapped Cauchy"),
  simulate_and_compare(theta_gc2,  "Generalized cardioid 2")
)

print(res_all)

winners <- do.call(rbind, lapply(split(res_all, res_all$Truth), function(z) {
  z[which.min(z$AIC), ]
}))

print(winners)

delta_table <- reshape(
  res_all[, c("Truth", "Model", "deltaAIC")],
  idvar = "Truth",
  timevar = "Model",
  direction = "wide"
)


## ============================================================
## 7. PLOT FUNCTION
## ============================================================

plot_all_fits <- function(theta,
                          truth_fun,
                          truth_args,
                          plot_title = "Simulation",
                          nstart_main = 15,
                          nstart_gc2 = 25) {
  
  grid <- seq(0, 2*pi, length.out = 1200)
  
  f_true <- do.call(truth_fun, c(list(theta = grid), truth_args))
  
  fit_mob <- fit_mobius_trig(theta, nstart = nstart_main)
  
  fit_car <- fit_generic(theta, neglog_cardioid, npar = 2, nstart = nstart_main)
  fit_vm  <- fit_generic(theta, neglog_vm,       npar = 2, nstart = nstart_main)
  fit_wc  <- fit_generic(theta, neglog_wc,       npar = 2, nstart = nstart_main)
  
  fit_ssvm <- fit_generic(theta, neglog_sineskew_vm, npar = 3, nstart = nstart_main)
  fit_ewc  <- fit_generic(theta, neglog_ext_wc,      npar = 3, nstart = nstart_main)
  fit_gc2  <- fit_generic(theta, neglog_gen_cardioid2,
                          npar = 5,
                          nstart = nstart_gc2)
  
  par_car <- c(
    mu = fit_car$par[1] %% (2*pi),
    rho = 0.5 * tanh(fit_car$par[2])
  )
  
  par_vm <- c(
    mu = fit_vm$par[1] %% (2*pi),
    kappa = exp(fit_vm$par[2])
  )
  
  par_wc <- c(
    mu = fit_wc$par[1] %% (2*pi),
    gamma = plogis(fit_wc$par[2])
  )
  
  par_ssvm <- c(
    mu = fit_ssvm$par[1] %% (2*pi),
    kappa = exp(fit_ssvm$par[2]),
    lambda = tanh(fit_ssvm$par[3])
  )
  
  par_ewc <- c(
    mu = fit_ewc$par[1] %% (2*pi),
    gamma = plogis(fit_ewc$par[2]),
    lambda = tanh(fit_ewc$par[3])
  )
  
  par_gc2 <- c(
    mu = fit_gc2$par[1] %% (2*pi),
    a1 = tanh(fit_gc2$par[2]),
    b1 = tanh(fit_gc2$par[3]),
    a2 = tanh(fit_gc2$par[4]),
    b2 = tanh(fit_gc2$par[5])
  )
  
  f_mob <- d_mobius_trig(
    grid,
    fit_mob$par["alpha"],
    fit_mob$par["omega"],
    fit_mob$par["beta"],
    fit_mob$par["rho"]
  )
  
  f_car <- d_cardioid(grid, par_car["mu"], par_car["rho"])
  
  f_vm <- d_vonmises(grid, par_vm["mu"], par_vm["kappa"])
  
  f_wc <- d_wrapped_cauchy(grid, par_wc["mu"], par_wc["gamma"])
  
  f_ssvm <- d_sineskew_vm(
    grid,
    par_ssvm["mu"],
    par_ssvm["kappa"],
    par_ssvm["lambda"]
  )
  
  f_ewc <- d_ext_wrapped_cauchy(
    grid,
    par_ewc["mu"],
    par_ewc["gamma"],
    par_ewc["lambda"]
  )
  
  f_gc2 <- d_gen_cardioid2(
    grid,
    par_gc2["mu"],
    par_gc2["a1"],
    par_gc2["b1"],
    par_gc2["a2"],
    par_gc2["b2"]
  )
  
  ymax <- max(c(f_true, f_mob, f_car, f_vm, f_wc, f_ssvm, f_ewc, f_gc2))
  
  hist(theta,
       breaks = 45,
       probability = TRUE,
       xlim = c(0, 2*pi),
       ylim = c(0, 1.10*ymax),
       main = plot_title,
       xlab = expression(theta),
       ylab = "Density",
       border = "white",
       col = "grey90",
       cex.main = 1.2,
       cex.lab = 1.1)
  
  lines(grid, f_true, lwd = 4, col = "black")
  lines(grid, f_mob,  lwd = 3, col = "red")
  lines(grid, f_vm,   lwd = 2.5, col = "blue")
  lines(grid, f_wc,   lwd = 2.5, col = "darkgreen")
  lines(grid, f_ssvm, lwd = 2.5, col = "orange")
  lines(grid, f_ewc,  lwd = 2.5, col = "purple")
  lines(grid, f_gc2,  lwd = 2.5, col = "brown")
  
  legend("topright",
         legend = c("True density",
                    "Mobius-trig",
                    "von Mises",
                    "Wrapped Cauchy",
                    "Sine-skewed vM",
                    "Extended WC",
                    "Gen. cardioid 2"),
         col = c("black", "red", "blue", "darkgreen",
                 "orange", "purple", "brown"),
         lwd = c(4, 3, 2.5, 2.5, 2.5, 2.5, 2.5),
         cex = 0.75,
         bty = "n")
}

## ============================================================
## 8. GENERATE PLOTS ONE BY ONE
## ============================================================

par(mfrow = c(1,1))

plot_all_fits(
  theta_vm,
  truth_fun = d_vonmises,
  truth_args = list(mu = 2.4, kappa = 4),
  plot_title = "Truth: von Mises"
)

plot_all_fits(
  theta_wc,
  truth_fun = d_wrapped_cauchy,
  truth_args = list(mu = 2.4, gamma = 0.65),
  plot_title = "Truth: Wrapped Cauchy"
)

plot_all_fits(
  theta_ssvm,
  truth_fun = d_sineskew_vm,
  truth_args = list(mu = 2.4, kappa = 2.2, lambda = 0.90),
  plot_title = "Truth: strongly sine-skewed von Mises"
)

plot_all_fits(
  theta_ewc,
  truth_fun = d_ext_wrapped_cauchy,
  truth_args = list(mu = 2.4, gamma = 0.75, lambda = 0.90),
  plot_title = "Truth: strongly extended wrapped Cauchy"
)

plot_all_fits(
  theta_gc2,
  truth_fun = d_gen_cardioid2,
  truth_args = list(mu = 2.4, a1 = 0.85, b1 = 0.55,
                    a2 = 0.25, b2 = -0.18),
  plot_title = "Truth: asymmetric generalized cardioid 2"
)

## ============================================================
## Simulate from true aFMM density
## ============================================================

theta_afmm <- r_from_density_reject(
  n = n,
  dfun = d_mobius_trig,
  alpha = 2.4,
  omega = 0.35,
  beta  = 1.1,
  rho   = 0.80
)

res_afmm <- simulate_and_compare(
  theta_afmm,
  "aFMM"
)

res_all2 <- rbind(res_all, res_afmm)

delta_table2 <- reshape(
  res_all2[, c("Truth", "Model", "deltaAIC")],
  idvar = "Truth",
  timevar = "Model",
  direction = "wide"
)
## ============================================================
## Plot: true aFMM density
## ============================================================

plot_all_fits(
  theta_afmm,
  truth_fun = d_mobius_trig,
  truth_args = list(
    alpha = 2.4,
    omega = 0.35,
    beta  = 1.1,
    rho   = 0.80
  ),
  plot_title = "Truth: aFMM"
)
print(delta_table2)
print(delta_table)