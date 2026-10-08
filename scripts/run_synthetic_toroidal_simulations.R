############################################################
## Structured synthetic toroidal simulations
############################################################
library(FMM)
FMM <- getFromNamespace("FMM", "FMM")
source("R/00_utils_torus.R")
source("R/02_perturbations_torus.R")
source("R/03_fourier2d.R")
source("R/04_vm4_torus.R")
source("R/05_synthetic_torus_generators.R")
source("R/06_toroidal_simulation_core.R")
source("R/S_FMM_adap.R")
if(!exists("S_fmm")) stop("S_fmm() is not available. Source or load the S-FMM implementation before running this script.")

set.seed(123)

B <- 100
n <- 300
kappa <- 25
ng <- 100

base1 <- generate_2modes1(ng)$Z
base2 <- generate_2modes2(ng)$Z

scenarios <- list(
  S1 = list(Z = base1, shift_sd = 0.02, alpha_mix = 0.05, noise_level = 0.01, shear_sd = 0, n_bumps = 0, bump_weight = 0, bump_sd = 0.1),
  S2 = list(Z = base1, shift_sd = 0.03, alpha_mix = 0.07, noise_level = 0.02, shear_sd = 0, n_bumps = 0, bump_weight = 0, bump_sd = 0.1),
  S3 = list(Z = base2, shift_sd = 0.02, alpha_mix = 0.05, noise_level = 0.01, shear_sd = 0, n_bumps = 0, bump_weight = 0, bump_sd = 0.1),
  S4 = list(Z = base2, shift_sd = 0.03, alpha_mix = 0.07, noise_level = 0.02, shear_sd = 0, n_bumps = 0, bump_weight = 0, bump_sd = 0.1)
)

all_results <- list()
summary_rows <- list()
neg_stats <- function(Z){
  
  dA <- (2*pi/nrow(Z)) * (2*pi/ncol(Z))
  
  signed_mass <- sum(Z) * dA
  neg_mass <- -sum(pmin(Z, 0)) * dA
  
  c(
    min_value    = min(Z, na.rm = TRUE),
    prop_neg     = mean(Z < 0, na.rm = TRUE),
    neg_mass     = neg_mass,
    signed_mass  = signed_mass,
    rel_neg_mass = neg_mass / signed_mass
  )
}
for(sc in names(scenarios)){
  cfg <- scenarios[[sc]]
  message("Running synthetic scenario ", sc)

  replicas <- generate_replicates2(
    cfg$Z,
    B = B,
    n = n,
    kappa = kappa,
    ng = ng,
    shift_sd = cfg$shift_sd,
    alpha_mix = cfg$alpha_mix,
    noise_level = cfg$noise_level,
    smooth_iter = 1,
    shear_sd = cfg$shear_sd,
    n_bumps = cfg$n_bumps,
    bump_weight = cfg$bump_weight,
    bump_sd = cfg$bump_sd
  )
  res <- run_toroidal_experiment(cfg$Z, replicas)
  
  ## ============================================================
  ## Additional density diagnostics
  ## ============================================================
  
  BB <- length(res$replicate_results)
  
  ISE_pos_mat <- matrix(
    NA_real_,
    nrow = BB,
    ncol = 3,
    dimnames = list(NULL, c("Fourier33", "S_FMM2", "VM4"))
  )
  
  neg_SFMM_mat <- matrix(
    NA_real_,
    nrow = BB,
    ncol = 5,
    dimnames = list(
      NULL,
      c("min_value", "prop_neg", "neg_mass",
        "signed_mass", "rel_neg_mass")
    )
  )
  
  neg_Fourier_mat <- neg_SFMM_mat
  
  for(i in seq_len(BB)){
    
    Ztrue <- replicas[[i]]$Z_b
    
    ZF <- res$replicate_results[[i]]$fits$fourier$fitted_matrix
    ZS <- res$replicate_results[[i]]$fits$sfmm$y_hat
    
    ## Projection to a valid nonnegative density
    ZF_pos <- as_density_grid(pmax(ZF, 0))
    ZS_pos <- as_density_grid(pmax(ZS, 0))
    
    ## ISE after positivity projection
    ISE_pos_mat[i, "Fourier33"] <-
      ISE_density_grid(Ztrue, ZF_pos)
    
    ISE_pos_mat[i, "S_FMM2"] <-
      ISE_density_grid(Ztrue, ZS_pos)
    
    ## VM4 is already nonnegative, so projected ISE = ordinary ISE
    ISE_pos_mat[i, "VM4"] <-
      res$summary$ISE_mat[i, "VM4"]
    
    ## Negativity of the original unconstrained fits
    neg_SFMM_mat[i, ] <- neg_stats(ZS)
    neg_Fourier_mat[i, ] <- neg_stats(ZF)
  }
  
  ## Save replicate-level diagnostics
  res$extra_diagnostics <- list(
    ISE_pos_mat = ISE_pos_mat,
    neg_SFMM_mat = neg_SFMM_mat,
    neg_Fourier_mat = neg_Fourier_mat
  )
  
  ## Save summaries
  res$summary$mean_ISE_pos <-
    colMeans(ISE_pos_mat, na.rm = TRUE)
  
  res$summary$sd_ISE_pos <-
    apply(ISE_pos_mat, 2, sd, na.rm = TRUE)
  
  res$summary$mean_neg_SFMM <-
    colMeans(neg_SFMM_mat, na.rm = TRUE)
  
  res$summary$sd_neg_SFMM <-
    apply(neg_SFMM_mat, 2, sd, na.rm = TRUE)
  
  res$summary$mean_neg_Fourier <-
    colMeans(neg_Fourier_mat, na.rm = TRUE)
  
  res$summary$sd_neg_Fourier <-
    apply(neg_Fourier_mat, 2, sd, na.rm = TRUE)
  
  res$summary$prop_replicates_negative_Fourier <-
    mean(neg_Fourier_mat[, "min_value"] < 0)
  
  ## NOW save the enriched result
  all_results[[sc]] <- res
  sm <- res$summary
  
  
  
 summary_rows[[sc]] <- data.frame(

  Scenario = sc,

  ## R2 against KDE
  R2_Fourier33_mean = sm$mean_R2["Fourier33"],
  R2_Fourier33_sd   = sm$sd_R2["Fourier33"],

  R2_S_FMM2_mean = sm$mean_R2["S_FMM2"],
  R2_S_FMM2_sd   = sm$sd_R2["S_FMM2"],

  R2_VM4_mean = sm$mean_R2["VM4"],
  R2_VM4_sd   = sm$sd_R2["VM4"],

  ## ISE against true perturbed surface Z_b
  ISE_Fourier33_mean = sm$mean_ISE["Fourier33"],
  ISE_Fourier33_sd   = sm$sd_ISE["Fourier33"],

  ISE_S_FMM2_mean = sm$mean_ISE["S_FMM2"],
  ISE_S_FMM2_sd   = sm$sd_ISE["S_FMM2"],

  ISE_VM4_mean = sm$mean_ISE["VM4"],
  ISE_VM4_sd   = sm$sd_ISE["VM4"],
  ## ISE after truncation at zero + renormalization
  ISEpos_Fourier33_mean = sm$mean_ISE_pos["Fourier33"],
  ISEpos_Fourier33_sd   = sm$sd_ISE_pos["Fourier33"],
  
  ISEpos_S_FMM2_mean = sm$mean_ISE_pos["S_FMM2"],
  ISEpos_S_FMM2_sd   = sm$sd_ISE_pos["S_FMM2"],
  
  ISEpos_VM4_mean = sm$mean_ISE_pos["VM4"],
  ISEpos_VM4_sd   = sm$sd_ISE_pos["VM4"],
  ## Existing structural measures
  CD_S_FMM2 = sm$CD_SFMM,
  CD_VM4 = sm$CD_VM4,
  AD_S_FMM2 = sm$AD_SFMM,
  ## Negativity diagnostics: S-FMM
  Min_S_FMM2 = sm$min_SFMM,
  
  MeanMin_S_FMM2 =
    sm$mean_neg_SFMM["min_value"],
  SDMin_S_FMM2 =
    sm$sd_neg_SFMM["min_value"],
  
  PropNeg_S_FMM2 =
    sm$mean_neg_SFMM["prop_neg"],
  SDPropNeg_S_FMM2 =
    sm$sd_neg_SFMM["prop_neg"],
  
  NegMass_S_FMM2 =
    sm$mean_neg_SFMM["neg_mass"],
  SDNegMass_S_FMM2 =
    sm$sd_neg_SFMM["neg_mass"],
  
  RelNegMass_S_FMM2 =
    sm$mean_neg_SFMM["rel_neg_mass"],
  SDRelNegMass_S_FMM2 =
    sm$sd_neg_SFMM["rel_neg_mass"],
  
  RepNeg_S_FMM2 =
    sm$prop_replicates_negative,
  
  ## Negativity diagnostics: Fourier
  MeanMin_Fourier33 =
    sm$mean_neg_Fourier["min_value"],
  SDMin_Fourier33 =
    sm$sd_neg_Fourier["min_value"],
  
  PropNeg_Fourier33 =
    sm$mean_neg_Fourier["prop_neg"],
  SDPropNeg_Fourier33 =
    sm$sd_neg_Fourier["prop_neg"],
  
  NegMass_Fourier33 =
    sm$mean_neg_Fourier["neg_mass"],
  SDNegMass_Fourier33 =
    sm$sd_neg_Fourier["neg_mass"],
  
  RelNegMass_Fourier33 =
    sm$mean_neg_Fourier["rel_neg_mass"],
  SDRelNegMass_Fourier33 =
    sm$sd_neg_Fourier["rel_neg_mass"],
  
  RepNeg_Fourier33 =
    sm$prop_replicates_negative_Fourier,

  row.names = NULL
)
 
}

summary_table <- do.call(rbind, summary_rows)
print(summary_table)

saveRDS(all_results, file = "outputs/tables/synthetic_toroidal_results.rds")
write.csv(summary_table, file = "outputs/tables/table_synthetic_toroidal.csv", row.names = FALSE)
