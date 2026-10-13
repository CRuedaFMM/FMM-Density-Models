############################################################
## Protein-based toroidal simulations: 1CAG and 1V9E
############################################################
library(FMM)
FMM <- getFromNamespace("FMM", "FMM")

source("R/S_FMM_adap.R")
source("R/00_utils_torus.R")
source("R/01_density_from_pdb.R")
source("R/02_perturbations_torus.R")
source("R/03_fourier2d.R")
source("R/04_vm4_torus.R")
source("R/06_toroidal_simulation_core.R")

## S_fmm must be available from the external S-FMM package/repository.
if(!exists("S_fmm")) stop("S_fmm() is not available. Source or load the S-FMM implementation before running this script.")

set.seed(123)

ng <- 100
kappa <- 25
B <- 100
n <- 300

scenarios <- list(
  S1 = list(pdb = "1CAG", perturbation = "noise", shift_sd = 0.02, alpha_mix = 0.05, noise_level = 0.02, smooth_iter = 1, shear_sd = 0, n_bumps = 0, bump_weight = 0, bump_sd = 0.1),
  S2 = list(pdb = "1CAG", perturbation = "noise", shift_sd = 0.02, alpha_mix = 0.05, noise_level = 0.50, smooth_iter = 1, shear_sd = 0, n_bumps = 0, bump_weight = 0, bump_sd = 0.1),
  S3 = list(pdb = "1CAG", perturbation = "structured", shift_sd = 0.02, alpha_mix = 0.05, noise_level = 0.01, smooth_iter = 1, shear_sd = 0.08, n_bumps = 2, bump_weight = 0.15, bump_sd = 0.30),
  S4 = list(pdb = "1V9E", perturbation = "noise", shift_sd = 0.02, alpha_mix = 0.05, noise_level = 0.02, smooth_iter = 1, shear_sd = 0, n_bumps = 0, bump_weight = 0, bump_sd = 0.1),
  S5 = list(pdb = "1V9E", perturbation = "noise", shift_sd = 0.02, alpha_mix = 0.05, noise_level = 0.50, smooth_iter = 1, shear_sd = 0, n_bumps = 0, bump_weight = 0, bump_sd = 0.1),
  S6 = list(pdb = "1V9E", perturbation = "structured", shift_sd = 0.02, alpha_mix = 0.05, noise_level = 0.01, smooth_iter = 1, shear_sd = 0.08, n_bumps = 2, bump_weight = 0.15, bump_sd = 0.30)
)

all_results <- list()
summary_rows <- list()

for(sc in names(scenarios)){
  cfg <- scenarios[[sc]]
  message("Running ", sc, " for ", cfg$pdb)

  Z <- get_Z_from_pdbs(pdb_ids = cfg$pdb, kappa = kappa, ng = ng)

  replicas <- generate_replicates2(
    Z,
    B = B,
    n = n,
    kappa = kappa,
    ng = ng,
    shift_sd = cfg$shift_sd,
    alpha_mix = cfg$alpha_mix,
    noise_level = cfg$noise_level,
    smooth_iter = cfg$smooth_iter,
    shear_sd = cfg$shear_sd,
    n_bumps = cfg$n_bumps,
    bump_weight = cfg$bump_weight,
    bump_sd = cfg$bump_sd
  )

  res <- run_toroidal_experiment(Z, replicas)
  all_results[[sc]] <- res

  sm <- res$summary

summary_rows[[sc]] <- data.frame(
  Scenario = sc,

  R2_Fourier33_mean = sm$mean_R2["Fourier33"],
  R2_Fourier33_sd   = sm$sd_R2["Fourier33"],

  R2_S_FMM2_mean = sm$mean_R2["S_FMM2"],
  R2_S_FMM2_sd   = sm$sd_R2["S_FMM2"],

  R2_VM4_mean = sm$mean_R2["VM4"],
  R2_VM4_sd   = sm$sd_R2["VM4"],

  CD_S_FMM2 = sm$CD_SFMM,
  CD_VM4 = sm$CD_VM4,
  AD_S_FMM2 = sm$AD_SFMM,
ISE_Fourier33_mean = sm$mean_ISE["Fourier33"],
ISE_Fourier33_sd   = sm$sd_ISE["Fourier33"],

ISE_S_FMM2_mean = sm$mean_ISE["S_FMM2"],
ISE_S_FMM2_sd   = sm$sd_ISE["S_FMM2"],

ISE_VM4_mean = sm$mean_ISE["VM4"],
ISE_VM4_sd   = sm$sd_ISE["VM4"],
  row.names = NULL
)


## Save each scenario immediately after completion
saveRDS(
  res,
  file = paste0("outputs/tables/protein_toroidal_results_", sc, ".rds")
)

write.csv(
  summary_rows[[sc]],
  file = paste0("outputs/tables/table_protein_toroidal_", sc, ".csv"),
  row.names = FALSE
)

}

summary_table <- do.call(rbind, summary_rows)
print(summary_table)

saveRDS(all_results, file = "outputs/tables/protein_toroidal_results.rds")
write.csv(summary_table, file = "outputs/tables/table_protein_toroidal.csv", row.names = FALSE)
