###############################################################
## Synthetic toroidal simulations
###############################################################

source("R/toroidal_utils.R")
source("R/toroidal_perturbations.R")
source("R/toroidal_fourier33.R")
source("R/toroidal_vm4.R")
source("R/toroidal_simulation_pipeline.R")

## IMPORTANT:
## Before running, source the core S-FMM code that defines S_fmm().
## Example:
## source("path/to/S_fmm_core.R")

set.seed(123)
B <- 100
n <- 300
ng <- 100
kappa <- 25

Z1 <- generate_2modes1(ng)$Z
Z2 <- generate_2modes2(ng)$Z

scenarios <- list(
  S1 = list(Z = Z1, shift_sd = 0.02, alpha_mix = 0.05, noise_level = 0.01),
  S2 = list(Z = Z1, shift_sd = 0.03, alpha_mix = 0.07, noise_level = 0.02),
  S3 = list(Z = Z2, shift_sd = 0.02, alpha_mix = 0.05, noise_level = 0.01),
  S4 = list(Z = Z2, shift_sd = 0.03, alpha_mix = 0.07, noise_level = 0.02)
)

results <- list()
summary_rows <- list()

for (sc in names(scenarios)) {
  cfg <- scenarios[[sc]]
  cat("Running", sc, "\n")
  reps <- generate_replicates2(cfg$Z, B = B, n = n, kappa = kappa, ng = ng,
                               shift_sd = cfg$shift_sd, alpha_mix = cfg$alpha_mix,
                               smooth_iter = 1, noise_level = cfg$noise_level,
                               shear_sd = 0, n_bumps = 0, bump_weight = 0, bump_sd = 0.1)
  res <- summarize_toroidal_replicates(cfg$Z, reps, N_vm = 3000, seed = 1)
  results[[sc]] <- res
  summary_rows[[sc]] <- one_row_summary(sc, res)
}

summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, "outputs/table7_synthetic_toroidal_simulations.csv", row.names = FALSE)
saveRDS(results, "outputs/synthetic_toroidal_simulations.rds")
print(summary_table)
