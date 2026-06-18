###############################################################
## Protein-based toroidal simulations: 1CAG and 1V9E
###############################################################

source("R/toroidal_utils.R")
source("R/toroidal_pdb.R")
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

scenarios <- list(
  S1 = list(pdb = "1CAG", type = "noise", noise_level = 0.02, shift_sd = 3, alpha_mix = 0.12, shear_sd = NA),
  S2 = list(pdb = "1CAG", type = "noise", noise_level = 0.50, shift_sd = 3, alpha_mix = 0.12, shear_sd = NA),
  S3 = list(pdb = "1CAG", type = "structured", noise_level = 0.01, shift_sd = 0.15, alpha_mix = 0.10, shear_sd = 0.08),
  S4 = list(pdb = "1V9E", type = "noise", noise_level = 0.02, shift_sd = 3, alpha_mix = 0.12, shear_sd = NA),
  S5 = list(pdb = "1V9E", type = "noise", noise_level = 0.50, shift_sd = 3, alpha_mix = 0.12, shear_sd = NA),
  S6 = list(pdb = "1V9E", type = "structured", noise_level = 0.01, shift_sd = 0.15, alpha_mix = 0.10, shear_sd = 0.08)
)

results <- list()
summary_rows <- list()

for (sc in names(scenarios)) {
  cfg <- scenarios[[sc]]
  cat("Running", sc, cfg$pdb, "\n")
  Z <- get_Z_from_pdb(cfg$pdb, kappa = kappa, ng = ng)
  if (cfg$type == "noise") {
    reps <- generate_replicates(Z, B = B, n = n, kappa = kappa, ng = ng,
                                shift_sd = cfg$shift_sd, alpha_mix = cfg$alpha_mix,
                                smooth_iter = 1, noise_level = cfg$noise_level)
  } else {
    reps <- generate_replicates2(Z, B = B, n = n, kappa = kappa, ng = ng,
                                 shift_sd = cfg$shift_sd, alpha_mix = cfg$alpha_mix,
                                 smooth_iter = 1, noise_level = cfg$noise_level,
                                 shear_sd = cfg$shear_sd, n_bumps = 1,
                                 bump_weight = 0.15, bump_sd = 0.30)
  }
  res <- summarize_toroidal_replicates(Z, reps, N_vm = 3000, seed = 1)
  results[[sc]] <- res
  summary_rows[[sc]] <- one_row_summary(sc, res)
}

summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, "outputs/table5_protein_toroidal_simulations.csv", row.names = FALSE)
saveRDS(results, "outputs/protein_toroidal_simulations.rds")
print(summary_table)
