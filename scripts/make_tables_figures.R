###############################################################
## Tables and figures post-processing
###############################################################

## This script is intentionally minimal. The main simulation scripts already
## write CSV and RDS files to outputs/.

if (file.exists("outputs/table5_protein_toroidal_simulations.csv")) {
  tab5 <- read.csv("outputs/table5_protein_toroidal_simulations.csv")
  print(tab5)
}

if (file.exists("outputs/table7_synthetic_toroidal_simulations.csv")) {
  tab7 <- read.csv("outputs/table7_synthetic_toroidal_simulations.csv")
  print(tab7)
}
