############################################################
## Sensitivity to S-FMM model size: S1
## 20 replicates: 1x1 vs 2x2 vs 3x3
############################################################

B_sens <- 20

sens_S1 <- matrix(
  NA_real_,
  nrow = B_sens,
  ncol = 6,
  dimnames = list(
    NULL,
    c(
      "R2_1x1", "ISE_1x1",
      "R2_2x2", "ISE_2x2",
      "R2_3x3", "ISE_3x3"
    )
  )
)

set.seed(123)

for(i in seq_len(B_sens)){
  
  cat("Sensitivity replicate", i, "of", B_sens, "\n")
  
  M     <- replicas[[i]]$M
  Ztrue <- replicas[[i]]$Z_b
  
  ## -------------------------
  ## S-FMM 1x1
  ## -------------------------
  fit_11 <- S_fmm(
    M = M,
    n_iter = 30,
    K = 1,
    L = 1
  )
  
  Z11 <- fit_11$y_hat
  
  sens_S1[i, "R2_1x1"] <-
    R2_2D(M, Z11)
  
  sens_S1[i, "ISE_1x1"] <-
    ISE_density_grid(Ztrue, Z11)
  
  
  ## -------------------------
  ## S-FMM 2x2
  ## already computed
  ## -------------------------
  Z22 <- res$replicate_results[[i]]$fits$sfmm$y_hat
  
  sens_S1[i, "R2_2x2"] <-
    R2_2D(M, Z22)
  
  sens_S1[i, "ISE_2x2"] <-
    ISE_density_grid(Ztrue, Z22)
  
  
  ## -------------------------
  ## S-FMM 3x3
  ## -------------------------
  fit_33 <- S_fmm(
    M = M,
    n_iter = 30,
    K = 3,
    L = 3
  )
  
  Z33 <- fit_33$y_hat
  
  sens_S1[i, "R2_3x3"] <-
    R2_2D(M, Z33)
  
  sens_S1[i, "ISE_3x3"] <-
    ISE_density_grid(Ztrue, Z33)
}


## Summary: mean and SD
sens_S1_summary <- rbind(
  Mean = colMeans(sens_S1, na.rm = TRUE),
  SD   = apply(sens_S1, 2, sd, na.rm = TRUE)
)

print(sens_S1_summary)

## Save
saveRDS(
  sens_S1,
  file = "outputs/tables/sensitivity_S1_1x1_2x2_3x3.rds"
)

write.csv(
  sens_S1_summary,
  file = "outputs/tables/sensitivity_S1_1x1_2x2_3x3.csv"
)
sens_S1_summary