###############################################################
## Ramachandran KDE surfaces from PDB structures
###############################################################

get_Z_from_pdb <- function(pdb_file, kappa = 25, ng = 100) {
  if (!requireNamespace("bio3d", quietly = TRUE)) stop("Package 'bio3d' is required.")
  if (!requireNamespace("ClusTorus", quietly = TRUE)) stop("Package 'ClusTorus' is required.")
  pdb <- bio3d::read.pdb(pdb_file)
  tor <- bio3d::torsion.pdb(pdb)
  dat <- data.frame(phi_deg = tor$phi, psi_deg = tor$psi)
  dat <- dat[is.finite(dat$phi_deg) & is.finite(dat$psi_deg), ]
  X <- cbind((dat$phi_deg*pi/180) %% (2*pi), (dat$psi_deg*pi/180) %% (2*pi))
  G <- grid_torus(ng)
  fhat <- ClusTorus::kde.torus(data = X, eval.point = G, concentration = kappa)
  normalize_mat(matrix(fhat, nrow = ng, ncol = ng))
}

get_Z_from_pdbs <- function(pdb_ids, kappa = 25, ng = 100,
                            groups = c("G", "X"), y_min = -Inf, y_max = Inf) {
  if (!requireNamespace("bio3d", quietly = TRUE)) stop("Package 'bio3d' is required.")
  if (!requireNamespace("ClusTorus", quietly = TRUE)) stop("Package 'ClusTorus' is required.")
  all_rows <- list(); kk <- 1
  for (id in pdb_ids) {
    cat("Reading", id, "...\n")
    pdb <- tryCatch(bio3d::read.pdb(id), error = function(e) NULL)
    if (is.null(pdb)) next
    tors <- bio3d::torsion.pdb(pdb)
    atoms <- pdb$atom
    xyz <- matrix(pdb$xyz, ncol = 3, byrow = TRUE)
    atoms$insert2 <- ifelse(is.na(atoms$insert), "", atoms$insert)
    ca_idx <- which(atoms$elety == "CA")
    if (length(ca_idx) == 0) next
    ca_tab <- atoms[ca_idx, c("chain", "resno", "insert2", "resid")]
    ca_tab$resid <- toupper(trimws(as.character(ca_tab$resid)))
    nres <- min(nrow(ca_tab), length(tors$phi), length(tors$psi))
    if (nres == 0) next
    centroid <- colMeans(xyz, na.rm = TRUE)
    ca_xyz <- xyz[ca_idx, , drop = FALSE]
    ca_dist <- sqrt(rowSums((ca_xyz - matrix(centroid, nrow = nrow(ca_xyz), ncol = 3, byrow = TRUE))^2))
    max_dist <- max(ca_dist, na.rm = TRUE)
    for (i in seq_len(nres)) {
      chain_i <- ca_tab$chain[i]; resno_i <- ca_tab$resno[i]
      insert_i <- ca_tab$insert2[i]; resid_i <- ca_tab$resid[i]
      sel <- which(atoms$chain == chain_i & atoms$resno == resno_i & atoms$insert2 == insert_i)
      N_i <- sel[atoms$elety[sel] == "N"][1]
      CA_i <- sel[atoms$elety[sel] == "CA"][1]
      C_i <- sel[atoms$elety[sel] == "C"][1]
      if (!is.na(N_i) && !is.na(CA_i) && !is.na(C_i) &&
          !is.na(tors$phi[i]) && !is.na(tors$psi[i])) {
        dist_to_centroid <- sqrt(sum((xyz[CA_i, ] - centroid)^2))
        y <- max_dist - dist_to_centroid
        grp <- NA_character_
        if (resid_i == "GLY") grp <- "G"
        if (!(resid_i %in% c("GLY", "PRO", "ILE", "VAL"))) grp <- "X"
        if (!is.na(grp)) {
          all_rows[[kk]] <- data.frame(group = grp, phi = tors$phi[i], psi = tors$psi[i], y = y)
          kk <- kk + 1
        }
      }
    }
  }
  if (length(all_rows) == 0) stop("No valid residues extracted.")
  dat <- do.call(rbind, all_rows)
  dat <- dat[dat$group %in% groups & is.finite(dat$phi) & is.finite(dat$psi) & dat$y >= y_min & dat$y <= y_max, ]
  if (nrow(dat) == 0) stop("No data after filtering.")
  X <- cbind((dat$phi*pi/180) %% (2*pi), (dat$psi*pi/180) %% (2*pi))
  G <- grid_torus(ng)
  fhat <- ClusTorus::kde.torus(data = X, eval.point = G, concentration = kappa)
  normalize_mat(matrix(fhat, nrow = ng, ncol = ng))
}
