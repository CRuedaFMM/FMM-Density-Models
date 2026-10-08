############################################################
## Ramachandran densities from PDB structures
############################################################

get_Z_from_pdb <- function(pdb_file, kappa = 25, ng = 100){
  requireNamespace("bio3d")
  requireNamespace("dplyr")
  requireNamespace("tibble")
  requireNamespace("ClusTorus")

  pdb <- bio3d::read.pdb(pdb_file)
  tor <- bio3d::torsion.pdb(pdb)

  dat <- tibble::tibble(phi_deg = tor$phi, psi_deg = tor$psi) |>
    dplyr::filter(!is.na(phi_deg), !is.na(psi_deg)) |>
    dplyr::mutate(phi = (phi_deg*pi/180) %% (2*pi),
                  psi = (psi_deg*pi/180) %% (2*pi)) |>
    dplyr::select(phi, psi)

  G <- grid_torus(ng)
  fhat <- ClusTorus::kde.torus(data = as.matrix(dat), eval.point = G, concentration = kappa)
  normalize_mat(matrix(fhat, nrow = ng, ncol = ng))
}

get_Z_from_pdbs <- function(pdb_ids,
                            kappa = 25,
                            ng = 100,
                            groups = c("G", "X"),
                            y_min = -Inf,
                            y_max = Inf){
  requireNamespace("bio3d")
  requireNamespace("dplyr")
  requireNamespace("ClusTorus")

  all_rows <- list(); kk <- 1

  for(id in pdb_ids){
    message("Reading ", id, " ...")
    pdb <- tryCatch(bio3d::read.pdb(id), error = function(e) NULL)
    if(is.null(pdb)) next

    tors <- bio3d::torsion.pdb(pdb)
    atoms <- pdb$atom
    xyz <- matrix(pdb$xyz, ncol = 3, byrow = TRUE)
    atoms$insert2 <- ifelse(is.na(atoms$insert), "", atoms$insert)

    ca_idx <- which(atoms$elety == "CA")
    if(length(ca_idx) == 0) next

    ca_tab <- atoms[ca_idx, c("chain", "resno", "insert2", "resid")]
    ca_tab$resid <- toupper(trimws(as.character(ca_tab$resid)))

    nres <- min(nrow(ca_tab), length(tors$phi), length(tors$psi))
    if(nres == 0) next

    centroid <- colMeans(xyz, na.rm = TRUE)
    ca_xyz <- xyz[ca_idx, , drop = FALSE]
    ca_dist <- sqrt(rowSums((ca_xyz - matrix(centroid, nrow = nrow(ca_xyz), ncol = 3, byrow = TRUE))^2))
    max_dist <- max(ca_dist, na.rm = TRUE)

    for(i in seq_len(nres)){
      sel <- which(atoms$chain == ca_tab$chain[i] &
                   atoms$resno == ca_tab$resno[i] &
                   atoms$insert2 == ca_tab$insert2[i])
      N_i <- sel[atoms$elety[sel] == "N"][1]
      CA_i <- sel[atoms$elety[sel] == "CA"][1]
      C_i <- sel[atoms$elety[sel] == "C"][1]

      if(!is.na(N_i) && !is.na(CA_i) && !is.na(C_i) &&
         !is.na(tors$phi[i]) && !is.na(tors$psi[i])){

        dist_to_centroid <- sqrt(sum((xyz[CA_i, ] - centroid)^2))
        y <- max_dist - dist_to_centroid

        resid_i <- ca_tab$resid[i]
        grp <- NA_character_
        if(resid_i == "GLY") grp <- "G"
        if(!(resid_i %in% c("GLY", "PRO", "ILE", "VAL"))) grp <- "X"

        if(!is.na(grp)){
          all_rows[[kk]] <- data.frame(
            pdb = id,
            chain = ca_tab$chain[i],
            resno = ca_tab$resno[i],
            resid = resid_i,
            group = grp,
            phi = tors$phi[i],
            psi = tors$psi[i],
            y = y,
            stringsAsFactors = FALSE
          )
          kk <- kk + 1
        }
      }
    }
  }

  if(length(all_rows) == 0) stop("No valid residues extracted from PDBs.")

  dat <- do.call(rbind, all_rows) |>
    dplyr::filter(group %in% groups,
                  !is.na(phi), !is.na(psi),
                  y >= y_min, y <= y_max) |>
    dplyr::mutate(phi = (phi*pi/180) %% (2*pi),
                  psi = (psi*pi/180) %% (2*pi)) |>
    dplyr::select(phi, psi)

  if(nrow(dat) == 0) stop("No data left after filtering.")

  G <- grid_torus(ng)
  fhat <- ClusTorus::kde.torus(data = as.matrix(dat), eval.point = G, concentration = kappa)
  normalize_mat(matrix(fhat, nrow = ng, ncol = ng))
}

sample_pairs_from_Z <- function(Z, n){
  Z <- normalize_mat(Z)
  nphi <- nrow(Z); npsi <- ncol(Z)
  idx <- sample(length(Z), size = n, replace = TRUE, prob = as.vector(Z))
  i <- ((idx - 1) %% nphi) + 1
  j <- ((idx - 1) %/% nphi) + 1
  phi <- seq(0, 2*pi, length.out = nphi + 1)[-(nphi + 1)][i]
  psi <- seq(0, 2*pi, length.out = npsi + 1)[-(npsi + 1)][j]
  data.frame(phi = phi, psi = psi)
}
