# Assemble a data set of peptides from heterogeneous organisms
# (sampled from IEDB entries)
make_heterogeneous_set <- function(nPos, nNeg,
                                   proteins, epitopes,
                                   tax_load_file,
                                   Orgs.to.Remove, hostIDs,
                                   rnd.seed,
                                   ncores = 1){
  require(data.table)

  # Join epitope/protein data
  jdf <- epitopes::prepare_join_df(epitopes, proteins)

  jdf <- epitopes::filter_epitopes(jdf,
                                   removeIDs     = Orgs.to.Remove,
                                   hostIDs       = hostIDs,
                                   tax_load_file = tax_load_file)

  # Randomise organism IDs
  set.seed(rnd.seed)
  orgIds <- sample(unique(jdf$sourceOrg_id))

  # Compose heterogeneous dataset
  for (i in seq_along(orgIds)){
    newSet <- jdf[sourceOrg_id %in% orgIds[1:i], ]
    cl     <- as.numeric(table(factor(newSet$Class,
                                      levels = c(-1, 1))))
    if (cl[1] >= nNeg & cl[2] >= nPos) break
  }

  if (cl[1] > nNeg){
    idx  <- sample(which(newSet$Class == -1), size = cl[1] - nNeg, replace = FALSE)
    newSet <- newSet[-idx, ]
  }

  if (cl[2] > nPos){
    idx  <- sample(which(newSet$Class == 1), size = cl[2] - nPos, replace = FALSE)
    newSet <- newSet[-idx, ]
  }

  # Generate dataset with features
  # Prepare windowed representation
  wdf <- epitopes::make_window_df(newSet, ncpus = ncores)
  wdf <- epitopes::calc_features(wdf, max.N = 2, ncpus = ncores)

  return(wdf)
}
