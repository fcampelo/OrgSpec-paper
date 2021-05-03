generate_datasets <- function(target.Org, hostIDs, removeIDs, rnd.seed,
                              heter.nPos, heter.nNeg,
                              epitopes, proteins, taxonomy,
                              coverage_threshold, identity_threshold,
                              ncpus){
  require(epitopes)
  require(dplyr)

  # Extract Organism-specific datasets (train / test)
  cat("\n Assembling OrgSpec Dataset\n")
  data.OrgSpec <- epitopes::make_OrgSpec_datasets(epitopes, proteins, taxonomy,
                                                  orgIDs             = target.Org,
                                                  hostIDs            = hostIDs,
                                                  save_folder        = "./data/splits",
                                                  coverage_threshold = coverage_threshold,
                                                  identity_threshold = identity_threshold,
                                                  ncpus              = ncpus)


  # Extract Heterogeneous and Hybrid training datasets
  cat("\n Assembling Heterogeneous Dataset\n")
  data.Heter <- epitopes::make_heterogeneous_dataset(epitopes, proteins, taxonomy,
                                                     nPos            = heter.nPos,
                                                     nNeg            = heter.nNeg,
                                                     removeIDs       = c(removeIDs, target.Org),
                                                     hostIDs         = hostIDs,
                                                     save_folder     = "./data/heterogeneous_data",
                                                     rnd.seed        = rnd.seed,
                                                     ncpus           = ncpus)

  # Assemble Hybrid data set
  data.Hybrid = rbind(data.OrgSpec[[1]]$wdf, data.Heter)
  saveRDS(data.Hybrid, file = "./data/heterogeneous_data/df_hybrid.rds")


  # Assemble dataset for full hold-out proteins
  cat("\n Assembling Hold-out Proteins Dataset\n")
  data.holdout_prots <- make_proteins_dataset(epitopes, proteins, taxonomy,
                                              prot_IDs = unique(data.OrgSpec[[2]]$wdf$Info_protein_id),
                                              orgIDs      = target.Org,
                                              hostIDs     = hostIDs,
                                              removeIDs   = removeIDs,
                                              save_folder = "./data/splits",
                                              ncpus       = ncpus)

  cat("\n Done!\n")
  invisible(TRUE)
}
