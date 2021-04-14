run_experiment <- function(rnd.seed, Tr_specs, ho_prots, ncpus){
  require(epitopes)

  for (i in 1:nrow(Tr_specs)){
    cat("\nTraining model on data:", Tr_specs$TrData[i], "\n")
    # Fit model
    my.model <- epitopes::fit_model(data.train = readRDS(Tr_specs$TrData[i]),
                                    data.test  = ho_prots,
                                    rnd.seed   = rnd.seed,
                                    ncpus      = ncpus,
                                    num.trees  = 200,
                                    importance = "impurity",
                                    replace    = FALSE)

    # Save results
    if(!dir.exists(Tr_specs$out_dir[i])) dir.create(Tr_specs$out_dir[i],
                                                    recursive = TRUE)

    saveRDS(my.model$rf_model, paste0(Tr_specs$out_dir[i], "/model.rds"))
    saveRDS(my.model$rf_preds, paste0(Tr_specs$out_dir[i], "/preds.rds"))
    saveRDS(my.model$rf_perf,  paste0(Tr_specs$out_dir[i], "/OOS_perf.rds"))

  }
}
