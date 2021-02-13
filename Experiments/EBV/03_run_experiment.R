library(tidymodels)
source("../01_general_scrips/fit_models.R")

rn.seed <- 20210107
fTe      <- "./data/splits/holdout_prots_w.rds"
Tr_specs <- data.frame(
  TrData  = c("./data/splits/01_training.rds",
             "./data/heterogeneous_data/df_heterogeneous.rds",
             "./data/heterogeneous_data/df_hybrid.rds"),
  out_dir = c("./output/TrOrgSpec",
             "./output/TrHeter",
             "./output/TrHybrid"))

t0 <- Sys.time()
res <- vector("list", nrow(Tr_specs))
for (i in 1:nrow(Tr_specs)){
  ti <- Sys.time()
  cat("\nTraining models on data set", i, "of", nrow(Tr_specs), "\n")
  res[[i]] <- fit_models(train_path = Tr_specs$TrData[i],
                         test_path = fTe, nfolds = 5, rn.seed = rn.seed)

  if(!dir.exists(Tr_specs$out_dir[i])) dir.create(Tr_specs$out_dir[i],
                                                  recursive = TRUE)
  saveRDS(res[[i]]$rf_model, paste0(Tr_specs$out_dir[i], "/model.rds"))
  saveRDS(res[[i]]$rf_preds, paste0(Tr_specs$out_dir[i], "/preds.rds"))
  saveRDS(res[[i]]$rf_perf, paste0(Tr_specs$out_dir[i], "/OOS_perf.rds"))

  cat("\nPerformance on test set:\n")
  print(res[[i]]$rf_perf)

  Ti <- difftime(Sys.time(), ti)
  Tt <- difftime(Sys.time(), t0)
  cat("\n===============================================",
      "\nData set", i, "of", nrow(Tr_specs), ": finished!",
      "\nExperiment started: ", sprintf("%s", t0),
      "\nModel time:", signif(Ti, 2), units(Ti),
      "\nTotal Elapsed time:", signif(Tt, 2), units(Tt),
      "\n===============================================\n")
}


