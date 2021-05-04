run_analysis <- function(rnd.seed, epitopes, data.train, pred.train.sets,
                         ho_prots, ho_peps, myres_files, preds_paths, ncpus,
                         n.boot = 999){
  
  require(dplyr)
  set.seed(rnd.seed)
  
  tr.epits <- unique(data.train$Info_epitope_id)
  pred.train.sets$RF_OrgSpec <- as.data.frame(epitopes) %>%
    dplyr::filter(epitope_id %in% tr.epits) %>%
    dplyr::select(Sequence = epit_seq)
  pred.train.sets$RF_Hybrid  <- pred.train.sets$RF_OrgSpec
  pred.train.sets$RF_Heter   <- data.frame(Sequence = "1234")
  
  # Consolidate results
  cat("\nConsolidating results...\n")
  my_results  <- gather_results(ho_prots    = ho_prots,
                                ho_peps     = ho_peps,
                                res_paths   = myres_files,
                                preds_paths = preds_paths)
  
  df <- my_results$results_by_peptide %>%
    left_join(dplyr::select(as.data.frame(epitopes), epitope_id, epit_seq),
              by = c("Info_epitope_id" = "epitope_id")) %>%
    dplyr::distinct()
  
  cat("\nCalculating performance...\n")
  myperf_pep <- calc_perf(df           = df,
                          res_names    = myres_files$name,
                          preds_names  = preds_paths$name,
                          pred.tr.sets = pred.train.sets,
                          ncpus        = ncpus,
                          nboot        = n.boot)
  
  # Prepare plot
  mp_pep <- make_plot3(df = myperf_pep, linepos = 3.5, plot_rows = 2)
  
  if(!dir.exists("./figures")) dir.create("./figures")
  ggplot2::ggsave("./figures/res_byPep.png", plot = mp_pep,
                  width = 8, height = 5, units = "in")
  
  
  
  # === hypothesis testing === #
  # 1) Considering full holdout set
  # Reference method: RF (organism specific)
  # Using results by peptide
  
  cat("\nCalculating p-values (part 1)... \n")
  
  # Estimate bootstrap distribution of differences UNDER THE NULL HYPOTHESES
  ref <- "RF_OrgSpec"
  diff_boot <- parallel::mclapply(1:n.boot,
                                  FUN = boot_diffs,
                                  df  = df,
                                  ref  = ref,
                                  mc.cores = ncpus,
                                  mc.set.seed = rnd.seed) %>%
    dplyr::bind_rows() %>%
    dplyr::rename_with(toupper)
  
  
  # Get actual observed differences
  diff_obs <- boot_diffs(NA, df = df, ref = ref) %>%
    dplyr::rename_with(toupper)
  
  diff_boot <- rbind(diff_obs, diff_boot)
  
  Pvals.raw <- diff_boot %>%
    dplyr::group_by(METHOD1, METHOD2) %>%
    dplyr::summarise(across(!starts_with("METHOD"),
                            function(x){
                              (sum(abs(x) >= abs(dplyr::first(x)))) / n()}),
                     .groups = "drop")
  
  # Correct p-values for MHT using Benjamini-Hochberg
  Pvals <- Pvals.raw %>%
    dplyr::mutate(across(!starts_with("METHOD"),
                         ~p.adjust(.x, method = "holm")),
                  NoLeak = FALSE)
  
  
  
  
  # 2) Considering only zero-leakage subsets
  cat("\nCalculating p-values (part 2)... \n")
  noleak_res <- lapply(pred.train.sets,
                       function(X){
                         filter(df, !(epit_seq %in% X$Sequence))})
  
  # Estimate bootstrap distribution of differences UNDER THE NULL HYPOTHESES
  diff_boot <- parallel::mclapply(1:n.boot,
                                  FUN = boot_diffs_nopair,
                                  df  = noleak_res,
                                  ref  = ref,
                                  mc.cores = ncpus,
                                  mc.set.seed = rnd.seed) %>%
    dplyr::bind_rows() %>%
    dplyr::rename_with(toupper)
  
  # Get actual observed differences
  diff_obs <- boot_diffs_nopair(NA, df = noleak_res, ref = ref) %>%
    dplyr::rename_with(toupper)
  
  diff_boot <- rbind(diff_obs, diff_boot)
  
  Pvals_noleak.raw <- diff_boot %>%
    dplyr::group_by(METHOD1, METHOD2) %>%
    dplyr::summarise(across(!starts_with("METHOD"),
                            function(x){
                              (sum(abs(x) >= abs(dplyr::first(x)))) / (n() - 1)}),
                     .groups = "drop")
  
  # Correct p-values for MHT using Benjamini-Hochberg
  Pvals_noleak <- Pvals_noleak.raw %>%
    dplyr::mutate(across(!starts_with("METHOD"),
                         ~p.adjust(.x, method = "holm")),
                  NoLeak = TRUE)
  
  Pvals <- rbind(Pvals, Pvals_noleak)
  Pvals.raw <- cbind(rbind(Pvals.raw, Pvals_noleak.raw), 
                     NoLeak = Pvals$NoLeak)
  
  predlist <- sort_predictions(my_results$results_by_position)
  
  saveRDS(object = c(list(myres       = my_results$results_by_position,
                          myres_pep   = my_results$results_by_peptide,
                          myperf_pep  = myperf_pep,
                          myplot_pep  = mp_pep,
                          Pvals_pep   = Pvals,
                          Pvals_pep.raw = Pvals.raw),
                     predlist),
          file = "./output/analysis.rds")
  
  cat("\nDone!\n")
}
