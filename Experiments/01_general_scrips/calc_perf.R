calc_perf <- function(df, res_names, preds_names){

  myEst <- myBoot(NA, df = df,
                  which.cols = grep("_class", names(df)),
                  which.prob.cols = grep("_prob", names(df))) %>%
    dplyr::bind_rows() %>%
    dplyr::select(Method, sens, spec, ppv, npv, mcc, auc, accuracy) %>%
    reshape2::melt(id.vars = "Method",
                   variable.name = "Metric", value.name = "Value") %>%
    dplyr::mutate(Method = gsub("_class", "", Method))

  myCI   <- lapply(1:999, FUN = myBoot,
                   df = df,
                   which.cols = grep("_class", names(df)),
                   which.prob.cols = grep("_prob", names(df))) %>%
    dplyr::bind_rows() %>%
    dplyr::select(Method, sens, spec, ppv, npv, mcc, auc, accuracy) %>%
    reshape2::melt(id.vars = "Method",
                   variable.name = "Metric", value.name = "Value") %>%
    dplyr::group_by(Method, Metric) %>%
    dplyr::summarise(Mean   = mean(Value, na.rm = TRUE),
                     StdErr = sd(Value, na.rm = TRUE)) %>%
                     #CIl    = quantile(Value, sig.level / 2, na.rm = TRUE),
                     #CIu    = quantile(Value, 1 - sig.level / 2, na.rm = TRUE)) %>%
    dplyr::mutate(Method = gsub("_class", "", Method))

  # Consolidate performance results
  myperf <- myEst %>%
    dplyr::left_join(myCI, by = c("Method", "Metric")) %>%
    dplyr::arrange(Metric) %>%
    dplyr::mutate(Type = grepl("RF_", Method),
                  Metric = factor(toupper(Metric),
                                  levels = c("ACCURACY", "SENS", "SPEC",
                                             "PPV", "NPV", "MCC", "AUC"),
                                  ordered = TRUE),
                  Method = factor(Method,
                                  levels = c(res_names, preds_names),
                                  labels = gsub("_", "-", c(res_names, preds_names)),
                                  ordered = TRUE))
  return(myperf)
}
