calc_perf <- function(df, res_names, preds_names, pred.tr.sets = NULL){
  
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
  
  if(!is.null(pred.tr.sets)){
    method.names <- names(pred.tr.sets)
    for (j in seq_along(method.names)){
      tmp <- df %>%
        filter(!(epit_seq %in% pred.tr.sets[[method.names[j]]]$Sequence))
      
      tmp1 <- myBoot(NA, df = tmp,
                     which.cols = grep(paste0(method.names[j], "_class"), names(df)),
                     which.prob.cols = grep(paste0(method.names[j], "_prob"), names(df))) %>%
        dplyr::bind_rows() %>%
        dplyr::select(Method, sens, spec, ppv, npv, mcc, auc, accuracy) %>%
        reshape2::melt(id.vars = "Method",
                       variable.name = "Metric", value.name = "Value") %>%
        dplyr::mutate(Method = gsub("_class", "", Method))
      
      tmp2   <- lapply(1:999, FUN = myBoot,
                       df = tmp,
                       which.cols = grep(paste0(method.names[j], "_class"), names(df)),
                       which.prob.cols = grep(paste0(method.names[j], "_prob"), names(df))) %>%
        dplyr::bind_rows() %>%
        dplyr::select(Method, sens, spec, ppv, npv, mcc, auc, accuracy) %>%
        reshape2::melt(id.vars = "Method",
                       variable.name = "Metric", value.name = "Value") %>%
        dplyr::group_by(Method, Metric) %>%
        dplyr::summarise(Mean   = mean(Value, na.rm = TRUE),
                         StdErr = sd(Value, na.rm = TRUE)) %>%
        dplyr::mutate(Method = gsub("_class", "", Method))
      
      if (j == 1){
        myEst_noleak <- tmp1
        myCI_noleak  <- tmp2
      } else {
        myEst_noleak <- rbind(myEst_noleak, tmp1)
        myCI_noleak <- rbind(myCI_noleak, tmp2)
      }
    }
    
    # Consolidate performance results
    myperf_noleak <- myEst_noleak %>%
      dplyr::left_join(myCI_noleak, by = c("Method", "Metric")) %>%
      dplyr::arrange(Metric) %>%
      dplyr::mutate(Type = grepl("RF_", Method),
                    Metric = factor(toupper(Metric),
                                    levels = c("ACCURACY", "SENS", "SPEC",
                                               "PPV", "NPV", "MCC", "AUC"),
                                    ordered = TRUE),
                    Method = factor(Method,
                                    levels = c(res_names, preds_names),
                                    labels = gsub("_", "-", c(res_names, preds_names)),
                                    ordered = TRUE),
                    NoLeak = TRUE)
    
    myperf$NoLeak <- FALSE
    
    myperf <- rbind(myperf, myperf_noleak)
  }
  
  return(myperf)
}
