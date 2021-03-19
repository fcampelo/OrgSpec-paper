# Bootstrap functions

# Use Bootstrap to calculate confidence intervals on performance:
myBoot <- function(i, df, which.cols, which.prob.cols = NULL){
  idx <- 1:nrow(df)
  if (!is.na(i)) idx <- sample.int(nrow(df), replace = TRUE)
  tmp <- as.data.frame(df[idx, ])
  if (!is.null(which.prob.cols)){
    res <- t(mapply(function(j, k){epitopes::calc_performance(truth = tmp$Class,
                                                              pred  = tmp[, j],
                                                              prob  = tmp[, k])},
                    which.cols, which.prob.cols))
  } else {
    res <- t(sapply(which.cols,
                    function(j){epitopes::calc_performance(truth = tmp$Class,
                                                           pred  = tmp[, j])}))
  }
  res <- lapply(as.data.frame(res), as.numeric)
  res$Method <- names(tmp)[which.cols]
  return(res)
}



# Create a single estimate of performance difference:
# - (i)  under the resampling-derived null hypothesis (if !is.na(i))
# - (ii) using the actual (observed) results (if is.na(i))
boot_diffs <- function(i, df, ref){
  ref.pred  <- grep(paste0(ref, "_class"), names(df))
  ref.prob  <- grep(paste0(ref, "_prob"), names(df))
  cmp.pred  <- grep("_class", names(df))
  cmp.pred  <- cmp.pred[-which(cmp.pred == ref.pred)]
  cmp.prob  <- grep("_prob", names(df))
  cmp.prob  <- cmp.prob[-which(cmp.prob == ref.prob)]
  
  
  if (!is.na(i)) {
    # Resample rows
    df <- df[sample.int(nrow(df), replace = TRUE), ]
    
    # Under the null hypothesis it makes no difference to
    # mix the reference and comparison data at approx. half the positions
    # (randomly selected)
    idx <- which(runif(nrow(df)) <= 0.5)
  }
  
  # Calculate performance differences between reference and all comparison
  # methods
  for (j in seq_along(cmp.pred)){
    yref <- df[, ref.pred]
    ycmp <- df[, cmp.pred[j]]
    pref <- df[, ref.prob]
    pcmp <- df[, cmp.prob[j]]
    if (!is.na(i)){
      yref[idx] <- ycmp[idx]
      ycmp[idx] <- df[idx, ref.pred]
      pref[idx] <- pcmp[idx]
      pcmp[idx] <- df[idx, ref.prob]
    }
    
    pref <- epitopes::calc_performance(truth = df$Class,
                                       pred  = yref,
                                       prob  = pref)
    pcmp <- epitopes::calc_performance(truth = df$Class,
                                       pred  = ycmp,
                                       prob  = pcmp)
    
    pdiff <- (pref - pcmp) %>%
      dplyr::mutate(Method1 = gsub("_class", "", names(df)[ref.pred],
                                   fixed = TRUE),
                    Method2 = gsub("_class", "", names(df)[cmp.pred[j]],
                                   fixed = TRUE)) %>%
      dplyr::select(Method1, Method2, sens, spec,
                    accuracy, ppv, npv, mcc, auc)
    
    if (j == 1){
      diffs <- pdiff
    } else {
      diffs <- rbind(diffs, pdiff)
    }
  }
  
  return(diffs)
}


# Create a single estimate of performance difference:
# - (i)  under the resampling-derived null hypothesis (if !is.na(i))
# - (ii) using the actual (observed) results (if is.na(i))
boot_diffs_nopair <- function(i, df, ref){
  require(dplyr)
  
  ref.pred  <- df[[ref]] %>% select(Info_epitope_id, Class, starts_with(ref)) %>%
    dplyr::rename("Pred" = !!paste0(ref, "_class"),
                  "Prob" = !!paste0(ref, "_prob"))
  df[[ref]] <- NULL
  preds <- names(df)
  df <- lapply(seq_along(df), 
               function(k){
                 df[[k]] %>%
                   select(Info_epitope_id, Class, starts_with(names(df)[k])) %>%
                   mutate(RepBy = (runif(nrow(df[[k]])) <= 0.5) * sample.int(size = nrow(df[[k]]), 
                                                                             n = nrow(ref.pred), 
                                                                             replace = TRUE),
                          RepBy = replace(RepBy, which(RepBy == 0), NA)) %>%
                   dplyr::rename("Pred" = !!paste0(names(df)[k], "_class"),
                                 "Prob" = !!paste0(names(df)[k], "_prob"))})
  names(df) <- preds
  
  # Calculate performance differences between reference and all comparison
  # methods
  for (k in seq_along(df)){
    ref.pred <- ref.pred %>%
      mutate(RepBy = (runif(nrow(ref.pred)) <= 0.5) * sample.int(size = nrow(ref.pred), 
                                                                n = nrow(df[[k]]), 
                                                                replace = TRUE),
             RepBy = replace(RepBy, which(RepBy == 0), NA))
    yref <- ref.pred[, c("Pred", "Prob", "Class")]
    ycmp <- df[[k]][, c("Pred", "Prob", "Class")]
    if (!is.na(i)){
      yref[which(!is.na(ref.pred$RepBy)), ] <- ycmp[ref.pred$RepBy[!is.na(ref.pred$RepBy)], c("Pred", "Prob", "Class")]
      ycmp[which(!is.na(df[[k]]$RepBy)), ]  <- ref.pred[df[[k]]$RepBy[!is.na(df[[k]]$RepBy)], c("Pred", "Prob", "Class")]
    }
    
    pref <- epitopes::calc_performance(truth = yref$Class,
                                       pred  = yref$Pred,
                                       prob  = yref$Prob)
    pcmp <- epitopes::calc_performance(truth = ycmp$Class,
                                       pred  = ycmp$Pred,
                                       prob  = ycmp$Prob)
    
    pdiff <- (pref - pcmp) %>%
      dplyr::mutate(Method1 = ref,
                    Method2 = names(df)[k]) %>%
      dplyr::select(Method1, Method2, sens, spec,
                    accuracy, ppv, npv, mcc, auc)
    
    if (k == 1){
      diffs <- pdiff
    } else {
      diffs <- rbind(diffs, pdiff)
    }
  }
  
  return(diffs)
}

