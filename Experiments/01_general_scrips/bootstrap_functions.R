# Bootstrap functions

# Use Bootstrap to calculate confidence intervals on performance:
myBoot <- function(i, df, which.cols){
  idx <- 1:nrow(df)
  if (!is.na(i)) idx <- sample.int(nrow(df), replace = TRUE)
  tmp <- as.data.frame(df[idx, ])
  res <- t(sapply(which.cols,
                  function(j){epitopes::calc_performance(truth = tmp$Class,
                                                         pred  = tmp[, j])}))
  res <- lapply(as.data.frame(res), as.numeric)
  res$Method <- names(tmp)[which.cols]
  return(res)
}



# Create a single estimate of performance difference:
# - (i)  under the resampling-derived null hypothesis (if !is.na(i))
# - (ii) using the actual (observed) results (if is.na(i))
boot_diffs <- function(i, df, ref.col, cmp.cols){
  if (!is.na(i)) {
    # Resample rows
    df <- df[sample.int(nrow(df), replace = TRUE), ]
    
    # Under the null hypothesis it makes no difference to 
    # mix the reference and comparison data at half the positions
    # (randomly selected)
    idx <- sample.int(nrow(df), size = ceiling(nrow(df) / 2), 
                      replace = FALSE)
  } 
  
  # Calculate perf differences between reference and all comparison
  # methods
  for (j in seq_along(cmp.cols)){
    yref <- df[, ref.col]
    ycmp <- df[, cmp.cols[j]]
    if (!is.na(i)){
      yref[idx] <- ycmp[idx]
      ycmp[idx] <- df[idx, ref.col]
    }
    
    pref <- epitopes::calc_performance(truth = df$Class,
                                       pred  = yref)
    pcmp <- epitopes::calc_performance(truth = df$Class,
                                       pred  = ycmp)
    
    pdiff <- (pref - pcmp) %>%
      dplyr::mutate(Method1 = gsub("_class", "", names(df)[ref.col],
                                   fixed = TRUE),
                    Method2 = gsub("_class", "", names(df)[cmp.cols[j]],
                                   fixed = TRUE)) %>%
      dplyr::select(Method1, Method2, sens, spec, 
                    accuracy, ppv, npv, mcc)
    
    if (j == 1){
      diffs <- pdiff
    } else {
      diffs <- rbind(diffs, pdiff)
    }
  }
  
  return(diffs)
}
