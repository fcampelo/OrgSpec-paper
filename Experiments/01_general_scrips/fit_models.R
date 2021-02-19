fit_models <- function(train_path, test_path = NULL, rn.seed = NULL){
  require(dplyr)
  require(tidymodels)
  
  # Load sets
  if (grepl(".csv", train_path)){
    mydata.tr <- read.csv(train_path, header = TRUE)
  } else {
    mydata.tr <- readRDS(train_path) %>%
      as.data.frame()
  }
  
  if (grepl(".csv", test_path)){
    mydata.te <- read.csv(test_path, header = TRUE)
  } else {
    mydata.te <- readRDS(test_path) %>%
      as.data.frame()
  }
  
  mydata.tr <- mydata.tr %>%
    dplyr::mutate(Class = as.factor(Class)) %>%
    dplyr::mutate(Info_UID = Info_protein_id) %>%
    dplyr::mutate(across(starts_with("Info"), as.character))
  
  mydata.te <- mydata.te %>%
    as.data.frame() %>%
    dplyr::mutate(Class = as.factor(Class)) %>%
    dplyr::mutate(across(starts_with("Info"), as.character))
  
  # add small decimal bias to 0-variance components to prevent accidental conversion to
  # integer.
  idx <- which(sapply(mydata.tr, function(x){length(unique(x))}) == 1 & !grepl("Info", names(mydata.tr)))
  if(length(idx)) mydata.tr[, idx] <- mydata.tr[, idx] + 1e-6
  idx <- which(sapply(mydata.te, function(x){length(unique(x))}) == 1 & !grepl("Info", names(mydata.te)))
  if(length(idx)) mydata.te[, idx] <- mydata.te[, idx] + 1e-6
  
  # Remove columns not present in both dataframes
  idx <- which(!names(mydata.tr) %in% names(mydata.te))
  if(length(idx) > 0) mydata.tr[, idx] <- NULL
  idx <- which(!names(mydata.te) %in% names(mydata.tr))
  if(length(idx) > 0) mydata.te[, idx] <- NULL
  
  if(!is.null(rn.seed)) set.seed(rn.seed)
  
  # Build recipe. No feature selection is performed.
  rec <- recipes::recipe(Class ~ ., data = mydata.tr) %>%
    update_role(starts_with("Info"), new_role = "ID")
  
  # Build model
  rf_mod <- rand_forest(trees = 200) %>%
    set_engine("ranger",
               importance  = "impurity",
               num.threads = parallel::detectCores() - 1,
               replace     = FALSE) %>%
    set_mode("classification")
  
  # Set up workflow
  wflow <- workflow() %>%
    add_recipe(rec) %>%
    add_model(rf_mod)
  
  rf_fit <- wflow %>%
    fit(mydata.tr)
  
  rf_class <- predict(rf_fit, mydata.te, type = "class")
  rf_probs <- predict(rf_fit, mydata.te, type = "prob")
  rf_preds <- mydata.te %>%
    dplyr::select(Info_UID, Info_center_pos, Class) %>%
    dplyr::bind_cols(.pred_prob = rf_probs$.pred_1, rf_class) %>%
    dplyr::mutate(Class       = as.numeric(as.character(Class)),
                  .pred_class = as.numeric(as.character(.pred_class)))
  
  multi_metric <- metric_set(accuracy, sens, spec, ppv, npv, mcc)
  
  outlist <- list(rf_CV_perf  = NULL,
                  rf_model    = rf_fit,
                  rf_preds    = rf_preds,
                  rf_perf     = rf_preds %>%
                    multi_metric(truth = as.factor(Class),
                                 estimate = as.factor(.pred_class)))
  
  return(outlist)
}
