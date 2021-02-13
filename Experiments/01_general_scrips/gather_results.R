gather_results <- function(mydata_path, res_paths, preds_paths){
  require(dplyr)
  
  # Preprocess results
  mydata  <- readRDS(mydata_path) %>%
    as.data.frame() %>%
    dplyr::mutate(Class = as.factor(Class)) %>%
    dplyr::mutate(across(starts_with("Info"), as.character))
  
  # Read our results
  for (i in seq_along(res_paths$name)){
    if(dir.exists(res_paths$path[i])){
      res <- readRDS(paste0(res_paths$path[i], "/preds.rds")) %>%
        dplyr::mutate(`.pred_class` = epitopes::smooth_predictions(as.numeric(as.character(`.pred_class`)),
                                                                   window_size = nchar(mydata$Info_window_seq[1]),
                                                                   type = "mode"),
                      Info_center_pos = as.numeric(Info_center_pos)) %>%
        dplyr::rename(!!paste0(res_paths$name[i], "_class") := .pred_class,
                      !!paste0(res_paths$name[i], "_prob")  := .pred_prob)
      
      if (!exists("myres")){
        myres <- res
      } else {
        myres <- myres %>%
          dplyr::left_join(res, by = c("Info_UID", "Info_center_pos", "Class"))
      }
    }
  }
  
  # Read results from literature predictors
  for (i in seq_along(preds_paths$name)){
    if(dir.exists(preds_paths$path[i]) || file.exists(preds_paths$path[i])){
      fx <- source(preds_paths$read.fun[i])$value
      res <- fx(res.path = preds_paths$path[i], prots.path = mydata_path)
      myres <- myres %>%
        dplyr::left_join(res, by = c("Info_UID", "Info_center_pos"))
    }
  }
  
  return(myres)
}
