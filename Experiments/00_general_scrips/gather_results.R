gather_results <- function(ho_prots, ho_peps, res_paths, preds_paths){
  require(dplyr)
  require(epitopes)

  # Preprocess results
  mydata  <- as.data.frame(ho_prots) %>%
    dplyr::mutate(Class = as.factor(Class)) %>%
    dplyr::mutate(across(starts_with("Info"), as.character))

  # Read results: proposed model + controls
  for (i in seq_along(res_paths$name)){
    if(file.exists(res_paths$file[i])){
      res <- readRDS(myres_files$file[i]) %>%
        dplyr::mutate(pred_class = smooth_predictions(as.numeric(as.character(pred_class)),
                                                      window_size = nchar(ho_prots$Info_window_seq[1]),
                                                      type    = "minsize",
                                                      minsize = 8),
                      Info_center_pos = as.numeric(Info_center_pos)) %>%
        dplyr::rename(!!paste0(myres_files$name[i], "_class") := pred_class,
                      !!paste0(myres_files$name[i], "_prob")  := pred_prob)

      if (i == 1){
        myres <- res
      } else {
        myres <- myres %>%
          dplyr::left_join(res, by = c("Info_UID", "Info_center_pos", "Class"))
      }
    }
  }

  # Read results from literature predictors
  for (i in seq_along(preds_paths$name)){
    # If path is a folder:
    if(dir.exists(preds_paths$path[i])){
      fn <- dir(preds_paths$path[i], full.names = TRUE,
                pattern = "html|csv|txt")

      for (j in seq_along(fn)){
        protID <- gsub(preds_paths$path[i], "", fn[j]) %>%
          gsub(pattern = "\\.[a-zA-Z]+|/", replacement = "")

        tmp <- do.call(preds_paths$read.fun[i],
                       args = list(res.file = fn[j],
                                   protID = protID,
                                   proteins = myres))
        if (j == 1){
          res <- tmp
        } else {
          res <- rbind(res, tmp)
        }
      }

      # If path is a file:
    } else if (file.exists(preds_paths$path[i])) {
      res <- do.call(preds_paths$read.fun[i],
                     args = list(res.file = preds_paths$path[i],
                                 proteins = myres))
    } else {
      res = data.frame(Info_UID = NA, Info_center_pos = NA)
    }

    myres <- myres %>%
      dplyr::left_join(res, by = c("Info_UID", "Info_center_pos"))

  }

  # ========================================================================== #
  # Consolidate results by peptide (rather than by amino acid position)
  myres2 <- as.data.frame(ho_peps) %>%
    dplyr::select(Info_epitope_id, Info_protein_id, Info_center_pos, Class) %>%
    dplyr::group_by(Info_epitope_id) %>%
    dplyr::summarise(Info_protein_id = first(Info_protein_id),
                     Info_start      = min(Info_center_pos),
                     Info_end        = max(Info_center_pos),
                     Class           = sign(sum(as.numeric(as.character(Class))) - 1e-3))


  for (i in seq_along(myres2$Info_epitope_id)){
    tmp <- myres %>%
      dplyr::filter(Info_UID        == myres2$Info_protein_id[i],
                    Info_center_pos >= myres2$Info_start[i],
                    Info_center_pos <= myres2$Info_end[i]) %>%
      dplyr::select(-Info_center_pos, -Class) %>%
      dplyr::summarise(Info_UID = first(Info_UID),
                       across(.cols = ends_with("class"),
                              .fns = ~ sign(sum(.x) - 1e-3)),
                       across(.cols = ends_with("prob"),
                              .fns = mean))
    if (i == 1){
      pepres <- tmp
    } else {
      pepres <- rbind(pepres, tmp)
    }
  }

  myres2 <- cbind(myres2, pepres) %>%
    dplyr::select(Info_epitope_id, Info_protein_id, everything(), -Info_UID)

  return(list(results_by_position = myres,
              results_by_peptide  = myres2))
}
