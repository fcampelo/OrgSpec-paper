gather_results_byPep <- function(mypeps_path, myres){
  
  myres2 <- readRDS(mypeps_path) %>%
    dplyr::select(Info_epitope_id, Info_protein_id, Info_center_pos, Class) %>%
    dplyr::group_by(Info_epitope_id) %>%
    dplyr::summarise(Info_protein_id = first(Info_protein_id),
                     Info_start      = min(Info_center_pos),
                     Info_end        = max(Info_center_pos),
                     Class           = sign(sum(as.numeric(as.character(Class))) - 1e-3))
  
  for (i in seq_along(myres2$Info_epitope_id)){
    tmp <- myres %>% 
      dplyr::filter(Info_UID == myres2$Info_protein_id[i],
                    Info_center_pos >= myres2$Info_start[i],
                    Info_center_pos <= myres2$Info_end[i]) %>%
      dplyr::select(-Info_UID, -Info_center_pos, -Class) %>%
      dplyr::summarise(across(.cols = ends_with("class"), 
                              .fns = ~ sign(sum(.x) - 1e-3)),
                       across(.cols = ends_with("prob"), 
                              .fns = mean))
    if (i == 1){
      pepres <- tmp
    } else {
      pepres <- rbind(pepres, tmp)
    }
  }
  
  return(cbind(myres2, pepres))
}
