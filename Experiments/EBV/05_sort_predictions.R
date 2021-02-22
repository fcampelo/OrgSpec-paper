library(dplyr)

# Read analysis output (for the by-AA predictions)
myres <- readRDS("./output/analysis.rds")$myres %>%
  dplyr::select(Info_UID, Info_center_pos, Class, 
                RF_OrgSpec_prob, RF_OrgSpec_class)

for (i in seq_along(unique(myres$Info_UID))){
  
  tmp <- myres %>%
    # Extract only the i-th protein
    dplyr::filter(Info_UID == unique(myres$Info_UID)[i]) %>%
    # Flag the start of each individual predicted epitope
    dplyr::mutate(IsBreak = (RF_OrgSpec_class == 1 & dplyr::lag(RF_OrgSpec_class, default = 1) == -1)) %>%
    dplyr::filter(RF_OrgSpec_class == 1) 
  
  if (nrow(tmp) > 0){
    tmp <- tmp %>%
      # Give each predicted epitope a number
      dplyr::mutate(EpNumber = cumsum(IsBreak)) %>%
      # Summarise the average probability of that epitope as the strength of its prediction
      dplyr::group_by(EpNumber) %>%
      dplyr::summarise(Info_UID  = first(Info_UID),
                       start_pos = min(Info_center_pos),
                       end_pos   = max(Info_center_pos),
                       length    = end_pos - start_pos + 1,
                       Prob      = mean(RF_OrgSpec_prob))
  }
  
  # Consolidate results
  if (i == 1){
    mypreds <- tmp
  } else {
    mypreds <- rbind(mypreds, tmp)
  }
}

mypreds <- mypreds %>%
  arrange(dplyr::desc(Prob)) %>%
  dplyr::select(-EpNumber)

saveRDS(mypreds, "./output/consolidated_preds.rds")
