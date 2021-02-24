sort_predictions <- function(myres){
  
  myres <- myres %>%
    dplyr::select(Info_UID, Info_center_pos, Class, 
                  RF_OrgSpec_prob, RF_OrgSpec_class)
  
  for (i in seq_along(unique(myres$Info_UID))){
    
    tmp <- myres %>%
      # Extract only the i-th protein
      dplyr::filter(Info_UID == unique(myres$Info_UID)[i]) %>%
      # Flag the start of each individual predicted epitope
      dplyr::mutate(IsBreak = (RF_OrgSpec_class == 1 & dplyr::lag(RF_OrgSpec_class, default = 1) == -1)) 
    
    tmp2 <- tmp %>%
      dplyr::filter(RF_OrgSpec_class == 1) 
    
    if (nrow(tmp2) > 0){
      tmp2 <- tmp2 %>%
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
      myprobs <- tmp
      mypreds <- tmp2
    } else {
      myprobs <- rbind(myprobs, tmp)
      mypreds <- rbind(mypreds, tmp2)
    }
  }
  
  mypreds <- mypreds %>%
    arrange(dplyr::desc(Prob)) %>%
    dplyr::select(-EpNumber)
  
  myprobs <- myprobs %>%
    dplyr::select(-IsBreak)
  
  return(list(mypreds = mypreds,myprobs = myprobs))
}
