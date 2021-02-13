read_SVMtrip <- function(res.path, prots.path){
  require(XML)
  require(dplyr)
  
  # Read proteins
  proteins <- readRDS(prots.path) %>%
    as.data.frame() %>%
    dplyr::select(!starts_with("feat_"))
  
  fn  <- dir(res.path, full.names = TRUE, pattern = ".html")
  ids <- gsub(".html", "", dir(res.path, full.names = FALSE, pattern = ".html"))
  
  for (i in seq_along(fn)){
    
    # Isolate protein and prepare variable for results
    prot  <- proteins %>% dplyr::filter(Info_UID == ids[i])
    prot$SVMtrip_prob  <- 0
    prot$SVMtrip_class <- -1
    
    # Read and preprocess results file
    tmp   <- XML::readHTMLTable(fn[i], header = TRUE)[[1]]
    
    if(!is.null(tmp)){
      pos <- strsplit(tmp$Location, split = " - ")
      
      for (j in seq_along(pos)){
        prot$SVMtrip_prob[pos[[j]][1]:pos[[j]][2]] <- as.numeric(tmp$Score[j])
        prot$SVMtrip_class[pos[[j]][1]:pos[[j]][2]] <- 1
      }
    }
    
    if (i == 1){
      myres <- prot
    } else {
      myres <- rbind(myres, prot)
    }
  }
  
  return(myres[, c("Info_UID", "Info_center_pos", "SVMtrip_prob", "SVMtrip_class")])
}
