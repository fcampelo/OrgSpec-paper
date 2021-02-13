read_ABCPred <- function(res.path, prots.path){
  require(XML)
  require(dplyr)

  # Read proteins
  proteins <- readRDS(prots.path) %>%
    as.data.frame() %>%
    dplyr::select(!starts_with("feat_"))

  fn <- dir(res.path, full.names = TRUE, pattern = ".html")

  for (i in seq_along(fn)){
    # Read and preprocess results file
    tmp   <- XML::readHTMLTable(fn[i], header = FALSE)
    preds <- tmp[[2]]
    names(preds) <- preds[1, ]
    preds <- preds[-c(1, which(is.na(preds[, 1]))), -5]
    preds[, c(1, 3, 4)] <- lapply(preds[, c(1, 3, 4)], as.numeric)

    # Isolate protein and prepare variable for results
    prot  <- proteins %>% dplyr::filter(Info_UID == tmp[[1]][1, 2])
    prot$ABCpred_class <- -1
    prot$ABCpred_prob  <- 0

    for (j in 1:nrow(preds)){
      stpos <- preds$`Start position`[j]
      idx   <- stpos:(stpos + nchar(preds$Sequence[j]) - 1)
      prot$ABCpred_class[idx] <- 1
      prot$ABCpred_prob[idx]  <- pmax(prot$ABCpred_prob[idx], preds$Score[j])
    }

    if (i == 1){
      myres <- prot
    } else {
      myres <- rbind(myres, prot)
    }
  }

  return(myres[, c("Info_UID", "Info_center_pos", "ABCpred_prob", "ABCpred_class")])
}
