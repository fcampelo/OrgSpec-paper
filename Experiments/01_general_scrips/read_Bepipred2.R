read_Bepipred2 <- function(res.path, ...){
  fn <- dir(res.path, full.names = TRUE, pattern = ".csv")

  for (i in seq_along(fn)){
    tmp <- read.csv(fn[i], header = TRUE, stringsAsFactors = FALSE)
    tmp <- tmp[, c("Entry", "Position", "EpitopeProbability")]

    if (i == 1){
      myres <- tmp
    } else {
      myres <- rbind(myres, tmp)
    }
  }
    names(myres) <- c("Info_UID", "Info_center_pos", "Bepipred2_prob")
    myres$Bepipred2_class <- -1 + 2 * (myres$Bepipred2_prob > 0.5)


  return(myres)
}
