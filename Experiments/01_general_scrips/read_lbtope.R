read_lbtope <- function(res.path, ...){

  myres  <- read.csv(res.path, header = FALSE, sep = "\t")
  pos    <- c(which(is.na(myres$V2)), nrow(myres) + 1)
  ids <- sapply(myres$V1[pos],
                function(x){
                  gsub(">", "", strsplit(x, split = " ")[[1]][4])},
                USE.NAMES = FALSE)

  myres$Info_UID <- ""
  myres$Info_center_pos <- 0
  for (i in 1:(length(pos) - 1)){
    idx <- (pos[i] + 1):(pos[i + 1] - 1)
    myres$Info_UID[idx] <- ids[i]
    myres$Info_center_pos[idx] <- seq_along(idx)
  }
  myres <- myres[-pos[1:(length(pos) - 1)], c(4, 5, 3)]
  names(myres) <- c("Info_UID", "Info_center_pos", "LBtope_prob")
  myres$LBtope_prob <- myres$LBtope_prob / 100
  myres$LBtope_class <- -1 + 2 * (myres$LBtope_prob > 0.5)

  return(myres)
}
