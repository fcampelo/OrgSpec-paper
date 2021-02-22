read_iBCEEL <- function(res.path, ...){
  require(XML)

  myres  <- XML::readHTMLTable(res.path, header = TRUE)[[1]]
  idvars <- strsplit(myres$`FASTA ID`, split = "pp")

  myres$Info_UID        <- sapply(idvars, function(x) x[1])
  myres$Info_UID        <- gsub("-", "_", myres$Info_UID, fixed = TRUE)
  myres$Info_center_pos <- sapply(idvars, function(x) x[2])
  myres$`iBCE-EL_prob`  <- pmax(0, pmin(1, myres$Prob))
  myres$`iBCE-EL_class` <- -1 + 2 * as.numeric(myres$`PIP or Non-PIP` == "BCE")

  myres$Info_center_pos <- as.numeric(myres$Info_center_pos)

  return(myres[, c("Info_UID", "Info_center_pos", "iBCE-EL_prob", "iBCE-EL_class")])
}
