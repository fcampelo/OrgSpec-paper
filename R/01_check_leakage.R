library(dplyr)
library(reshape2)
library(seqinr)

epits <- readRDS("../Experiments/00_general_datasets/00_epitopes_20201006.rds")
prots <- readRDS("../Experiments/00_general_datasets/00_proteins_20201007.rds")

# retrieve and consolidate training datasets for the literature predictors
## 1) ABCPred
con <- file("../predictors_training_data/data_abcpred/redundant_bcipep", 
            open = "r")

abcpred.data <- data.frame(Field = character(), Value = character())

while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  myVector <- unlist(strsplit(oneLine, "[[:space:]]"))
  myVector <- gsub('\"', "", myVector)
  if (tolower(myVector[1]) %in% c("sequence", "immunogenicity")){
    tmp <- data.frame(Field = tolower(myVector[[1]]), 
                      Value = myVector[[length(myVector)]])
    if (nrow(abcpred.data) == 0){
      abcpred.data <- tmp
    } else {
      abcpred.data <- rbind(abcpred.data, tmp)
    }
  }
}
close(con)

abcpred.data <- data.frame(
  Sequence = abcpred.data$Value[abcpred.data$Field == "sequence"],
  Class    = abcpred.data$Value[abcpred.data$Field == "immunogenicity"])

abcpred.data$Class[abcpred.data$Class != "No"] <- "Yes"


## 1) Bepipred2
bepipred_ids <- seqinr::read.fasta("../predictors_training_data/data_bepipred2/iedb_linear_epitopes.fasta")
bepipred_ids <- unlist(lapply(bepipred_ids, 
                              function(x){paste(toupper(x), collapse = "")}))
bepipred_ids <- data.frame(Sequence = unname(bepipred_ids),
                           IDstring = names(bepipred_ids))
bepipred_ids$Class <- sapply(strsplit(bepipred_ids$IDstring, split = "_"),
                             function(x){gsub("ID", "", x[[1]])})
bepipred_ids$ID <- sapply(strsplit(bepipred_ids$IDstring, split = "_"),
                          function(x){x[[2]]})

bepipred.data <- bepipred_ids %>%
  dplyr::left_join(epits, by = c("ID" = "epitope_id")) %>%
  dplyr::select(Sequence = epit_seq, Class)


## 3) iBCE-EL
ibceel_ids <- c(seqinr::read.fasta("../predictors_training_data/data_ibceel/B-positive.txt"),
                seqinr::read.fasta("../predictors_training_data/data_ibceel/B-negative.txt"))
ibceel_ids <- unlist(lapply(ibceel_ids, 
                            function(x){paste(toupper(x), collapse = "")}))
ibceel.data <- data.frame(Sequence = unname(ibceel_ids),
                          IDstring = names(ibceel_ids)) %>%
  dplyr::mutate(Class = sapply(strsplit(IDstring, split = "_"),
                               function(x){gsub("ID", "", x[[1]])})) %>%
  dplyr::select(-IDstring)


## 4) LBtope
lbtope.data <- read.csv("../predictors_training_data/data_lbtope/LBtope_Variable_Positive_epitopes.txt", 
                        header = FALSE) %>%
  dplyr::rename(Sequence = V1) %>%
  dplyr::mutate(Class = NA)


## 5) SVMtrip
svmtrip.data <- read.csv("../predictors_training_data/data_svmtrip/negativeset/negativeset_20AA", 
                         header = FALSE) %>%
  dplyr::rename(Sequence = V1) %>%
  dplyr::mutate(Class = "Negative")

svmtrip.data <- rbind(svmtrip.data,
                      read.csv("../predictors_training_data/data_svmtrip/posistiveset/positiveset_20AA", 
                               header = FALSE) %>%
                        dplyr::rename(Sequence = V1) %>%
                        dplyr::mutate(Class = "Positive"))


# Check for the presence of epitopes in the holdout sets of each organism
# in the training sets of the existing methods
data.list <- list(abcpred  = abcpred.data,
                  bepipred = bepipred.data,
                  ibcel    = ibceel.data,
                  lbtope   = lbtope.data,
                  svmtrip  = svmtrip.data)


orgs <- dir("../Experiments/")
idx  <- which(!grepl("0", orgs))
org.list <- data.frame(Name = orgs[idx],
                       path = dir("../Experiments", full.names = TRUE)[idx])

for (i in 1:nrow(org.list)){
  cat("\n Processing data for", org.list$Name[i])
  tr.data  <- read.csv(paste0(org.list$path[i], "/data/splits/01_training.csv"))
  tr.ids   <- unique(tr.data$Info_epitope_id)
  tr.seqs  <- epits$epit_seq[epits$epitope_id %in% tr.ids]
    
  ho.data  <- read.csv(paste0(org.list$path[i], "/data/splits/02_holdout.csv"))
  ho.ids   <- unique(ho.data$Info_epitope_id)
  ho.seqs  <- epits$epit_seq[epits$epitope_id %in% ho.ids]
  ho.prots <- prots$TSeq_sequence[prots$UID %in% unique(ho.data$Info_protein_id)]

  leak_count <- lapply(data.list,
                       function(x){length(which(ho.seqs %in% x$Sequence))})

  # Count differently for SVMtrip (due to windowing of training set)
  leak_count$svmtrip <- length(unlist(lapply(data.list$svmtrip$Sequence,
                                             function(x){grep(x, ho.prots)})))

  orgspec_leak <- length(which(ho.seqs %in% tr.seqs))

  if (i == 1){
    df <- cbind(Org = org.list$Name[i],
                holdout_size = length(ho.seqs),
                as.data.frame(leak_count) / length(ho.seqs),
                OrgSpec = orgspec_leak / length(ho.seqs))
  } else {
    df <- rbind(df,
                cbind(Org = org.list$Name[i],
                      holdout_size = length(ho.seqs),
                      as.data.frame(leak_count) / length(ho.seqs),
                      OrgSpec = orgspec_leak / length(ho.seqs)))
  }
  cat("\n")
  print(df)
}
