library(dplyr)
library(epitopes)
library(seqinr)

# Predictions table:
X     <- readRDS("./output/analysis.rds")
prots <- readRDS("../00_general_datasets/00_proteins_20201007.rds")

X$mypreds$Seq <- mapply(
  function(prot, st, en){
    substr(prots$TSeq_sequence[prots$UID == prot], st, en)
  },
  X$mypreds$Info_UID,
  X$mypreds$start_pos,
  X$mypreds$end_pos)

names(X$mypreds) <- c("Protein", "Start pos", "End pos", "Length", "Average Prob.", "Sequence")
X$mypreds[2:4] <- lapply(X$mypreds[2:4], as.integer)
X <- X$mypreds %>%
  dplyr::mutate(CandidateID = paste0("Pred", sprintf("%02d", seq_along(Protein))))

saveRDS(X, "./output/TrOrgSpec/pred_peptides.rds")
