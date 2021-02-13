# Pathogen: Hep C
# Host: humans
OrgID   <- 11102
hostIDs <- c(9606, 10000113, 10000119, 10001673)

# Load holdout data
wdf <- readRDS("./data/splits/02_holdout.rds")

# Save folder
save_folder <- "./data/splits/"


library(epitopes)
library(dplyr)
library(data.table)
library(seqinr)


# Load basic epitope/protein data as well as split 3
proteins      <- readRDS("../00_general_datasets/00_proteins_20201007.rds")
epitopes      <- readRDS("../00_general_datasets/00_epitopes_20201006.rds")
tax_load_file <- "../00_general_datasets/00_taxonomy_20201007.rds"


# Assemble filtered dataset
jdf <- epitopes::prepare_join_df(epitopes, proteins) %>%
  epitopes::filter_epitopes(orgIDs = OrgID,
                            hostIDs = hostIDs,
                            tax_load_file = tax_load_file)

# Extract holdout proteins
res_prot <- proteins[UID %in% wdf$Info_protein_id]

# Windowed representation of holdout proteins
wres_prot <- epitopes::make_window_df(res_prot,
                                      window_size = nchar(wdf$Info_window_seq[1]))

# Extract labels
lres_prot <- epitopes::label_proteins(res_prot, jdf, set.positive = "mode")
lres_prot <- lres_prot[, list(Info_UID, Info_center_pos, Class)]

# Join labels and calculate features
wres_prot <- wres_prot %>%
  left_join(lres_prot, by = c("Info_UID", "Info_center_pos")) %>%
  epitopes::calc_features(max.N = 2)


# Save reserve set proteins as FASTA and CSV
write.csv(wres_prot,
          paste0(save_folder, "/holdout_prots_w.csv"),
          row.names = FALSE)
saveRDS(wres_prot, paste0(save_folder, "/holdout_prots_w.rds"))
seqinr::write.fasta(as.list(gsub("[BJXZ]", "", wres_prot$Info_window_seq)),
                    names = paste(wres_prot$Info_UID, wres_prot$Info_center_pos, sep = "pp"),
                    file.out = paste0(save_folder, "/holdout_prots_w.fasta"),
                    as.string = TRUE)

write.csv(res_prot[, c("TSeq_accver","TSeq_taxid","TSeq_orgname","TSeq_sequence","UID")],
          paste0(save_folder, "/holdout_prots.csv"),
          row.names = FALSE)
seqinr::write.fasta(as.list(res_prot$TSeq_sequence), names = res_prot$UID,
                    file.out = paste0(save_folder, "/holdout_prots.fasta"),
                    as.string = TRUE, nbchar = 10000)

