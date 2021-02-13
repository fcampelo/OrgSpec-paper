# Install from github
# > devtools::install_github("fcampelo/epitopes", 
#                            ref = "v0.4.11-OrgSpec-paper")
#
# ============================================================================ #
#
# The 3 relevant data files for this script are:
# - "00_epitopes_20201006.rds"
# - "00_proteins_20201007.rds"
# - "00_taxonomy_20201007.rds"
#
# These were generated using:
#
# > epitopes::get_IEDB(save_folder = "XYZ")
# > epitopes <- epitopes::get_LBCE(data_folder = "XYZ",
# >                      ncpus = parallel::detectCores() - 1,
# >                      save_folder = "YZX")
# > proteins <- epitopes::get_proteins(uids = unique(epitopes$protein_id),
# >                      save_folder = "ZYX")
# > epitopes::get_taxonomy(uids = unique(epitopes$sourceOrg_id),
# >                      save_folder = "ZYX")
#
# This script also requires that BLAST be installed locally in the machine.

# Pathogen: Hep C
# Host: humans
OrgID   <- 11102
hostIDs <- c(9606, 10000113, 10000119, 10001673)

library(epitopes)
library(seqinr)
library(data.table)

epitopes      <- readRDS("../00_general_datasets/00_epitopes_20201006.rds")
proteins      <- readRDS("../00_general_datasets/00_proteins_20201007.rds")
tax_load_file <- "../00_general_datasets/00_taxonomy_20201007.rds"

# Join proteins into epitope dataset
jdf <- epitopes::prepare_join_df(epitopes, proteins,
                                 min_epit = 8, max_epit = 25,
                                 only_exact = FALSE,
                                 pos.mismatch.rm = "all", # If epitope is not where it should be in the protein sequence, remove.
                                 set.positive = "mode")   # In case of regions with conflicting labels epitope class is decided by simple majority vote.

# Filter data
jdf <- epitopes::filter_epitopes(jdf,
                                 orgIDs  = OrgID,
                                 hostIDs = hostIDs,
                                 tax_load_file = tax_load_file)

# Prepare windowed representation
wdf <- epitopes::make_window_df(jdf)

# Calculate features
wdf <- epitopes::calc_features(wdf, max.N = 2)

if(!dir.exists("./data/splits")) dir.create("./data/splits", recursive = TRUE)
saveRDS(wdf, "./data/splits/00full_set.rds")

prots <- proteins[UID %in% unique(jdf$protein_id), ]

# Run BLAST to determine splits
if(!dir.exists("./data/BLASTp")) dir.create("./data/BLASTp")
fn <- "./data/BLASTp/prots.fasta"
seqinr::write.fasta(as.list(prots$TSeq_sequence), names = prots$UID,
                    file.out = fn, as.string = TRUE, nbchar = 100000)

# Run BLASTp to determine protein similarities
system(paste0("makeblastdb -in ", fn, " -dbtype prot"))
system(paste0("blastp -query ", fn, " -db ", fn, " -seg no ",
              "-outfmt '6 qseqid sseqid length qlen slen nident pident qcovhsp mismatch gaps qstart qend sstart send evalue score' > ",
              fn, "-BLAST"))

# Split data into training/test/final validation sets, splitting by protein ID
# and using BLASTp to minimise similarities across splits
splits <- epitopes::split_epitope_data(wdf, split_level = "prot",
                                       split_perc = c(75, 25),
                                       split_names = c("01_training", "02_holdout"),
                                       save_folder = "./data/splits",
                                       blast_file = "./data/BLASTp/prots.fasta-BLAST",
                                       coverage_threshold = 80, identity_threshold = 80)

write.csv(splits[[1]]$wdf, "./data/splits/01_training.csv", row.names = FALSE)
saveRDS(splits[[1]]$wdf, "./data/splits/01_training.rds")
write.csv(splits[[2]]$wdf, "./data/splits/02_holdout.csv", row.names = FALSE)
saveRDS(splits[[2]]$wdf, "./data/splits/02_holdout.rds")