# Assemble a data set of peptides from heterogeneous organisms
# (sampled from IEDB entries) for the tests

rn.seed <- 20210103

# Pathogen: Epstein-Barr virus
# Host: humans
# Org to remove: T. cruzi (to prevent skewing of the data)
target.Org     <- 10376
hostIDs        <- c(9606, 10000113, 10000119, 10001673)
Orgs.to.Remove <- 5693

# Desired number of entries for each class:
nPos <- 2000
nNeg <- 2000

save_folder   <- "./data/heterogeneous_data/"
orgData       <- readRDS("./data/splits/01_training.rds")

library(epitopes)
library(data.table)
library(dplyr)
source("../01_general_scrips/make_heterogeneous_set.R")

# Load epitope/protein data
proteins      <- readRDS("../00_general_datasets/00_proteins_20201007.rds")
epitopes      <- readRDS("../00_general_datasets/00_epitopes_20201006.rds")
tax_load_file <- "../00_general_datasets/00_taxonomy_20201007.rds"

wdf <- make_heterogeneous_set(hostIDs = hostIDs,
                              Orgs.to.Remove = c(target.Org, Orgs.to.Remove),
                              nPos = nPos, nNeg = nNeg,
                              proteins = proteins,
                              epitopes = epitopes,
                              tax_load_file = tax_load_file,
                              rnd.seed = rn.seed)


# Save to file
if(!dir.exists(save_folder)) dir.create(save_folder)
saveRDS(wdf, paste0(save_folder, "df_heterogeneous.rds"))

### Compose and save hybrid set
Hwdf <- rbind(wdf, orgData)
saveRDS(Hwdf, paste0(save_folder, "df_hybrid.rds"))