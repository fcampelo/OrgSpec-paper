# Assemble data sets for experiment:
# - Organism Specific: Train and Test
# - Heterogeneous and Hybrid

target.Org   <- 6282     # O. volvulus
rnd.seed     <- 20210103 # Random seed
nPos <- nNeg <- 3000     # For the heterogeneous dataset
ncpus        <- 4        # For parallel execution

hostIDs <- c(9606, 10000113, 10000119, 10001673) # Human IEDB IDs

# Ignore epitopes coming from T. cruzi (to prevent skewing of the data), plants
# and vertebrates when composing the heterogeneous dataset
removeIDs <- c(5693, 33090, 7742)

# Basic datasets
epitopes <- readRDS("../00_general_datasets/00_epitopes_20201006.rds")
proteins <- readRDS("../00_general_datasets/00_proteins_20201007.rds")
taxonomy <- readRDS("../00_general_datasets/00_taxonomy_20201007.rds")

# Generate datasets
source("../00_general_scrips/generate_datasets.R")
generate_datasets(target.Org, hostIDs, removeIDs, rnd.seed,
                  nPos, nNeg,
                  epitopes, proteins, taxonomy,
                  coverage_threshold = 80, identity_threshold = 80,
                  ncpus = ncpus)
