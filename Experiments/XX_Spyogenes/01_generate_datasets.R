# Assemble data sets for experiment:
# - Organism Specific: Train and Test
# - Heterogeneous and Hybrid

target.Org   <- 1314     # S. pyogenes
rnd.seed     <- 20210103 # Random seed
nPos <- nNeg <- 2000     # For the heterogeneous dataset
ncpus        <- 4        # For parallel execution

hostIDs <- c(9606, 10000113, 10000119, 10001673) # Human IEDB IDs

# Ignore epitopes coming from T. cruzi (to prevent skewing of the data), plants
# and vertebrates when composing the heterogeneous dataset
removeIDs <- c(5693, 33090, 7742)

# Basic datasets
epitopes <- readRDS("../00_general_datasets/00_epitopes_20201006.rds")
proteins <- readRDS("../00_general_datasets/00_proteins_20201007.rds")
taxonomy <- readRDS("../00_general_datasets/00_taxonomy_20201007.rds")

### ============== DATA REBALANCING BY SUB-SAMPLING ============== ###
library(dplyr)
library(epitopes)
jdf <- prepare_join_df(epitopes, proteins, only_exact = FALSE)
jdf <- filter_epitopes(jdf, orgIDs = target.Org, hostIDs = hostIDs, 
                       tax_list = taxonomy)

# Original data balance: 0.9607 x 0.0392
table(jdf$Class)
table(jdf$Class)/sum(table(jdf$Class))

# Find proteins containing only (or mostly) negative examples
x <- jdf %>% 
  dplyr::group_by(protein_id) %>% 
  dplyr::summarise(neg = sum(Class == -1), pos = sum(Class == 1)) %>%
  dplyr::arrange(desc((neg+1)/(pos+1)))
print(x, n = 20)

# Remove observations from proteins containing only negative examples 
ids <- x[x$pos == 0, ]$protein_id
epitopes <- epitopes[-which(epitopes$protein_id %in% ids), ]
x <- x[x$pos > 0, ]
print(x, n = 20)

# Manually remove observations from proteins with very high neg/pos ratios, 
# as long as they don't contain many positive observations.
torm <- x$protein_id[c(1:13, 16:18)]
epitopes <- epitopes[-which(epitopes$protein_id %in% torm), ]
x <- x[!(x$protein_id %in% torm), ]

# Final data balance: 0.7542 x 0.2458
colSums(x[, -1])
colSums(x[, -1]) / sum(colSums(x[, -1]))

### ============================================================== ###

# Generate datasets
source("../00_general_scripts/generate_datasets.R")
generate_datasets(target.Org, hostIDs, removeIDs, rnd.seed,
                  nPos, nNeg,
                  epitopes, proteins, taxonomy,
                  coverage_threshold = 80, identity_threshold = 80,
                  ncpus = ncpus)
