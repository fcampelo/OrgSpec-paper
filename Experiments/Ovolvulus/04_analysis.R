# Pathogen: Onchocerca volvulus
OrgID   <- 6282

library(data.table)
library(reshape2)
library(dplyr)
library(epitopes)
library(ggplot2)
library(ggthemes)
library(xtable)
library(tidymodels)
library(tikzDevice)
source("../01_general_scrips/gather_results.R")
source("../01_general_scrips/gather_results_byPep.R")
source("../01_general_scrips/calc_perf.R")
source("../01_general_scrips/bootstrap_functions.R")
source("../01_general_scrips/make_plot2.R")

# Set PRNG seed
set.seed(20210107)

# Basic data
proteins <- readRDS("../00_general_datasets/00_proteins_20201007.rds")
epitopes <- readRDS("../00_general_datasets/00_epitopes_20201006.rds")
tax_path <- "../00_general_datasets/00_taxonomy_20201007.rds"

# Paths
mydata_path <- "./data/splits/holdout_prots_w.rds"
mypeps_path <- "./data/splits/02_holdout.rds"

res_paths   <- data.frame(
  name = c("RF_OrgSpec", "RF_Hybrid", "RF_Heter"),
  path = c("./output/TrOrgSpec/",
           "./output/TrHybrid/",
           "./output/TrHeter/"))

preds_paths <- data.frame(
  name     = c("ABCpred", "Bepipred2", "iBCE-EL", "LBtope", "SVMtrip"),
  path     = c("./output/ABCpred/",
               "./output/Bepipred2/",
               "./output/iBCE-EL/PIP-EL.html",
               "./output/lbtope/downloadresult.txt",
               "./output/SVMtrip/"),
  read.fun = c("../00_general_datasets/read_ABCPred.R",
               "../00_general_datasets/read_Bepipred2.R",
               "../00_general_datasets/read_iBCEEL.R",
               "../00_general_datasets/read_lbtope.R",
               "../00_general_datasets/read_SVMtrip.R"))
# =========

# Consolidate predictions by position
myres  <- gather_results(mydata_path = mydata_path,
                         res_paths = res_paths,
                         preds_paths = preds_paths)

# Consolidate predictions by peptide
myres_pep <- gather_results_byPep(mypeps_path = mypeps_path, 
                                  myres       = myres)
myperf_pep <- calc_perf(df          = myres_pep, 
                        res_names   = res_paths$name, 
                        preds_names = preds_paths$name)

# Prepare plots
mp_pep <- make_plot2(df = myperf_pep, linepos = 3.5, plot_rows = 2)

if(!dir.exists("./figures")) dir.create("./figures")
ggplot2::ggsave("./figures/res_Ov_byPep.png", plot = mp_pep, 
                width = 8, height = 5, units = "in")

tikz("./figures/res_Ov_byPep.tex",
     width = 8, height = 5)
mp_pep
dev.off()

# === hypothesis testing === #
# Reference method: RF (organism specific)
# Using results by peptide

# Estimate bootstrap distribution of differences UNDER THE NULL HYPOTHESES
ref.col   <- grep("RF_OrgSpec_class", names(myres_pep))
cmp.cols  <- grep("_class", names(myres_pep))
cmp.cols  <- cmp.cols[-which(cmp.cols == ref.col)]
diff_boot <- parallel::mclapply(1:10000,
                                FUN = boot_diffs,
                                df  = myres_pep,
                                ref.col  = ref.col,
                                cmp.cols = cmp.cols,
                                mc.cores = parallel::detectCores() - 1,
                                mc.set.seed = 20210201) %>%
  dplyr::bind_rows() %>%
  dplyr::rename_with(toupper)

# Get actual observed differences
diff_obs <- boot_diffs(NA, df = myres_pep, ref.col  = ref.col,
                       cmp.cols = cmp.cols) %>%
  dplyr::rename_with(toupper)

diff_boot <- rbind(diff_obs, diff_boot)

Pvals <- diff_boot %>%
  dplyr::group_by(METHOD1, METHOD2) %>%
  dplyr::summarise(across(!starts_with("METHOD"),
                          function(x){
                            (sum(abs(x) >= abs(dplyr::first(x))) - 1) / (n() - 1)}), 
                   .groups = "drop") %>%
  # Correct p-values for MHT using Benjamini-Hochberg
  dplyr::mutate(across(!starts_with("METHOD"),
                       ~p.adjust(.x, method = "BH")))

saveRDS(object = list(myres       = myres,
                      myres_pep   = myres_pep,
                      myperf_pep  = myperf_pep,
                      myplot_pep  = mp_pep,
                      Pvals_pep   = Pvals),
        file = "./output/analysis.rds")