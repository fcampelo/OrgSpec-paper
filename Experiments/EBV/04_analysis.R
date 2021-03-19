# Pathogen: Epstein-Barr Virus
OrgID   <- 10376

library(data.table)
library(reshape2)
library(dplyr)
library(epitopes)
library(ggplot2)
library(ggthemes)
library(tikzDevice)
source("../01_general_scrips/gather_results.R")
source("../01_general_scrips/gather_results_byPep.R")
source("../01_general_scrips/calc_perf.R")
source("../01_general_scrips/bootstrap_functions.R")
source("../01_general_scrips/make_plot3.R")
source("../01_general_scrips/sort_predictions.R")

# Set PRNG seed
set.seed(20210107)

# Get training data from each predictor
epits     <- readRDS("../00_general_datasets/00_epitopes_20201006.rds")
tr.epits  <- unique(readRDS("./data/splits/01_training.rds")$Info_epitope_id)
pred.tr.sets <- readRDS("../../predictors_training_data/predictor_training_seqs.rds")
pred.tr.sets$RF_OrgSpec <- epits[epitope_id %in% tr.epits, .(Sequence = epit_seq)]
pred.tr.sets$RF_Hybrid  <- pred.tr.sets$RF_OrgSpec
pred.tr.sets$RF_Heter   <- data.frame(Sequence = "1234")

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
  read.fun = c("../01_general_scrips/read_ABCPred.R",
               "../01_general_scrips/read_Bepipred2.R",
               "../01_general_scrips/read_iBCEEL.R",
               "../01_general_scrips/read_lbtope.R",
               "../01_general_scrips/read_SVMtrip.R"))
# =========

# Consolidate predictions by position
myres  <- gather_results(mydata_path = mydata_path,
                         res_paths = res_paths,
                         preds_paths = preds_paths)

# Consolidate predictions by peptide
myres_pep <- gather_results_byPep(mypeps_path = mypeps_path,
                                  myres       = myres)

df <- myres_pep %>%
  left_join(epits[, .(Info_epitope_id = epitope_id, epit_seq)], 
            by = "Info_epitope_id") %>%
  dplyr::distinct()

myperf_pep <- calc_perf(df           = df,
                        res_names    = res_paths$name,
                        preds_names  = preds_paths$name,
                        pred.tr.sets = pred.tr.sets)

# Prepare plot
mp_pep <- make_plot3(df = myperf_pep, linepos = 3.5, plot_rows = 2)

if(!dir.exists("./figures")) dir.create("./figures")
ggplot2::ggsave("./figures/res_Ov_byPep.png", plot = mp_pep,
                width = 8, height = 5, units = "in")

tikz("./figures/res_Ov_byPep.tex",
     width = 8, height = 5)
mp_pep
dev.off()

# === hypothesis testing === #
# 1) Considering full holdout set
# Reference method: RF (organism specific)
# Using results by peptide

# Estimate bootstrap distribution of differences UNDER THE NULL HYPOTHESES
ref <- "RF_OrgSpec"
nBoot.pval <- 1000
diff_boot <- parallel::mclapply(1:nBoot.pval,
                                FUN = boot_diffs,
                                df  = df,
                                ref  = ref,
                                mc.cores = parallel::detectCores() - 1,
                                mc.set.seed = 20210201) %>%
  dplyr::bind_rows() %>%
  dplyr::rename_with(toupper)

# Get actual observed differences
diff_obs <- boot_diffs(NA, df = myres_pep, ref = ref) %>%
  dplyr::rename_with(toupper)

diff_boot <- rbind(diff_obs, diff_boot)

Pvals <- diff_boot %>%
  dplyr::group_by(METHOD1, METHOD2) %>%
  dplyr::summarise(across(!starts_with("METHOD"),
                          function(x){
                            (sum(abs(x) >= abs(dplyr::first(x)))) / (n() - 1)}),
                   .groups = "drop") %>%
  # Correct p-values for MHT using Benjamini-Hochberg
  dplyr::mutate(across(!starts_with("METHOD"),
                       ~p.adjust(.x, method = "BH")),
                NoLeak = FALSE)


# 2) Considering only zero-leakage subsets
noleak_res <- lapply(pred.tr.sets,
                     function(X){
                       filter(df, !(epit_seq %in% X$Sequence))})

# Estimate bootstrap distribution of differences UNDER THE NULL HYPOTHESES
diff_boot <- parallel::mclapply(1:nBoot.pval,
                                FUN = boot_diffs_nopair,
                                df  = noleak_res,
                                ref  = ref,
                                mc.cores = parallel::detectCores() - 1,
                                mc.set.seed = 20210201) %>%
  dplyr::bind_rows() %>%
  dplyr::rename_with(toupper)

# Get actual observed differences
diff_obs <- boot_diffs_nopair(NA, df = noleak_res, ref = ref) %>%
  dplyr::rename_with(toupper)

diff_boot <- rbind(diff_obs, diff_boot)

Pvals_noleak <- diff_boot %>%
  dplyr::group_by(METHOD1, METHOD2) %>%
  dplyr::summarise(across(!starts_with("METHOD"),
                          function(x){
                            (sum(abs(x) >= abs(dplyr::first(x)))) / (n() - 1)}),
                   .groups = "drop") %>%
  # Correct p-values for MHT using Benjamini-Hochberg
  dplyr::mutate(across(!starts_with("METHOD"),
                       ~p.adjust(.x, method = "BH")),
                NoLeak = TRUE)

Pvals <- rbind(Pvals, Pvals_noleak)

predlist <- sort_predictions(myres)

saveRDS(object = c(list(myres       = myres,
                        myres_pep   = myres_pep,
                        myperf_pep  = myperf_pep,
                        myplot_pep  = mp_pep,
                        Pvals_pep   = Pvals,
                        nBoot.pval  = nBoot.pval),
                   predlist),
        file = "./output/analysis.rds")
