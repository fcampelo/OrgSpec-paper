# Run experiment (train models using organism-specific, heterogeneous and
# hybrid data sets, evaluate performance on hold-out set, calculate performance
# indices and generate individual plots)

OrgID    <- 1314      # Pathogen: S. pyogenes
rnd.seed <- 20210429
ncpus    <- 7

# Load training and holdout data
epitopes        <- readRDS("../00_general_datasets/00_epitopes_20201006.rds")
data.train      <- readRDS("./data/splits/01_training.rds")
pred.train.sets <- readRDS("../../predictors_training_data/predictor_training_seqs.rds")
ho_prots        <- readRDS("./data/splits/prots_df.rds")
ho_peps         <- readRDS("./data/splits/02_holdout.rds")

invisible(sapply(dir("../00_general_scripts/",
                     pattern = ".R",
                     full.names = TRUE),
                 source))

# ============================================================================ #
# Train models and get hold-out predictions

# Training sets
Tr_specs <- data.frame(
  TrData  = c("./data/splits/01_training.rds",
              "./data/heterogeneous_data/df_heterogeneous.rds",
              "./data/heterogeneous_data/df_hybrid.rds"),
  out_dir = c("./output/TrOrgSpec",
              "./output/TrHeter",
              "./output/TrHybrid"))

run_experiment(rnd.seed, Tr_specs, ho_prots, ncpus)

# ============================================================================ #
# Perform analysis

# Results paths
myres_files <- data.frame(name = c("RF_OrgSpec", "RF_Hybrid", "RF_Heter"),
                          file = c("./output/TrOrgSpec/preds.rds",
                                   "./output/TrHybrid/preds.rds",
                                   "./output/TrHeter/preds.rds"))

preds_paths <- data.frame(name     = c("ABCpred", "Bepipred2", "iBCE-EL",
                                       "LBtope", "SVMtrip"),
                          path     = c("./output/ABCpred/",
                                       "./output/Bepipred2/",
                                       "./output/iBCE-EL/",
                                       "./output/lbtope/",
                                       "./output/SVMtrip/"),
                          read.fun = c("read_abcpred", "read_bepipred2",
                                       "read_ibceel", "read_lbtope",
                                       "read_svmtrip"))

run_analysis(rnd.seed, epitopes, data.train, pred.train.sets,
             ho_prots, ho_peps, myres_files, preds_paths, ncpus)


