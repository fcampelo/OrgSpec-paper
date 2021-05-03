library(dplyr)
library(tidyr)
library(ranger)
library(epitopes)
library(ggplot2)

# Get test dirs
dirs <- dir("./", pattern = "[1-9]+_")

res <- vector("list", length(dirs))

for (i in seq_along(dirs)){
  # Load model trained for organism i
  res[[i]]$model <- readRDS(paste0(dirs[i], "/output/TrOrgSpec/model.rds"))
  res[[i]]$preds <- vector("list", length(dirs))
  
  for (j in seq_along(dirs)){
    cat("\nChecking model", i, "on dataset", j, "\n")
    # Load and prepare holdout set of proteins of organism j
    ho_prots <- readRDS(paste0(dirs[j],"/data/splits/prots_df.rds")) %>%
      as.data.frame() %>%
      dplyr::mutate(across(!starts_with("feat_"), as.character),
                    across(starts_with("feat_"), function(x) {x + 1e-9}),
                    Class = as.factor(Class))
    
    # Predict with model trained for organism i
    test.pred <- stats::predict(res[[i]]$model,
                                data = dplyr::select(ho_prots,
                                                     starts_with("feat_"),
                                                     Class))
    # Process predictions
    rf_probs <- test.pred$predictions[, "1"]
    rf_class <- ifelse(rf_probs >= 0.5, 1, -1)
    rf_preds <- ho_prots %>%
      dplyr::select(Info_UID, Info_center_pos, Class) %>%
      dplyr::bind_cols(pred_prob = rf_probs) %>%
      dplyr::mutate(Class      = as.numeric(as.character(Class)),
                    pred_class = as.numeric(as.character(rf_class)))
    
    # Calculate performance
    myperf <- epitopes::calc_performance(truth = ho_prots$Class,
                                         pred  = as.factor(rf_class),
                                         prob  = rf_probs,
                                         posValue = "1",
                                         negValue = "-1")
    
    res[[i]]$preds[[j]]$preds <- rf_preds
    res[[i]]$preds[[j]]$perf  <- myperf
    
    if (i == 1 && j == 1){
      df <- cbind(Model = dirs[i],
                  Data  = dirs[j],
                  res[[i]]$preds[[j]]$perf)
    } else {
      df <- rbind(df, 
                  cbind(Model = dirs[i],
                        Data  = dirs[j],
                        res[[i]]$preds[[j]]$perf))
    }
    
  }
}

df <- df %>%
  dplyr::select(Model, Data, 
                ACCURACY = accuracy, 
                AUC = auc, 
                MCC = mcc) %>%
  tidyr::pivot_longer(!(Model:Data), names_to = "Metric", values_to = "Value") %>%
  ungroup() %>%
  mutate(Model  = gsub("[0-9]+_", "", Model),
         Model.facets = paste0("Model: ", Model),
         Data   = gsub("[0-9]+_", "", Data),
         IsSame = ifelse(Model == Data, TRUE, NA))

ggplot(df, aes(x = Data, y = Value, 
               colour = !is.na(IsSame))) + 
  geom_point(size = 3, show.legend = FALSE) +
  geom_segment(aes(x = IsSame * .5, xend = IsSame * 3.5, 
                   y    = IsSame * Value,
                   yend = IsSame * Value),
               col = "#555555", lty = 2, alpha = 0.5,
               show.legend = FALSE) +
  geom_point(size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("#1b9e77", "#d95f02")) +
  scale_y_continuous(expand = c(.055,.05), breaks = seq(-1, 1, by = 0.2)) +
  facet_grid(Metric ~ Model.facets, scales = "free") + 
  theme_light() + 
  theme(strip.text  = element_text(colour = "black", face = "bold"),
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))

ggsave(plot = last_plot(), filename = "../figures/performance_cross.png",
       width = 6, height = 6, units = "in")
ggsave(plot = last_plot(), filename = "../figures/performance_cross.tiff",
       width = 6, height = 6, units = "in")

saveRDS(list(res = res, df = df), "../output/cross_performance.rds")
