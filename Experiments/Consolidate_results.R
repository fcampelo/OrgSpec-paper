# Consolidate all results

library(dplyr)
library(ggplot2)
library(ggforce)
library(ggnewscale)
library(reshape2)
library(kableExtra)

if(!dir.exists("../figures")) dir.create("../figures")

# Method and performance metric labels (ordered, for plotting)
methods <- c("ABCpred", "Bepipred2", "iBCE-EL", "LBtope", "SVMtrip",
             "RF-Heter","RF-Hybrid","RF-OrgSpec")

metrics <- c("ACCURACY", "SENS", "PPV", "NPV", "AUC", "MCC")

# ============================================================================ #
# Consolidate results

dirs <- dir("./", pattern = "[1-9]+_")
preds <- vector("list", length(dirs))
names(preds) <- sapply(dirs, function(x)strsplit(x, split = "_")[[1]][2])

for (i in seq_along(dirs)){
  mydata <- readRDS(paste0(dirs[i], "/output/analysis.rds"))
  
  # Extract predictions data
  preds[[i]]$mypreds <- mydata$mypreds
  preds[[i]]$myprobs <- mydata$myprobs
  
  # Extract performance data
  tmp <- mydata$myperf_pep %>%
    dplyr::mutate(Organism = names(preds)[i]) %>%
    dplyr::filter(!(Metric %in% c("SPEC", "F1"))) # <-- not needed for the paper
  
  # Extract p-values
  tmp2 <- mydata$Pvals_pep %>%
    dplyr::select(-SPEC) %>%                 # <-- not needed for the paper
    dplyr::mutate(Organism = names(preds)[i])
  
  # Extract ROC curves
  cl.cols <- grep("_class", names(mydata$myres_pep))
  pr.cols <- grep("_prob", names(mydata$myres_pep))
  for (j in seq_along(cl.cols)){
    myperf <- epitopes::calc_performance(truth = mydata$myres_pep$Class,
                                         pred = mydata$myres_pep[, cl.cols[j]],
                                         prob = mydata$myres_pep[, pr.cols[j]],
                                         ret.as.list = TRUE)
    if (j == 1){
      tmp3 <- data.frame(Organism = names(preds)[i],
                         Method = gsub("_class", "", names(mydata$myres_pep)[cl.cols[j]]),
                         TPR    = myperf$tpr,
                         FPR    = myperf$fpr)
    } else {
      tmp3 <- rbind(tmp3,
                    data.frame(Organism = names(preds)[i],
                               Method = gsub("_class", "", names(mydata$myres_pep)[cl.cols[j]]),
                               TPR    = myperf$tpr,
                               FPR    = myperf$fpr))
    }
  }
  
  tmp4 <- mydata$Pvals_pep.raw %>%
    dplyr::select(-SPEC) %>%                 # <-- not needed for the paper
    dplyr::mutate(Organism = names(preds)[i])
  
  if (i == 1){
    df    <- tmp
    pvals <- tmp2
    rocs  <- tmp3
    pvals.raw <- tmp4
  } else {
    df    <- rbind(df, tmp)
    pvals <- rbind(pvals, tmp2)
    rocs  <- rbind(rocs, tmp3)
    pvals.raw <- rbind(pvals.raw, tmp4)
  }
}

# prepare data frame with p-values
pvals <- pvals %>%
  ungroup() %>%
  dplyr::select(-METHOD1) %>%
  dplyr::rename(Method = METHOD2) %>%
  reshape2::melt(id.vars = c("Method", "Organism", "NoLeak"),
                 variable.name = "Metric",
                 value.name = "pValue") %>%
  dplyr::mutate(Method = gsub("_", "-", Method, fixed = TRUE),
                pValue = round(pValue, 3))

# Add p-values to results df.
df <- df %>%
  dplyr::mutate(Method = as.character(Method),
                Metric = as.character(Metric)) %>%
  dplyr::left_join(pvals, by = c("Organism", "Method", "NoLeak", "Metric")) %>%
  dplyr::mutate(Method = factor(Method, levels = methods, ordered = TRUE),
                Metric = factor(Metric, levels = metrics, ordered = TRUE),
                pValue = pValue) %>%
  group_by(Metric, Organism) %>%
  mutate(Result    = ifelse(Mean >= first(Mean), "better", "worse"))

# Other cosmetic changes for better plotting
df$Result[df$pValue > 0.05] <- "non-signif"
df$Result[df$Method == "RF-OrgSpec"] <- "ref. method"
df$Result <- factor(df$Result,
                    levels = c("ref. method", "non-signif", "better", "worse"),
                    ordered = TRUE)
df$pValue[df$pValue <= .01] <- paste0("<0.01")
df$pValue[df$pValue >= 0.9] <- paste0(">0.9")


# ============================================================================ #
# Generate plots with performance metrics and comparisons

# Plot ACC, AUC, MCC
mp <- ggplot(filter(df, NoLeak == FALSE, Metric %in% c("ACCURACY", "AUC", "MCC")),
             aes(x = Method, y = Value,
                 ymin = Value - StdErr, ymax = Value + StdErr,
                 colour  = Result)) +
  geom_text(aes(label = pValue), size = 2, col = "#222222", nudge_x = -0.4) +
  geom_pointrange(fatten = 1.5, size = .75, show.legend = FALSE) +
  scale_color_manual(values = c("#555555", "#7570c3", "#1b9e77", "#d95f02")) +
  ylab("") + xlab("") +
  facet_grid(Organism ~ Metric, scales = "free") +
  geom_vline(xintercept = 5.35, size = .2, lty = 2) +
  scale_y_continuous(expand = c(.055,.05), breaks = seq(-1, 1, by = 0.1)) +
  scale_x_discrete(expand = c(.1,.05)) +
  coord_flip() + theme_light() +
  theme(strip.text  = element_text(colour = "black", face = "bold"),
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        plot.margin = margin(2, 2, 2, 0))

ggsave(plot = mp, filename = "../figures/res_all_pathogens_1.png",
       width = 7.5, height = 4.8, units = "in")
ggsave(plot = mp, filename = "../figures/res_all_pathogens_1.tiff",
       width = 7.5, height = 4.8, units = "in")

# Plot SENS, PPV, NPV
mp <- mp %+% dplyr::filter(df, NoLeak == FALSE,
                           Metric %in% c("SENS", "PPV", "NPV")) +
  ylab("Estimated performance") +
  scale_color_manual(values = c("#555555", "#7570c3", "#d95f02", "#1b9e77"))

ggsave(plot = mp, filename = "../figures/res_all_pathogens_2.png",
       width = 7.5, height = 4.8, units = "in")
ggsave(plot = mp, filename = "../figures/res_all_pathogens_2.tiff",
       width = 7.5, height = 4.8, units = "in")

# ============================================================================ #
# Summary results table
x <- df %>%
  select(Organism, Method, Metric, Value, StdErr, NoLeak) %>%
  mutate(Value  = round(Value, 2), 
         StdErr = round(StdErr, 2),
         Field  = paste0("SSS", Value, "PM ", StdErr, "SSS")) %>%
  group_by(Organism, Method, Metric) %>%
  summarise(Field = ifelse(Organism == "HepC",
                           paste0(first(Field), " (", last(Field), ")"),
                           as.character(first(Field)))) %>%
  ungroup() %>%
  distinct() %>%
  tidyr::pivot_wider(names_from = c(Metric),
                     values_from = c(Field))

# ============================================================================ #
# Plot ROC curves
rocs$Method <- factor(gsub("_", "-", rocs$Method), levels = methods, ordered = TRUE)

mp <- ggplot(rocs, aes(x = FPR, y = TPR, colour = Method)) +
  geom_line() + geom_abline(slope = 1, lty = 2) + theme_light() +
  scale_colour_brewer(type = "div") +
  facet_grid(. ~ Organism, scales = "free") +
  theme_light() +
  theme(strip.text  = element_text(colour = "black", face = "bold"),
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        legend.background = element_rect(fill = "#ffffff00",
                                         colour = "#444444"))

ggsave(plot = mp, filename = "../figures/ROC_all.png",
       width = 8, height = 4, units = "in")
ggsave(plot = mp, filename = "../figures/ROC_all.tiff",
       width = 8, height = 4, units = "in")

# O. volvulus ROC
mp <- ggplot(dplyr::filter(rocs, Organism == "Ovolvulus"),
             aes(x = FPR, y = TPR, colour = Method)) +
  geom_line() + geom_abline(slope = 1, lty = 2) + theme_light() +
  scale_colour_brewer(type = "div") +
  theme_light() +
  theme(strip.text  = element_text(colour = "black", face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = c(.85,.3),
        legend.background = element_rect(fill = "#ffffff00",
                                         colour = "#444444")) +
  guides(colour = guide_legend(override.aes = list(fill = NA)))

ggsave(plot = mp, filename = "../figures/Ov_ROC.png",
       width = 5, height = 4.5, units = "in")
ggsave(plot = mp, filename = "../figures/Ov_ROC.tiff",
       width = 5, height = 4.5, units = "in")


# ============================================================================ #
# Generate protein prediction map examples for paper
for (i in seq_along(preds)){
  tmp <- preds[[i]]$myprobs %>%
    dplyr::mutate(Prediction = (1 + RF_OrgSpec_class) / 2,
                  Class      = (1 + Class) / 2,
                  PosPred    = Prediction / Prediction,
                  IsNew      = factor(Prediction & Class,
                                      levels = c(TRUE, FALSE, NA),
                                      labels = c("Known epitope", NA, "New epitope"),
                                      exclude = NULL)) %>%
    dplyr::select(-RF_OrgSpec_class)
  
  
  mp <- tmp %>%
    ggplot(aes(x = Info_center_pos,
               y = Class)) +
    geom_point(size = 2, pch = 73, col = "#99dd99") +
    geom_line(aes(y = Prediction),
              lwd = .1, show.legend = FALSE) +
    geom_point(data = subset(tmp, !is.na(IsNew)),
               aes(y = PosPred, colour = IsNew),
               size = .6, pch = 15) +
    scale_colour_manual(values = c("#118811", "#992222")) +
    new_scale_color() +
    geom_line(aes(y = RF_OrgSpec_prob, colour = RF_OrgSpec_prob),
              lwd = .1, show.legend = FALSE) +
    scale_colour_gradient(low = "#55bbbb",high = "#cc8888") +
    scale_y_continuous(breaks = c(0, .5, 1), limits = c(0, 1.05)) +
    ylab("Prediction") + xlab("") +
    theme_light() +
    theme(strip.text  = element_text(colour = "black", face = "bold"),
          axis.text.x = element_text(size = 8),
          plot.margin = margin(3, 3, 3, 3),
          legend.title = element_blank(),
          legend.position = "none",
          legend.background = element_rect(colour = "#aaaaaa"))
  
  for (j in seq(1, ceiling(length(unique(tmp$Info_UID)) / 12))){
    x <- mp + facet_wrap_paginate(Info_UID ~ .,
                                  scales = "free",
                                  nrow = 6, ncol = 2, page = j)
    ggsave(plot = x,
           filename = paste0("../figures/prot_", dirs[i], "-", j, ".png"),
           width = 10, height = 12, units = "in")
    ggsave(plot = x,
           filename = paste0("../figures/prot_", dirs[i], "-", j, ".tiff"),
           width = 10, height = 12, units = "in")
  }
}



# Raw p-values table
pvals.raw%>%
  ungroup() %>%
  dplyr::filter(NoLeak == FALSE) %>%
  dplyr::mutate(Comparison = paste0(METHOD1, " vs. ", METHOD2),
                Comparison = gsub("_", "-", Comparison, fixed = TRUE)) %>%
  dplyr::select(Organism, Comparison, everything(), -METHOD1, -METHOD2, -NoLeak) %>%
  dplyr::mutate(across(where(is.numeric), ~round(.x, digits = 3))) %>%
  kableExtra::kable(format = "latex", longtable = TRUE)