# Consolidate all results

library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggnewscale)
library(ggforce)
library(reshape2)
library(tikzDevice)

dirs <- c("EBV", "HepC", "Ovolvulus", "Pvivax")
preds <- vector("list", length(dirs))
names(preds) <- dirs

for (i in seq_along(dirs)){
  mydata     <- readRDS(paste0("../Experiments/", 
                               dirs[i], "/output/analysis.rds")) 
  
  # Extract prediction data
  preds[[i]]$mypreds <- mydata$mypreds
  preds[[i]]$myprobs <- mydata$myprobs
  
  # Extract performance data
  tmp <- mydata$myperf_pep %>%
    dplyr::mutate(Organism = dirs[i]) %>%
    dplyr::filter(!(Metric %in% c("SENS", "SPEC")))
  
  # Extract p-values
  tmp2 <- mydata$Pvals_pep %>%
    dplyr::mutate(Organism = dirs[i])
  
  if (i == 1){
    df    <- tmp
    pvals <- tmp2
  } else {
    df    <- rbind(df, tmp)
    pvals <- rbind(pvals, tmp2)
  }
}

# prepare p-value data frame
pvals <- pvals %>%
  ungroup() %>%
  dplyr::select(-METHOD1) %>%
  dplyr::rename(Method = METHOD2) %>%
  reshape2::melt(id.vars = c("Method", "Organism"), 
                 variable.name = "Metric", 
                 value.name = "pValue") %>%
  dplyr::mutate(Method = gsub("_", "-", Method, fixed = TRUE),
                pValue = signif(pValue, 2))

# Method and perf. metric labels for plotting
methods <- c("RF-OrgSpec","RF-Hybrid","RF-Heter",
             "ABCpred", "Bepipred2", "iBCE-EL", "LBtope", "SVMtrip")

metrics <- c("ACCURACY", "PPV", "NPV", "MCC")

# Add p-values to results df.
df <- df %>%
  dplyr::mutate(Method = as.character(Method),
                Metric = as.character(Metric)) %>%
  dplyr::left_join(pvals, by = c("Organism", "Method", "Metric")) %>%
  dplyr::mutate(Method = factor(Method, levels = methods, ordered = TRUE),
                Metric = factor(Metric, levels = metrics, ordered = TRUE),
                pValue = pValue) %>%
  group_by(Metric) %>%
  mutate(nudge_dir = ifelse(Mean < mean(Mean), 1, -1)) %>%
  group_by(Metric, Organism) %>%
  mutate(Result = ifelse(Mean >= first(Mean), "better", "worse"))

df$Result[df$pValue > 0.05] <- "non-signif"
df$Result[df$Method == "RF-OrgSpec"] <- "ref. method"
df$Result <- factor(df$Result, 
                    levels = c("ref. method", "non-signif", "better", "worse"),
                    ordered = TRUE)
df$pValTex <- df$pValue
df$pValTex[df$pValTex == 0] <- "$<0.0001$"
df$pValue[df$pValue == 0] <- "<0.0001"

# Define colour map
mycols <- c("#555555", "#7570c3", "#d95f02", "#1b9e77")

mp <- ggplot(df, aes(x       = Method,
                     y       = Value, 
                     ymin    = Value - StdErr, 
                     ymax    = Value + StdErr,
                     colour  = Result)) +
  geom_pointrange(fatten = 2, size = 1, show.legend = FALSE) +
  scale_color_manual(values = mycols) +
  geom_text(aes(label = pValue), size = 2.5, col = "#222222",
            nudge_x = -0.4, nudge_y = df$nudge_dir * 0.06) +
  ylab("Estimated performance") + xlab("") +
  facet_grid(Organism ~ Metric, scales = "free") +
  geom_vline(xintercept = 3.35, size = .2, lty = 2) + 
  coord_flip() + 
  theme_light() +
  theme(strip.text  = element_text(colour = "black", face = "bold"),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        plot.margin = margin(2, 2, 2, 0))

ggsave(plot = mp, filename = "../figures/res_all_pathogens.png", 
       width = 10, height = 10, units = "in")

# Save as tikz figure
tikz("../figures/res_all_pathogens.tex",
     width = 7.5, height = 9)
ggplot(df, aes(x       = Method,
               y       = Value, 
               ymin    = Value - StdErr, 
               ymax    = Value + StdErr,
               colour  = Result)) +
  geom_pointrange(fatten = 2, size = 1, show.legend = FALSE) +
  scale_color_manual(values = mycols) +
  geom_text(aes(label = pValTex), size = 2.5, col = "#222222",
            nudge_x = -0.4, nudge_y = df$nudge_dir * 0.08) +
  ylab("Estimated performance") + xlab("") +
  facet_grid(Organism ~ Metric, scales = "free") +
  geom_vline(xintercept = 3.35, size = .2, lty = 2) + 
  coord_flip() + 
  theme_light() +
  theme(strip.text  = element_text(colour = "black", face = "bold"),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        plot.margin = margin(2, 2, 2, 0))
dev.off()


# Generate protein prediction map examples for paper
# (version 3)
for (i in seq_along(preds)){
  tmp <- preds[[i]]$myprobs %>%
    dplyr::mutate(Prediction = (1 + RF_OrgSpec_class) / 2,
                  Class      = (1 + Class) / 2,
                  PosPred    = Prediction / Prediction,
                  IsNew      = factor(Prediction & Class,
                                      levels = c(TRUE, FALSE, NA),
                                      labels = c("Known epitope", NA, "New epitope"),
                                      exclude = NULL),
                  Info_UID = paste0(dirs[i], ": ", Info_UID)) %>%
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
    ylim(0, 1.05) + ylab("Prediction") + xlab("Protein Position") +
    theme_light() +
    theme(strip.text  = element_text(colour = "black", face = "bold"),
          axis.text.x = element_text(size = 8),
          plot.margin = margin(3, 3, 3, 3),
          legend.title = element_blank(),
          legend.position = "none",
          legend.background = element_rect(colour = "#aaaaaa"))
  
  for (j in seq(1, ceiling(length(unique(tmp$Info_UID)) / 6))){
    x <- mp + facet_wrap_paginate(Info_UID ~ ., 
                                  scales = "free",
                                  nrow = 2, ncol = 3, page = j)
    ggsave(plot = x,
           filename = paste0("../figures/prot_", dirs[i], "-", j, ".png"),
           width = 10, height = 3, units = "in")
  }
}
