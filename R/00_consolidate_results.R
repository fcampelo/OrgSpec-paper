# Consolidate all results

library(dplyr)
library(ggplot2)
library(ggthemes)
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
for (i in seq_along(preds)){
  tmp <- preds[[i]]$myprobs %>%
    dplyr::mutate(IsNew = as.numeric(RF_OrgSpec_class == 1 & is.na(Class)),
                  IsOld = as.numeric(RF_OrgSpec_class == 1 & (Class == 1)),
                  Info_UID = paste0(dirs[i], ":", Info_UID))
  tmp$IsNew[which(!tmp$IsNew)] <- NA
  tmp$IsOld[which(!tmp$IsOld)] <- NA
  tmp$RF_OrgSpec_class[which(tmp$RF_OrgSpec_class == -1)] <- 0
  tmp$Class[which(tmp$Class == -1)] <- NA
  
  # generate 6 panels/figure, split multiple figs if needed
  for (j in 1:ceiling(length(unique(preds[[i]]$myprobs$Info_UID)) / 3)){
    st <- 1 + (j - 1) * 3
    en <- st + 2
    mp <- tmp %>%
      dplyr::filter(Info_UID %in% unique(Info_UID)[st:en]) %>%
      ggplot(aes(x = Info_center_pos)) +
      geom_point(aes(y = 1.03 * Class), size = .2, pch = 20, alpha = .3) + 
      geom_line(aes(y = RF_OrgSpec_class), lwd = .1) + 
      geom_point(aes(y = IsOld), col = "#628DFC", size = .5, pch = 20) +
      geom_point(aes(y = IsNew), col = "#FC8D62", size = 1, pch = 20) +
      geom_line(aes(y = RF_OrgSpec_prob, colour = RF_OrgSpec_prob), 
                lwd = .1, show.legend = FALSE) +
      scale_colour_viridis_c(direction = -1, begin = .25, end = .75) +
      facet_wrap(Info_UID ~ ., scales = "free", ncol = 3) +
      ylim(0, 1.05) + ylab("Prediction") + xlab("Protein Position") +
      theme_light() + 
      theme(strip.text  = element_text(colour = "black", face = "bold"),
            axis.text.x = element_text(size = 8),
            plot.margin = margin(2, 2, 2, 2),
            legend.position = "bottom")
    
    ggsave(plot = mp, 
           filename = paste0("../figures/prot_", dirs[i], "-", j, ".png"), 
           width = 10, height = 2, units = "in")
  }
}










# # Generate protein prediction map examples for paper
library(ggnewscale)
for (i in seq_along(preds)){
  tmp <- preds[[i]]$myprobs %>%
    dplyr::mutate(IsNew = as.numeric(RF_OrgSpec_class == 1 & is.na(Class)))
  tmp$IsNew[which(!tmp$IsNew)] <- NA
  tmp$RF_OrgSpec_class[which(tmp$RF_OrgSpec_class == -1)] <- 0
  
  mp <- tmp %>%
    ggplot(aes(x = Info_center_pos,
               y = RF_OrgSpec_class)) +
    geom_point(aes(colour = as.factor(RF_OrgSpec_class)),
               size = .5, pch = 20, show.legend = FALSE) +
    scale_colour_manual(values = c("#66C2A5", "#628DFC")) +
    geom_line(lwd = .1, show.legend = FALSE) +
    geom_point(aes(y = 1.02 * IsNew, colour = NULL),
               col = "#FC8D62", size = .5, pch = 20, show.legend = FALSE) +
    new_scale_color() +
    geom_line(aes(y = RF_OrgSpec_prob, colour = RF_OrgSpec_prob),
              lwd = .1, show.legend = FALSE) +
    scale_colour_viridis_c(direction = -1, begin = .25, end = .75) +
    facet_wrap(Info_UID ~ ., scales = "free") +
    ylim(0, 1.05) + ylab("Prediction") + xlab("Protein Position") +
    theme_light() +
    theme(strip.text  = element_text(colour = "black", face = "bold"),
          axis.text.x = element_text(size = 8),
          plot.margin = margin(2, 2, 2, 2),
          legend.position = "bottom") +
    labs(col = "Epitope prob.") +
    ggtitle(paste0("Predicted epitopes: ", dirs[i]))
  
  print(mp)
}
# 
