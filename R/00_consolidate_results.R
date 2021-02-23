# Consolidate all results

library(dplyr)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(tikzDevice)

dirs <- c("EBV", "HepC", "Ovolvulus", "Pvivax")
preds <- vector("list", length(dirs))
names(preds) <- dirs

for (i in seq_along(dirs)){
  mydata     <- readRDS(paste0("../Experiments/", 
                               dirs[i], "/output/analysis.rds")) 
  preds[[i]] <- readRDS(paste0("../Experiments/", dirs[i], 
                               "/output/consolidated_preds.rds"))
  
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

ggsave(plot = mp, filename = "../figures/res_all_pathogens.png", 
       width = 10, height = 10, units = "in")

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
  mp <- ggplot(preds[[i]]$myprobs,
               aes(x = Info_center_pos,
                   y = RF_OrgSpec_class,
                   colour = RF_OrgSpec_prob)) +
    geom_point(size = .2, pch = 20) + 
    geom_line(lwd = .1, show.legend = FALSE) +
    geom_line(aes(y = RF_OrgSpec_prob), 
              lwd = .1, show.legend = FALSE) +
    scale_colour_viridis_c(direction = -1) +
    facet_wrap(Info_UID ~ ., scales = "free") +
    geom_hline(yintercept = 0.5, lty = 3, lwd = .2, show.legend = FALSE) +
    ylim(-1, 1.2) + 
    theme(strip.text  = element_text(colour = "black", face = "bold"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          plot.margin = margin(2, 2, 2, 0))
  print(mp)
}


# To plot protein results (not ready yet)
# 
# plt3 <- ggplot(myres, aes(x = Info_center_pos, y = RF_os)) +
# geom_point(aes(y = Class), colour = "#FF0000", alpha = 0.5, pch = 20) +
#   geom_line() +
#   facet_wrap(Info_UID ~ ., scales = "free", ncol = 3) +
#   scale_y_discrete() +
#   xlab("Protein position") + ylab("") +
#   theme_light() +
#   theme(strip.text = element_text(colour = "black"))
# 
# ggplot(dplyr::filter(myres, Info_UID == "A0A044V9S3"), 
#        aes(x = Info_center_pos, y = RF_OrgSpec_class)) + 
#   geom_line(lty = 2, alpha = .5) + 
#   geom_point(aes(colour = as.factor(RF_OrgSpec_class)), show.legend = FALSE) + 
#   ylab("Prediction") + xlab("Position") + 
#   scale_y_continuous(breaks = c(-1,1), labels = c("Neg.", "Pos.")) + 
#   theme_light()
# 
# ggplot(myprobs, 
#        aes(x = Info_center_pos, 
#            y = RF_OrgSpec_class, 
#            colour = RF_OrgSpec_prob)) + 
#   geom_point(size = .2) + geom_line(lwd = .1) +
#   facet_wrap(Info_UID ~ ., scales = "free", ncol = 3) + 
#   ylim(-1, 1) + scale_colour_continuous(type = "viridis")

