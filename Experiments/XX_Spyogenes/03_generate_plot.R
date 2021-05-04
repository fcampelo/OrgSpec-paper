# generate plot

library(dplyr)
library(ggplot2)
library(reshape2)

# Load results data
mydata <- readRDS("./output/analysis.rds")

# Method and performance metric labels (ordered, for plotting)
methods <- c("ABCpred", "Bepipred2", "iBCE-EL", "LBtope", "SVMtrip",
             "RF-Heter","RF-Hybrid","RF-OrgSpec")
metrics <- c("ACCURACY", "SENS", "PPV", "NPV", "AUC", "MCC")

preds <- list(1)

# Extract predictions data
preds$mypreds <- mydata$mypreds
preds$myprobs <- mydata$myprobs

# Extract performance data
df <- mydata$myperf_pep %>%
  dplyr::filter(!(Metric %in% c("SPEC", "F1"))) # <-- not needed for the paper

# Extract p-values
pvals <- mydata$Pvals_pep %>% 
  dplyr::select(-SPEC)

# Extract ROC curves
cl.cols <- grep("_class", names(mydata$myres_pep))
pr.cols <- grep("_prob", names(mydata$myres_pep))
for (j in seq_along(cl.cols)){
  myperf <- epitopes::calc_performance(truth = mydata$myres_pep$Class,
                                       pred = mydata$myres_pep[, cl.cols[j]],
                                       prob = mydata$myres_pep[, pr.cols[j]],
                                       ret.as.list = TRUE)
  if (j == 1){
    rocs <- data.frame(Method = gsub("_class", "", 
                                     names(mydata$myres_pep)[cl.cols[j]]),
                       TPR    = myperf$tpr,
                       FPR    = myperf$fpr)
  } else {
    rocs <- rbind(rocs,
                  data.frame(Method = gsub("_class", "", 
                                           names(mydata$myres_pep)[cl.cols[j]]),
                             TPR    = myperf$tpr,
                             FPR    = myperf$fpr))
  }
}

pvals <- pvals %>%
  ungroup() %>%
  dplyr::select(-METHOD1) %>%
  dplyr::rename(Method = METHOD2) %>%
  reshape2::melt(id.vars = c("Method", "NoLeak"),
                 variable.name = "Metric",
                 value.name = "pValue") %>%
  dplyr::mutate(Method = gsub("_", "-", Method, fixed = TRUE),
                pValue = round(pValue, 3))

df <- df %>%
  dplyr::mutate(Method = as.character(Method),
                Metric = as.character(Metric)) %>%
  dplyr::left_join(pvals, by = c("Method", "NoLeak", "Metric")) %>%
  dplyr::mutate(Method = factor(Method, levels = methods, ordered = TRUE),
                Metric = factor(Metric, levels = metrics, ordered = TRUE),
                pValue = pValue) %>%
  group_by(Metric) %>%
  mutate(Result = ifelse(Mean >= first(Mean), "better", "worse"))

df$Result[df$pValue > 0.05] <- "non-signif"
df$Result[df$Method == "RF-OrgSpec"] <- "ref. method"
df$Result <- factor(df$Result,
                    levels = c("ref. method", "non-signif", "better", "worse"),
                    ordered = TRUE)
df$pValue[df$pValue <= .01] <- paste0("<0.01")
df$pValue[df$pValue >= 0.9] <- paste0(">0.9")


# Generate the plot 
mp <- ggplot(filter(df, NoLeak == FALSE),
             aes(x = Method, y = Value,
                 ymin = Value - StdErr, ymax = Value + StdErr,
                 colour  = Result)) +
  geom_text(aes(label = pValue), size = 2, col = "#222222", nudge_x = -0.4) +
  geom_pointrange(fatten = 1.5, size = .75, show.legend = FALSE) +
  scale_color_manual(values = c("#555555", "#7570c3", "#1b9e77", "#d95f02")) +
  ylab("Estimated performance") + xlab("") +
  facet_wrap(. ~ Metric, scales = "free", nrow = 2) +
  geom_vline(xintercept = 5.35, size = .2, lty = 2) +
  scale_y_continuous(expand = c(.055,.05), breaks = seq(-1, 1, by = 0.1)) +
  scale_x_discrete(expand = c(.1,.05)) +
  coord_flip() + theme_light() +
  theme(strip.text  = element_text(colour = "black", face = "bold"),
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7))

ggsave(plot = mp, filename = "../../figures/res_Spyogenes.png",
       width = 6.5, height = 4, units = "in")
ggsave(plot = mp, filename = "../../figures/res_Spyogenes.tiff",
       width = 6.5, height = 4, units = "in")


mydata$Pvals_pep.raw %>% 
  dplyr::filter(NoLeak == FALSE) %>%
  dplyr::select(-SPEC, -METHOD1, -NoLeak) %>%
  dplyr::mutate(across(where(is.numeric), ~round(.x, digits = 3))) %>%
  kableExtra::kable(format = "latex")

