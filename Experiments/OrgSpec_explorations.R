library(ranger)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(Rtsne)

# Format to export figures
fig.format = ".png"

# ========================= Check Feature Importance ========================= #
# Load and prepare data
dirs <- c("./01_EBV/", "./02_HepC/", 
          "./03_Ovolvulus/")
orgs <- sapply(dirs, function(x) strsplit(x, split = "/|_")[[1]][3], 
               USE.NAMES = FALSE)

feat_imp <- vector("list", length(dirs))

for (i in seq_along(dirs)){
  names(feat_imp)[i]    <- orgs[i]
  
  mod <- readRDS(paste0(dirs[i], "/output/TrOrgSpec/model.rds"))
  feat_imp[[i]]$OrgSpec <- ranger::importance(mod)
  
  mod <- readRDS(paste0(dirs[i], "/output/TrHeter/model.rds"))
  feat_imp[[i]]$Heter <- importance(mod)
  
  # Ignore hybrid models, let's focus on main differences of OrgSpec vs. Heterogeneous  
  # models[[i]]$Hybrid   <- readRDS(paste0(dirs[i], "/output/TrHybrid/model.rds"))
  # feat_imp[[i]]$Hybrid <- importance(models[[i]]$Hybrid)
}
rm(mod)

# Consolidate feature importance dataframe
df <- tibble(Organism = character(), Model = character(), 
             Feature = character(), Type = character(), 
             Importance = numeric())

Feat_group <- tibble(Feat = names(feat_imp[[i]]$Heter),
                     Type = rep("AAdescr.", length(feat_imp[[i]]$Heter)))
Feat_group$Type[grep("feat_seq_entropy", Feat_group$Feat)] <- "Entropy"
Feat_group$Type[grep("feat_molecular_weight", Feat_group$Feat)] <- "Mol. weight"
Feat_group$Type[grep("atoms$", Feat_group$Feat)] <- "Atoms"
Feat_group$Type[grep("feat_Perc_[a-zA-Z]+$", Feat_group$Feat)] <- "Freq-Types"
Feat_group$Type[grep("feat_Perc_[A-Z]$", Feat_group$Feat)] <- "Freq-1AA"
Feat_group$Type[grep("feat_Perc_[A-Z][A-Z]$", Feat_group$Feat)] <- "Freq-2AA"
Feat_group$Type[grep("feat_CT", Feat_group$Feat)] <- "CT"

for (i in seq_along(feat_imp)){
  for (j in seq_along(feat_imp[[i]])){
    tmp <- tibble(Feature    = names(feat_imp[[i]][[j]]),
                  Importance = as.numeric(feat_imp[[i]][[j]]) / sum(feat_imp[[i]][[j]]),
                  Organism = names(feat_imp)[i],
                  Model    = names(feat_imp[[i]])[j]) %>%
      left_join(Feat_group, by = c("Feature" = "Feat"))
    df <- df %>%
      bind_rows(tmp)
  }
}
rm(tmp)

df <- df %>%
  mutate(Feature = gsub("feat_", "", Feature))

# Compare distribution of normalised importance values 
df %>%
  group_by(Organism, Model) %>%
  arrange(desc(Importance), .by_group = TRUE) %>%
  mutate(Rank = 1:n(),
         Type = factor(Type)) %>%
  ungroup() %>%
  ggplot(aes(x = Rank, y = Importance, colour = Type)) + 
  geom_point(pch = 16, alpha = .6, stroke = 0) + 
  #scale_y_log10() +
  facet_grid(Organism ~ Model) + 
  theme_light() +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(colour = "black"),
        legend.position = "bottom") + 
  scale_color_brewer(type = "qual", drop = FALSE) + 
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) + 
  ggtitle("Normalised feature importance")

ggsave(last_plot(), width = 9, height = 7,
       filename = paste0("../figures/01 Normalised feature importance", 
                         fig.format))

# Check feature importance in heterogeneous vs. orgspec models
df %>% 
  pivot_wider(id_cols = c("Organism", "Feature", "Type"), 
              names_from = "Model",
              values_from = "Importance") %>%
  mutate(Type = factor(Type)) %>%
  ggplot(aes(x = OrgSpec, y = Heter, colour = Type)) + 
  geom_point(pch = 16, alpha = .6, stroke = 0) + 
  geom_smooth(aes(colour = NULL), col = 1, lty = 2, method = "lm", se = FALSE)+
  scale_x_log10() + scale_y_log10() +
  facet_wrap(Organism ~ ., nrow = 1) + 
  theme_light() +
  theme(strip.text = element_text(colour = "black"),
        legend.position = "bottom") + 
  scale_color_brewer(type = "qual", drop = FALSE) + 
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) + 
  ggtitle("Normalised feature importance: OrgSpec vs. Heterogeneous models") + 
  xlab("Importance: OrgSpec") + ylab("Importance: Heterogeneous")

ggsave(last_plot(), width = 9, height = 5,
       filename = paste0("../figures/02 Feature importance OrgSpec vs Heter", 
                         fig.format))

df2 <- df %>% 
  pivot_wider(id_cols = c("Model", "Feature", "Type"), 
              names_from = "Organism",
              values_from = "Importance") 

# Compare feature importance in different OrgSpec models

ggplot(df2, aes(x = EBV, y = Ovolvulus, colour = Type)) + 
  geom_point(pch = 16, alpha = .6, stroke = 0, show.legend = FALSE) + 
  geom_smooth(aes(colour = NULL), col = 1, lty = 2, method = "lm", se = FALSE) +
  scale_x_log10() + scale_y_log10() + 
  theme_light() +
  scale_color_brewer(type = "qual", drop = FALSE) + 
  guides(colour = guide_legend(override.aes = list(size=5))) + 
  #ggtitle("Normalised feature importance: EBV vs. O. volvulus") + 
  xlab("Importance: EBV OrgSpec") + ylab("Importance: O. volvulus OrgSpec")
ggsave(last_plot(), width = 8, height = 3,
       filename = paste0("../figures/03a Feature importance EBV vs Ovolvulus OrgSpec", 
                         fig.format))

ggplot(df2, aes(x = EBV, y = HepC, colour = Type)) + 
  geom_point(pch = 16, alpha = .6, stroke = 0, show.legend = FALSE) + 
  geom_smooth(aes(colour = NULL), col = 1, lty = 2, method = "lm", se = FALSE) +
  scale_x_log10() + scale_y_log10() + 
  theme_light() +
  scale_color_brewer(type = "qual", drop = FALSE) + 
  guides(colour = guide_legend(override.aes = list(size=5))) + 
  #ggtitle("Normalised feature importance: EBV vs. Hep C") + 
  xlab("Importance: EBV OrgSpec") + ylab("Importance: HepC OrgSpec")
ggsave(last_plot(), width = 8, height = 3,
       filename = paste0("../figures/03b Feature importance EBV vs HepC OrgSpec", 
                         fig.format))

ggplot(df2, aes(x = Ovolvulus, y = HepC, colour = Type)) + 
  geom_point(pch = 16, alpha = .6, stroke = 0) + 
  geom_smooth(aes(colour = NULL), col = 1, lty = 2, method = "lm", se = FALSE) +
  scale_x_log10() + scale_y_log10() + 
  theme_light() +
  scale_color_brewer(type = "qual", drop = FALSE) + 
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) + 
  #ggtitle("Normalised feature importance: O. volvulus vs. Hep C") + 
  xlab("Importance: O. volvulus OrgSpec") + ylab("Importance: HepC OrgSpec") + 
  theme(legend.position = c(.1, .5), 
        legend.background = element_blank())
ggsave(last_plot(), width = 8, height = 3,
       filename = paste0("../figures/03c Feature importance Ovolvulus vs HepC OrgSpec", 
                         fig.format))


# Extract top K features for each Organism and Model type
K = 30
dfK <- df %>%
  mutate(Type = factor(Type)) %>%
  group_by(Organism, Model) %>%
  arrange(desc(Importance), .by_group = TRUE) %>%
  mutate(Rank = 1:n()) %>%
  ungroup()

df3 <- dfK %>%
  filter(Rank <= K) %>%
  select(-Importance)

df3 %>%
  pivot_wider(id_cols = c("Rank", "Organism"), 
              names_from = "Model",
              values_from = "Feature") %>%
  print(n = Inf)

# Plot top K features
ggplot(df3, aes(x = Model, y = Rank, label = Feature, 
                group = Feature)) + 
  geom_point(aes(colour = Type), size = .2) +
  geom_line(aes(colour = Type), alpha = .4, show.legend = FALSE) +
  geom_label(aes(fill = Type), alpha = .6, size = 3, show.legend = FALSE) +
  scale_color_brewer(type = "qual", drop = FALSE) +
  scale_fill_brewer(type = "qual", drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_y_reverse(breaks = 1:K) +
  facet_grid(. ~ Organism) + 
  theme_light() + 
  theme(panel.grid = element_blank(),
        strip.text = element_text(colour = "black"),
        legend.position = "top")
ggsave(last_plot(), width = 9, height = 12,
       filename = paste0("../figures/04a Top Features by Org x Model Type", 
                         fig.format))

# Another way to look at the same data
ggplot(df3, aes(x = Organism, y = Rank, label = Feature, 
                group = Feature)) + 
  geom_point(aes(colour = Type), size = .2) +
  geom_line(aes(colour = Type), alpha = .4, show.legend = FALSE) +
  geom_label(aes(fill = Type), alpha = .6, size = 3, show.legend = FALSE) +
  scale_color_brewer(type = "qual", drop = FALSE) +
  scale_fill_brewer(type = "qual", drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_y_reverse(breaks = 1:K) +
  facet_grid(. ~ Model) + 
  theme_light() + 
  theme(panel.grid = element_blank(),
        strip.text = element_text(colour = "black"),
        legend.position = "top")
ggsave(last_plot(), width = 9, height = 12,
       filename = paste0("../figures/04b Top Features by Model Type x Org", 
                         fig.format))

# Occurrence among top K for OrgSpec and heterogeneous models
for (k in 1:30){
  tmp <- dfK %>%
    filter(Rank <= k) %>%
    select(-Importance) %>%
    group_by(Model, Feature) %>%
    summarise(Count = n(),
              Type  = first(Type), 
              topK  = k)
  if (k == 1){
    df.ani <- tmp
  } else {
    df.ani <- df.ani %>% bind_rows(tmp)
  }
}

df.ani %>%
  filter(topK %in% c(5, 10, 20, 30)) %>%
  pivot_wider(id_cols = c("Feature", "Type", "topK"), 
              names_from = Model, values_from = Count) %>%
  mutate(Heter = ifelse(is.na(Heter), 0, Heter),
         OrgSpec = ifelse(is.na(OrgSpec), 0, OrgSpec),
         topK    = factor(sprintf("Top %02d features", topK),
                          ordered = TRUE)) %>%
  ggplot(aes(x = Heter, y = OrgSpec)) + 
  geom_point(aes(colour = Type), size = 4) + 
  xlim(-.5, 3.5) + ylim(-.5, 3.5) +
  facet_wrap(topK ~ ., ncol = 2) +
  geom_text_repel(aes(label = Feature, colour = Type), 
                  force = 2,
                  max.overlaps = 50, show.legend = FALSE) +
  scale_color_brewer(type = "qual", drop = FALSE) +
  scale_fill_brewer(type = "qual", drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_light() + 
  theme(strip.text = element_text(colour = "black"),
        legend.position = "top") +
  xlab(paste0("Feature Occurrences among top K: Heter models")) + 
  ylab(paste0("Feature Occurrences among top K: OrgSpec models"))
ggsave(last_plot(), width = 15, height = 12,
       filename = paste0("../figures/04c - feature occurrence by model", 
                         fig.format))




# ========================= Check Neighbourhoods ========================= #

recalculate <- FALSE

if (!recalculate){
  proj.data <- readRDS("./proj_data.rds")
  tsne.df   <- readRDS("./proj_data_tsne.rds")
  tsne.all  <- readRDS("./proj_data_tsne_all.rds")
} else {
  for (i in seq_along(dirs)){
    tmp <- readRDS(paste0(dirs[i], "/data/splits/01_training.rds")) %>%
      as_tibble() %>%
      select(-starts_with("Info")) %>%
      mutate(Organism = orgs[i]) %>%
      select(Organism, starts_with("feat"), Class)
    
    if (i == 1) {
      alldata <- tmp
    } else {
      alldata <- bind_rows(alldata, tmp)
    }
  }
  rm(tmp)
  
  # Add a heterogeneous dataset
  alldata <- bind_rows(
    alldata, 
    readRDS("./03_Ovolvulus/data/heterogeneous_data/df_heterogeneous.rds") %>%
      as_tibble() %>%
      select(-starts_with("Info")) %>%
      mutate(Organism = "Heterogeneous") %>%
      select(Organism, starts_with("feat"), Class))
  
  # Remove CTs and dipeptide features
  alldata <- alldata[, c(1:83, 427:446, 847)]
  
  proj.data <- alldata %>%
    mutate(across(starts_with("feat_"), ~(.x - min(.x))/(1e-12 + diff(range(.x)))))
  
  idx <- duplicated(proj.data %>% select(starts_with("feat_")))
  
  proj.data <- proj.data[which(!idx), ]
  
  rm(alldata)
  saveRDS(proj.data, "./proj_data.rds")
  
  for (i in seq_along(unique(proj.data$Organism))){
    for (perp in c(5,25,100)){
      cat("\n", orgs[i], ": Perp = ", perp, "\n")
      tsne <- proj.data %>%
        filter(Organism == orgs[i]) %>%
        select(starts_with("feat_")) %>%
        Rtsne(dims        = 2, 
              theta       = 0.1,
              perplexity  = perp, 
              verbose     = TRUE,
              max_iter    = 2000,
              num_threads = 7)
      tmp <- as_tibble(tsne$Y) %>%
        mutate(Organism = orgs[i],
               Class    = as.factor(filter(proj.data, 
                                           Organism == orgs[i])$Class),
               Perplexity = perp)
      if (i == 1 & perp == 5) {
        tsne.df <- tmp
      } else {
        tsne.df <- rbind(tsne.df, tmp)
      }
    }
  }
  
  saveRDS(tsne.df, "./proj_data_tsne.rds")
  
  # unified tSNE projection 
  tsne.all <- proj.data %>%
    select(starts_with("feat_")) %>%
    Rtsne(dims        = 2, 
          theta       = 0.1,
          perplexity  = 100, 
          verbose     = TRUE,
          max_iter    = 2000,
          num_threads = 6)
  
  tsne.all <- as_tibble(tsne.all$Y) %>%
    mutate(Organism = proj.data$Organism,
           Class    = as.factor(proj.data$Class))
  
  saveRDS(tsne.all, "./proj_data_tsne_all.rds")
  
}


tsne.df %>% 
  mutate(Perplexity = factor(paste0("Perpl. = ", Perplexity), 
                             levels = paste0("Perpl. = ", 
                                             unique(tsne.df$Perplexity))),
         Class = ifelse(Class == 1, "Positive", "Negative")) %>%
  ggplot(aes(x = V1, y = V2, colour = Class)) + 
  geom_point(pch = 15, alpha = .1, size = 2, stroke = 0) + 
  facet_wrap(Organism ~ Perplexity, scales = "free") + 
  scale_color_brewer(type = "qual") +
  guides(colour = guide_legend(override.aes = list(alpha=1)))+
  theme_light() + 
  theme(strip.text = element_text(colour = "black"))
ggsave(last_plot(), width = 12, height = 12,
       filename = paste0("../figures/05 tSNE by Org with col=Class", 
                         fig.format))

tsne.all %>%
  mutate(Organism = relevel(factor(Organism), ref = "Heterogeneous"),
         Class = ifelse(Class == 1, "Positive", "Negative")) %>%
  ggplot(aes(x = V1, y = V2, 
             colour = Class)) + 
  geom_point(pch = 15, alpha = .1, size = 2, stroke = 0) +
  scale_color_brewer(type = "qual") +
  guides(colour = guide_legend(override.aes = list(alpha=1)))+
  theme_light() + 
  theme(strip.text = element_text(colour = "black")) + 
  facet_wrap(Organism ~ ., nrow = 2) + 
  ggtitle("t-SNE projection", 
          subtitle = "(Heterogeneous contains examples from many distinct pathogens, including samples from EBV, HepC and Ovolvulus)")
ggsave(last_plot(), width = 9, height = 9,
       filename = paste0("../figures/06 tSNE unified projection", 
                         fig.format))

tsne.df %>% 
  filter(Perplexity == 100) %>%
  mutate(Class = ifelse(Class == 1, "Positive", "Negative")) %>%
  ggplot(aes(x = V1, y = V2)) + 
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", contour = FALSE) +
  geom_point(pch = 16, stroke = 0, size = .2) +
  facet_grid(Organism ~ Class) + 
  scale_fill_distiller(palette=4, direction=1) +
  theme_light() + 
  theme(strip.text = element_text(colour = "black"))
ggsave(last_plot(), width = 9, height = 10,
       filename = paste0("../figures/07 tSNE densities", 
                         fig.format))
