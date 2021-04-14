make_plot3 <- function(df, linepos = 3.5, plot_rows = 2,
                       mycols = c("#555555", "#1b9e77", "#d95f02", "#7570b3")) {
  require(ggplot2)

  df$vnudge = -0.2 + 0.4*df$NoLeak
  mp <- ggplot(df, aes(x = Method,
                       y = Mean, ymin = Mean - StdErr, ymax = Mean + StdErr,
                       colour = Type)) +
    scale_color_manual(values = mycols) +
    geom_pointrange(aes(alpha = NoLeak),
                    fatten = 2, size = 1,
                    position = position_nudge(x = df$vnudge)) +
    scale_alpha_discrete(range = c(0.35, 0.9)) +
    guides(alpha = guide_legend("Zero-leakage"), colour = FALSE) +
    ylab("Performance estimates and standard errors") + xlab("") +
    coord_flip() +
    theme_light() +
    theme(strip.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = c(.85, .3)) +
    facet_wrap(. ~ Metric, nrow = plot_rows, scales = "free_x")

  return(mp)
}
