make_plot3 <- function(df, linepos = 3.5, plot_rows = 2,
                       mycols = c("#555555", "#1b9e77", "#d95f02", "#7570b3")) {
  require(ggplot2)
  require(ggthemes)
  
  # Offset MCC to facilitate faceting
  df[grep("MCC", df$Metric), c("Value", "Mean")] <-
    0.5 * (1 + df[grep("MCC", df$Metric), c("Value", "Mean")])
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
    annotate("segment", y = 0, yend = 1, x = linepos, xend = linepos,
             col = "#bbbbbb") +
    facet_wrap_custom(. ~ Metric, nrow = plot_rows, scales = "free_x",
                      scale_overrides = list(
                        scale_override(6,
                                       scale_y_continuous(breaks = c(0, .25, .5, .75, 1),
                                                          labels = c("-1.00", "-0.50", "0.00", "0.50", "1.00"))))) +
    ylim(0, 1)
  
  return(mp)
}



# Adapted from https://fishandwhistle.net/post/2018/modifying-facet-scales-in-ggplot2/

scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}


CustomFacetWrap <- ggplot2::ggproto(
  "CustomFacetWrap", ggplot2::FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggplot2::ggproto_parent(ggplot2::FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)


facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) ||
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}



CustomFacetGrid <- ggplot2::ggproto(
  "CustomFacetGrid", ggplot2::FacetGrid,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggplot2::ggproto_parent(ggplot2::FacetGrid, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)


facet_grid_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_grid
  facet_super <- ggplot2::facet_grid(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) ||
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetGrid,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}
