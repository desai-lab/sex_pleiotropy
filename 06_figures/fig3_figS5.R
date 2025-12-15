library(ggplot2)
library(viridis)
library(dplyr)
library(ggpubr)
library(patchwork)
library(grid)
library(gridExtra)

theme_gge <- function (base_size = 10, base_family = "sans")
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title = element_text(size = 7),
      panel.background = element_rect(fill = "white", colour = NA),
      #panel.border = element_rect(fill = NA, colour = "black"),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.25),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      axis.ticks = element_blank(),
      legend.key = element_rect(fill = "white", colour = NA),
      complete = TRUE)
}

clone_dat <- read.csv("clone_dat.csv")

# Create empty lists to store plots/models
combined_plots_list <- list()
lm_results_list <- list()

# List of count variables
count_vars <- c("tot_hits", "multi_count", "single_count")

# Loop over each variable
for (count_var in count_vars) {
  # try first plot:
  fit_mut_asex <- ggplot(clone_dat[clone_dat$regime=='asexual',], aes(x = !!sym(count_var), y = mc_gain, col = env)) +
    facet_wrap( ~ env) + theme_gge() + geom_point() +
    geom_hline(yintercept = 0) +
    scale_color_manual(values = viridis(6))
  
  fit_mut_sex <- ggplot(clone_dat[clone_dat$regime=='sexual',], aes(x = !!sym(count_var), y = mc_gain, col = env)) +
    facet_wrap( ~ env) + theme_gge() + geom_point() +
    geom_hline(yintercept = 0) +
    scale_color_manual(values = viridis(6))
  
  # calculate quantile of mutations for asexuals:
  clone_dat_asex <- dplyr::filter(clone_dat, regime == 'asexual', !!sym(count_var) > 1)
  clone_dat_sex  <- dplyr::filter(clone_dat, regime == 'sexual', !!sym(count_var) > 1)
  
  # write function to flag quantiles
  quantile_flagger <- function(x, qs = NULL){
    
    # use global qs if provided, otherwise compute from x
    if (is.null(qs)) {
      qs <- as.vector(quantile(x[[count_var]]))
    }
    
    cuts <- cut(x[[count_var]], qs, include.lowest = TRUE, right = FALSE)
    levels(cuts) <- c('1','2','3','4')
    x$cuts <- as.numeric(as.vector(cuts))
    
    # z-score fitness within this subset (env)
    x$mc_gain_z <- as.vector(scale(x$mc_gain))
    
    return(x)
  }
  
  # global mutation-count quartiles per regime
  qs_asex <- quantile(clone_dat_asex[[count_var]])
  qs_sex  <- quantile(clone_dat_sex[[count_var]])
  
  # ASEXUAL
  # calculate quantiles:
  clone_dat_asex_l  <- split(clone_dat_asex, clone_dat_asex$env)
  clone_dat_asex_l2 <- lapply(clone_dat_asex_l, quantile_flagger, qs = qs_asex)
  clone_dat_asex_l3 <- do.call(rbind, clone_dat_asex_l2)
  
  # re-level environments:
  clone_dat_asex_l3$env <- factor(clone_dat_asex_l3$env, levels = c('YPD','SC30C','SC37C','SC_pH7.3','SC_0.2M_NaCl','lowP'))
  
  # SEXUAL
  # calculate quantiles:
  clone_dat_sex_l  <- split(clone_dat_sex, clone_dat_sex$env)
  clone_dat_sex_l2 <- lapply(clone_dat_sex_l, quantile_flagger, qs = qs_sex)
  clone_dat_sex_l3 <- do.call(rbind, clone_dat_sex_l2)
  
  # re-level environments:
  clone_dat_sex_l3$env <- factor(clone_dat_sex_l3$env, levels = c('YPD','SC30C','SC37C','SC_pH7.3','SC_0.2M_NaCl','lowP'))
  
  # calculate color scaling for mutation number across clones:
  # scale to viridis and go full range
  mut_hits_sorted <- sort(c(clone_dat_asex_l3[[count_var]],clone_dat_sex_l3[[count_var]]))
  
  custom_labels <- c(
    'YPD' = 'Home (YPD)', 
    'SC30C' = '30°C', 
    'SC37C' = '37°C', 
    'SC_pH7.3' = 'pH 7.3', 
    'SC_0.2M_NaCl' = '+NaCl', 
    'lowP' = '-P'
  )
  quartile_labels <- as_labeller(c('1'="0-25%", '2'="25-50%", '3'="50-75%", '4'="75-100%"))
  
  # ASEXUAL
  fit_mut_z_asex <- ggplot(
    clone_dat_asex_l3,
    aes(x = env, y = mc_gain_z, col = !!sym(count_var), group = clone_id)
  ) +
    facet_wrap(~ cuts, ncol = 4, labeller = quartile_labels) +
    geom_line(alpha = 0.3) +          # 50% opacity lines
    geom_point() +                    # points in same color as line
    theme_gge() +
    geom_hline(yintercept = 0) +
    scale_color_viridis(
      name   = "Mutation Count",
      limits = c(1, ceiling(max(clone_dat_asex_l3[[count_var]])/10) * 10),
      guide  = guide_colorbar(title.position = "top")
    ) +
    labs(title = "Asexual") +
    theme(
      axis.text.x  = element_blank(),
      axis.line.x  = element_blank(),
      axis.title.x = element_blank(),
      legend.position = c(0, 1.25),
      legend.direction = "horizontal",
      legend.title = element_text(size = 10, hjust = 0.5),
      axis.title.y = element_blank(),
      plot.title   = element_text(size = 10, hjust = 0.02, vjust = -10, face = "bold")
    )
  
  # SEXUAL
  fit_mut_z_sex <- ggplot(
    clone_dat_sex_l3,
    aes(x = env, y = mc_gain_z, col = !!sym(count_var), group = clone_id)
  ) +
    facet_wrap(~ cuts, ncol = 4) +
    geom_line(alpha = 0.3) +          # 50% opacity lines
    geom_point() +                    # points in same color as line
    theme_gge() +
    geom_hline(yintercept = 0) +
    scale_color_viridis(
      limits = c(1, ceiling(max(clone_dat_asex_l3[[count_var]])/10) * 10)
    ) +
    scale_x_discrete(labels = custom_labels) +
    labs(x = "Assay Environment", title = "Sexual") +
    theme(
      axis.text.x   = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "none",
      strip.text    = element_blank(),
      axis.title.x  = element_text(size = 12, margin = margin(t = 10)),
      axis.title.y  = element_blank(),
      plot.title    = element_text(size = 10, hjust = 0.02, vjust = -4, face = "bold")
    )
  
  # Create a grob for the horizontal line and label
  line_and_label <- wrap_elements(
    gridExtra::grid.arrange(
      grid::textGrob("Mutation Count Quartile", gp = grid::gpar(fontsize = 10)),
      grid::segmentsGrob(x0 = 0.125, x1 = 0.875, y0 = -9, y1 = -9, gp = grid::gpar(col = "black", lwd = 1)), 
      heights = c(0.9, 0.1)
    )
  )
  
  # Combine the plots vertically
  mut_quartile_plot_main <- line_and_label / (fit_mut_z_asex / fit_mut_z_sex) +
    plot_layout(heights = c(0.02, 1))
  
  y_label <- wrap_elements(full = grid::textGrob("Fitness Gain (z-scored)", rot = 90, gp = gpar(fontsize = 12)))
  mut_quartile_plot <- y_label + mut_quartile_plot_main +
    plot_layout(widths = c(0.05, 1)) +
    plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))
  
  # now: plot average pleiotropic z-score versus mutation number as single biplot:
  pleio_means_asex <- dplyr::filter(clone_dat_asex_l3, env != 'YPD') %>% group_by(clone_id, kk_pop_id, regime) %>%
    summarise(pleio_mean_z = mean(mc_gain_z),
              hits = unique(!!sym(count_var)))
  
  pleio_means_sex <- dplyr::filter(clone_dat_sex_l3, env != 'YPD') %>% group_by(clone_id, kk_pop_id, regime) %>%
    summarise(pleio_mean_z = mean(mc_gain_z),
              hits = unique(!!sym(count_var)))
  
  # plot:
  pleio_reg_asex <- ggplot(pleio_means_asex, aes(x = hits, y = pleio_mean_z)) +
    geom_point(col = 'steelblue') +
    theme_gge() +
    theme(axis.title.y = element_text(size = 12, vjust = 3), plot.title = element_text(size = 10, hjust = 0.05, vjust = -4, face = "bold")) +
    scale_y_continuous(limits = c(-2.25,2.25), breaks = seq(-2,2,1)) +
    stat_smooth(method = 'lm', formula = y ~ x, color = 'darkgrey', alpha = 0.2) +
    labs(title = "Asexual", x = NULL, y = "Mean Pleiotropic Fitness Gain (z-scored)")
  
  pleio_reg_sex <- ggplot(pleio_means_sex, aes(x = hits, y = pleio_mean_z)) +
    geom_point(col = 'darkorange2') +
    theme_gge() +
    theme(axis.title.y = element_text(size = 12, vjust = 3), axis.title.x = element_text(size = 12, margin = margin(t = 5)), plot.title = element_text(size = 10, hjust = 0.05, vjust = -4, face = "bold")) +
    scale_y_continuous(limits = c(-2.25,2.25), breaks = seq(-2,2,1)) +
    stat_smooth(method = 'lm', formula = y ~ x, color = 'darkgrey', alpha = 0.2, linetype = 'dotted') +
    labs(title = "Sexual", x = NULL, y = "Mean Pleiotropic Fitness Gain (z-scored)")
  
  # Combine the plots vertically with patchwork
  fitness_mut_plot <- pleio_reg_asex / pleio_reg_sex +
    plot_layout(guides = "collect", axis_titles = "collect") +
    labs(x = "Mutation Count", y = "Mean Pleiotropic Fitness Gain (z-scored)") +
    plot_annotation(
      theme = theme(plot.margin = margin(5, 5, 5, 5))
    )

  # Wrap each plot as an independent element
  wrapped_left_plot <- wrap_elements(mut_quartile_plot)
  wrapped_right_plot <- wrap_elements(fitness_mut_plot)
  
  # Combine the wrapped plots side by side
  combined_plot <- (wrapped_left_plot | wrapped_right_plot) +
    plot_layout(widths = c(3, 1)) +  
    plot_annotation(tag_levels = 'A', tag_prefix = '', tag_suffix = '') & 
    theme(plot.tag = element_text(face = 'bold', size = 25), plot.tag.position = c(0.01, 0.98), plot.margin = margin(5, 5, 5, 5))

  # Store the plot in the list
  combined_plots_list[[count_var]] <- combined_plot
  
  # Linear regression models
  lm_asex <- lm(pleio_mean_z ~ hits, data = pleio_means_asex)
  lm_sex <- lm(pleio_mean_z ~ hits, data = pleio_means_sex)
  print(paste("Results for", count_var, "- Asexual"))
  print(summary(lm_asex))
  print(paste("Results for", count_var, "- Sexual"))
  print(summary(lm_sex))
  
  # Store the models if needed
  lm_results_list[[count_var]] <- list(
    asexual = lm_asex,
    sexual = lm_sex
  )
}

total_hit_plot <- combined_plots_list[["tot_hits"]]
ggsave(total_hit_plot, filename = "figure3.pdf", width = 12, height = 6, units = "in", dpi = 300)
ggsave(total_hit_plot, filename = "figure3.jpg", width = 12, height = 6, units = "in", dpi = 300)

# Extract single-hit and multi-hit plots
single_hit_plot <- combined_plots_list[["single_count"]]
multi_hit_plot <- combined_plots_list[["multi_count"]]

# Label each section
label_single_hit <- wrap_elements(grid::textGrob("Single-hit Mutations", gp = grid::gpar(fontsize = 14, fontface = "bold")))
label_multi_hit <- wrap_elements(grid::textGrob("Multi-hit Mutations", gp = grid::gpar(fontsize = 14, fontface = "bold")))

# Ensure plot_annotation applies ONLY to the main plots, not labels
single_hit_annotated <- wrap_elements(single_hit_plot + plot_annotation(tag_levels = "A"))
multi_hit_annotated <- wrap_elements(multi_hit_plot + plot_annotation(tag_levels = list(c("C", "D"))))

# Stack each label with its corresponding plot (left + right subplots retain labels)
single_hit_section <- wrap_elements((label_single_hit / single_hit_annotated) + plot_layout(heights = c(0.05, 1)))
multi_hit_section <- wrap_elements((label_multi_hit / multi_hit_annotated) + plot_layout(heights = c(0.05, 1)))

# Combine both sections
figureS4 <- (
  single_hit_section /
    multi_hit_section
) + plot_layout(heights = c(1, 1))

# Save the combined figure
ggsave(figureS4, filename = "figureS5.pdf", width = 12, height = 13, units = "in", dpi = 300)
ggsave(figureS4, filename = "figureS5.jpg", width = 12, height = 13, units = "in", dpi = 300)
