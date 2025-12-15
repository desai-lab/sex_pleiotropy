library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)
library(Cairo)

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
set.seed(24)

# Define kk_pop_ids to remove
excluded_kk_pop_ids <- c("a1.1_H06", "a3_H02", "a3_A03")

# load aggregated clone data:
CD4 <- read.csv("parsed_clone_fitness_data.csv", T)
CD4$anc <- 0
CD4$anc[grep('MJM',paste0(CD4$kk_pop_id))] <- 1

# remove outliers that are impossible to accurately estimate:
CD4 <- CD4[abs(CD4$fitness_gain_avg)<0.2,]

CD4 <- CD4 %>%
  filter(!kk_pop_id %in% excluded_kk_pop_ids)

# mean center the populations around the total env-specific means:
env_means_clones <- dplyr::group_by(CD4, env) %>% summarise(env_mean = mean(fitness_gain_avg))
CD4b <- merge(CD4, env_means_clones, by = 'env')
ranges2 <- dplyr::group_by(CD4b[CD4b$anc==0,], env, regime) %>%
  summarise(mn = min(fitness_gain_avg - env_mean), mx = max(fitness_gain_avg - env_mean))
CD4b$mc_avg <- CD4b$fitness_gain_avg - CD4b$env_mean # calculate mean-centered fitness gain

CD4b <- CD4b %>%
  mutate(env = factor(env,
                      levels = c('YPD', 'SC30C', 'SC37C', 'SC_pH7.3', 'SC_0.2M_NaCl', 'lowP'),
                      labels = c('YPD', '30°C', '37°C', 'pH 7.3', '+NaCl', '-P')))

# define ranges for plotting lines:
ranges <- dplyr::group_by(CD4b[CD4b$fitness_gain_avg<0.2,], env, regime) %>% summarise(mn = min(fitness_gain_avg),mx = max(fitness_gain_avg))

clone_pointplot <- ggplot(CD4b %>% filter(!is.na(regime), fitness_gain_avg < 0.2), 
                        aes(x = regime, y = fitness_gain_avg, color = regime)) +
  # Jittered points with slightly smaller size
  geom_jitter(width = 0.2, size = 1.8, alpha = 0.6) +
  # Horizontal line for the median, wider across the x-axis
  stat_summary(fun = median, geom = "crossbar", width = 0.8, size = 0.6) +  # Wider median line
  # Vertical line showing the range (min to max)
  geom_errorbar(stat = "summary", fun.min = min, fun.max = max, width = 0.1, size = 0.5) +
  annotate("rect", ymin = -Inf, ymax = 0, xmin = -Inf, xmax = Inf, alpha = 0.1, fill = "gray50") +
  geom_hline(yintercept = 0, lty = 'dashed',col='gray40',alpha=0.5) +
  facet_wrap(~ env, scales = 'fixed', nrow = 1, strip.position = "bottom") +  # Side-by-side subplots
  scale_y_continuous(limits = c(-0.15,0.15), breaks = seq(-0.2,0.15,0.05)) +
  theme_gge() +
  theme(
    legend.position = c(0.90, 0.95),
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    strip.background = element_blank(),  # Remove background for facet labels
    strip.text = element_text(size = 14, margin = margin(t = 30, b = 5)),  # Make facet labels bigger
    panel.spacing = unit(1, "lines")  # Add spacing between panels if needed
  ) +
  scale_color_manual(values = c('steelblue', 'darkorange2')) +  # Colors for jittered points and median lines
  guides(color = guide_legend(override.aes = list(shape = 'circle', linetype = 0, size = 4, alpha = 1))) +  # Use lines in the legend
  ylab("Fitness Gain") + xlab("")

# fitness_plot <- clone_pointplot +
#   labs(tag = "a") +
#   plot_annotation(
#     theme = theme(plot.margin = margin(5, 5, 5, 5))
#   ) & 
#   theme(plot.tag = element_text(face = 'bold', size = 25), plot.tag.position = c(0.01, 0.98), plot.margin = margin(5, 5, 5, 5))
# 
# fitness_plot

# ggsave(fitness_plot, filename = 'figureS3a.pdf', width = 8, height = 6, units = "in", dpi = 300)
# ggsave(fitness_plot, filename = 'figureS3a.jpg', width = 8, height = 6, units = "in", dpi = 300)

clone_pcplot <- ggplot(CD4b, aes(x = env, y = fitness_gain_avg, col = regime, group = sample_id)) + geom_line(alpha= 0.5, size = 0.75) + theme_gge() +
  scale_color_manual(values = c('steelblue','darkorange2')) +
  scale_y_continuous(limits = c(-0.15,0.15), breaks = seq(-0.2,0.15,0.05)) +
  scale_x_discrete(expand = c(0.075, 0)) +  # Reduce white space on either side
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.position = 'none'
  ) +
  geom_hline(yintercept = 0, lty = 'dashed',col='gray40',alpha=0.5) +
  xlab('') + ylab('Fitness Gain') #+
#coord_flip()
clone_pcplot

# ggsave(clone_pcplot, filename = 'figureS3b.pdf', width = 8, height = 5, units = 'in', dpi = 300)
# ggsave(clone_pcplot, filename = 'figureS3b.jpg', width = 8, height = 5, units = 'in', dpi = 300)

# Combine the two plots vertically
fitness_plot <- clone_pointplot / clone_pcplot + 
  plot_layout(heights = c(1, 1))+  # Adjust relative heights of the plots
  plot_annotation(tag_levels = 'A', tag_prefix = '', tag_suffix = '') & # Add "a)" and "b)" to the plots
  theme(
    plot.tag = element_text(face = 'bold', size = 25),
    plot.tag.position = c(0.01, 1.02)
  )

fitness_plot

# ggsave(fitness_plot, filename = 'figureS3ab.pdf', width = 8, height = 10, units = 'in', dpi = 300)
# ggsave(fitness_plot, filename = 'figureS3ab.jpg', width = 8, height = 10, units = 'in', dpi = 300)

#### VARIANCE CALCULATIONS AND PLOTS ####
# calculate summary statistics for plotting:
#c_var_e <- dplyr::group_by(CD4b[CD4b$anc==0,], regime, kk_pop_id) %>% summarise(var_e = sqrt(var(fitness_gain_avg)))
c_var_e <- dplyr::group_by(CD4b[CD4b$anc==0,], regime, kk_pop_id, env) %>%
  summarise(fitness_gain_avg_pop = mean(fitness_gain_avg)) %>%
  group_by(regime,kk_pop_id) %>%
  summarise(var_e = sqrt(var(fitness_gain_avg_pop)))
c_var_e$regime <- factor(c_var_e$regime, levels = c("asexual", "sexual"))

var_e_plot1 <- ggplot(c_var_e, aes(x = regime, y = var_e, color = as.character(regime))) + 
  geom_jitter(aes(fill = regime), width = 0.1, shape = 21, size = 3, color = "black", stroke = 1) +
  scale_fill_manual(values = c("sexual" = "darkorange2", "asexual" = "steelblue")) +
  stat_summary(fun = median, geom = "crossbar", width = 0.7, size = 0.8, color = "black") +
  labs(y = expression(SD[fitness~gain]), x = NULL, title = "Across-Environment") +
  theme_gge() +
  scale_x_discrete(labels = c("asexual" = "Asexual", "sexual" = "Sexual")) +
  theme(
    legend.position = 'none', 
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 14), 
    axis.text.x = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 14, hjust = 0.5, vjust = 3)
  )

var_e_plot1

#ggsave(var_e_plot1, filename = 'figure2c.pdf', width = 2.5, height = 5, units = "in", dpi = 300)

# average within pops; calculate t-test:
with(c_var_e, t.test(var_e ~ regime, var.equal = FALSE))

## Making each point a clone instead of a population (average of 2 clones)
c_var_e2 <- dplyr::group_by(CD4b[CD4b$anc==0,], regime, sample_id) %>%
  summarise(var_e = sqrt(var(fitness_gain_avg)))
c_var_e2$regime <- factor(c_var_e2$regime, levels = c("asexual", "sexual"))

var_e_plot2 <- ggplot(c_var_e2, aes(x = regime, y = var_e, color = as.character(regime))) + 
  geom_jitter(aes(fill = regime), width = 0.1, shape = 21, size = 3, color = "black", stroke = 1) +
  scale_fill_manual(values = c("sexual" = "darkorange2", "asexual" = "steelblue")) +
  stat_summary(fun = median, geom = "crossbar", width = 0.7, size = 0.8, color = "black") +
  labs(y = expression(SD[fitness~gain]), x = NULL, title = "Across-Environment") +
  theme_gge() +
  scale_x_discrete(labels = c("asexual" = "Asexual", "sexual" = "Sexual")) +
  theme(
    legend.position = 'none', 
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 14), 
    axis.text.x = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 14, hjust = 0.5, vjust = 3)
  )

var_e_plot2

# ggsave(var_e_plot2, filename = 'across_env_clones.pdf', width = 2.5, height = 5, units = "in", dpi = 300)

# average within pops; calculate t-test:
with(c_var_e2, t.test(var_e ~ regime, var.equal = FALSE))


c_var_g <- dplyr::group_by(CD4b[CD4b$anc==0,], regime, env) %>% summarise(var_g = sqrt(var(fitness_gain_avg)))
c_var_g$env <- factor(c_var_g$env, levels = c('YPD', '30°C', '37°C', 'pH 7.3', '+NaCl', '-P'))
c_var_g$regime <- factor(c_var_g$regime, levels = c('asexual','sexual'))
clone_var_g <- ggplot(c_var_g, aes(x = regime, y = var_g, group = env, col = env)) +
  geom_line(color = "black") +
  geom_point(aes(fill = regime), shape = 21, size = 3, color = "black", stroke = 1) + 
  labs(y = expression(SD[fitness~gain]), x = NULL, title = "Across-Clone") +
  theme_gge() +
  theme(
    legend.position = 'none', 
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 14), 
    axis.text.x = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 14, hjust = 0.5, vjust = 3)
  ) +
  scale_fill_manual(values = c("sexual" = "darkorange2", "asexual" = "steelblue")) +
  scale_y_continuous(limits = c(0,0.04), breaks = seq(0,0.04,0.01)) +
  scale_x_discrete(labels = c("asexual" = "Asexual", "sexual" = "Sexual")) +
  geom_text(
    data = c_var_g %>% filter(regime == "sexual"),  # Filter only sexual points
    aes(label = env, y = var_g + ifelse(env == "pH 7.3", 0, ifelse(env == "37°C", 0, ifelse(env == "-P", 0, 0)))),  # Use env as the label
    position = position_nudge(x = 0.15),
    hjust = 0,
    size = 3,  # Text size
    color = "black"  # Text color
  )

#ggsave(clone_var_g, filename = 'figure2b.pdf', width = 2.5, height = 5, units = "in", dpi = 300)

std_dev_plot <- clone_var_g + plot_spacer() + var_e_plot2 +
  plot_layout(ncol = 1, heights = c(1, 0.1, 1)) +
  plot_annotation(tag_levels = list(c("C", "D")), tag_prefix = '', tag_suffix = '') &
  theme(
    plot.tag = element_text(face = 'bold', size = 25),
    plot.tag.position = c(0.1, 1.01),
    plot.margin = margin(5, 5, 5, 5),
  )

std_dev_plot

# ggsave(std_dev_plot, filename = 'figureS3cd.pdf', width = 3, height = 10, units = "in", dpi = 300)
# ggsave(std_dev_plot, filename = 'figureS3cd.jpg', width = 3, height = 10, units = "in", dpi = 300)

wrapped_fitness_plot <- wrap_elements(full = fitness_plot)
wrapped_std_dev_plot <- wrap_elements(full = std_dev_plot)

# Combine the two plots
combined_plot <- wrapped_fitness_plot + wrapped_std_dev_plot + 
  plot_layout(ncol = 2, widths = c(8, 3))

combined_plot

ggsave(combined_plot, filename = 'figureS3.pdf', width = 11, height = 10, units = 'in', dpi = 300)
ggsave(combined_plot, filename = 'figureS3.jpg', width = 11, height = 10, units = 'in', dpi = 300)