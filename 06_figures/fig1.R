library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)
library(tidyr)
library(lawstat)
library(Cairo)
library(magick)

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

# Define kk_pop_ids to remove -- these are the 6 clones with suspected whole-genome duplications in fig S7
excluded_kk_pop_ids <- c("a1.1_H06", "a3_H02", "a3_A03")

# load aggregated population fitness data
PD4 <- read.csv(file = 'parsed_pop_fitness_data.csv',T)
PD4$regime <- 'asexual'
PD4$regime[grep('a3',PD4$kk_pop_id)] <- 'sexual'
PD4$env <- factor(PD4$env, levels = c('YPD','SC30C','SC37C','SC_pH7.3','SC_0.2M_NaCl','lowP'))
PD4 <- PD4 %>%
  mutate(env = factor(env, 
                      levels = c('YPD', 'SC30C', 'SC37C', 'SC_pH7.3', 'SC_0.2M_NaCl', 'lowP'),
                      labels = c('YPD', '30°C', '37°C', 'pH 7.3', '+NaCl', '-P')))
PD4 <- PD4 %>%
  filter(!kk_pop_id %in% excluded_kk_pop_ids)

# define ranges for plotting lines:
ranges <- dplyr::group_by(PD4[PD4$fitness_gain_mu_pop<0.2,], env, regime) %>% summarise(mn = min(fitness_gain_mu_pop),mx = max(fitness_gain_mu_pop))

pop_pointplot <- ggplot(PD4[PD4$fitness_gain_mu_pop < 0.2, ], 
                        aes(x = regime, y = fitness_gain_mu_pop, color = regime)) +
  # Jittered points with slightly smaller size
  geom_jitter(width = 0.2, size = 1.8, alpha = 0.6) +
  # Horizontal line for the median, wider across the x-axis
  stat_summary(fun = median, geom = "crossbar", width = 0.8, size = 0.6, show.legend = FALSE) +  # Wider median line
  # Vertical line showing the range (min to max)
  geom_errorbar(stat = "summary", fun.min = min, fun.max = max, width = 0.1, size = 0.5, show.legend = FALSE) +
  annotate("rect", ymin = -Inf, ymax = 0, xmin = -Inf, xmax = Inf, alpha = 0.1, fill = "gray50") +
  geom_hline(yintercept = 0, lty = 'dashed',col='gray40',alpha=0.5) +
  facet_wrap(~ env, scales = 'fixed', nrow = 1, strip.position = "bottom") +  # Side-by-side subplots
  scale_y_continuous(limits = c(-0.15,0.15), breaks = seq(-0.2,0.15,0.05)) +
  theme_gge() +
  theme(
    legend.position = c(0.85, 0.95),
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    strip.background = element_blank(),  # Remove background for facet labels
    strip.text = element_text(size = 14, margin = margin(t = 10, b = 5)),  # Make facet labels bigger
    panel.spacing = unit(1, "lines")  # Add spacing between panels if needed
  ) +
  scale_color_manual(values = c('steelblue', 'darkorange2')) +  # Colors for jittered points and median lines
  ylab("Fitness Gain") + xlab("") +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 5, alpha = 1)))

fitness_plot <- pop_pointplot

# Test each distribution for normality
shapiro_results <- PD4 %>%
  filter(fitness_gain_mu_pop < 0.2) %>%
  group_by(env, regime) %>%
  summarise(
    n = n(),
    W = shapiro.test(fitness_gain_mu_pop)$statistic,
    p_value = shapiro.test(fitness_gain_mu_pop)$p.value
  ) %>%
  arrange(env, regime) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))
print(shapiro_results)

# Brunner-Munzel Tests 
bm_results <- PD4 %>%
  filter(fitness_gain_mu_pop < 0.2) %>%
  group_by(env) %>%
  summarise(
    n_sexual = sum(regime == "sexual"),
    n_asexual = sum(regime == "asexual"),
    test = list(brunner.munzel.test(
      fitness_gain_mu_pop[regime == "sexual"],
      fitness_gain_mu_pop[regime == "asexual"]
    ))
  ) %>%
  mutate(
    stat = sapply(test, function(x) x$statistic),
    df = sapply(test, function(x) x$parameter),
    p_value = sapply(test, function(x) x$p.value),
    p_adj_BH = p.adjust(p_value, method = "BH")
  ) %>%
  select(env, n_sexual, n_asexual, stat, df, p_value, p_adj_BH)

print(bm_results)

#### VARIANCE CALCULATIONS AND PLOTS ####

var_g <- dplyr::group_by(PD4[PD4$fitness_gain_mu_pop<0.2,], regime, env) %>% summarise(var_g = sqrt(var(fitness_gain_mu_pop)))
var_g$env <- factor(var_g$env, levels = c('YPD', '30°C', '37°C', 'pH 7.3', '+NaCl', '-P'))
var_g$regime <- factor(var_g$regime, levels = c('asexual','sexual'))
pop_var_g_p1 <- ggplot(var_g, aes(x = regime, y = var_g, group = env, col = env)) + 
  geom_line(color = "black") +
  geom_point(aes(fill = regime), shape = 21, size = 3, color = "black", stroke = 1) + 
  labs(y = expression(SD[fitness~gain]), x = NULL) +
  theme_gge() +
  theme(
    legend.position = 'none',
    legend.title = element_blank(),
    legend.text = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 14), 
    axis.text.x = element_text(size = 12, color = "black"),
    axis.line.x = element_blank(),
    plot.margin = margin(t=5, b = 25, r = 15, l = 10)
  ) +
  scale_fill_manual(values = c("sexual" = "darkorange2", "asexual" = "steelblue")) +
  scale_y_continuous(limits = c(0,0.04), breaks = seq(0,0.04,0.01)) +
  scale_x_discrete(labels = c("asexual" = "Asexual", "sexual" = "Sexual")) +
  geom_text(
    data = var_g %>% filter(regime == "sexual"),  # Filter only sexual points
    aes(label = env, y = var_g + ifelse(env == "pH 7.3", 0.0005, ifelse(env == "37°C", -0.0005, ifelse(env == "-P", -0.0005, 0)))),  # Use env as the label
    position = position_nudge(x = 0.15),
    hjust = 0,
    size = 3,  # Text size
    color = "black"  # Text color
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21, stroke = 0, size = 5)))

# Test each var_g (sexual and asexual) for normality with Shapiro-Wilk
asexual_var_g <- var_g %>% filter(regime == "asexual") %>% pull(var_g)
sexual_var_g  <- var_g %>% filter(regime == "sexual")  %>% pull(var_g)
shapiro_asexual <- shapiro.test(asexual_var_g)
shapiro_sexual  <- shapiro.test(sexual_var_g)
shapiro_var_g_results <- tibble(
  regime = c("asexual", "sexual"),
  W = c(shapiro_asexual$statistic, shapiro_sexual$statistic),
  p_value = c(shapiro_asexual$p.value, shapiro_sexual$p.value)
) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))
print(shapiro_var_g_results)

# Paired samples t test for var_g values
var_g_wide <- var_g %>%
  pivot_wider(names_from = regime, values_from = var_g)
t_test_var_g <- t.test(var_g_wide$asexual, var_g_wide$sexual, paired = TRUE)
print(t_test_var_g)

wrapped_fitness_plot <- wrap_elements(full = fitness_plot)
wrapped_pop_var_g_p1 <- wrap_elements(full = pop_var_g_p1)

# Read in the schematic PDF and make it a plot
schematic_img  <- image_read_pdf("fig1_schematic.pdf")[1]
schematic_grob <- grid::rasterGrob(schematic_img, interpolate = TRUE)
schematic_plot <- wrap_elements(schematic_grob)

# left column: schematic (top) over B (bottom)
# right column: A spanning both rows
combined_plot <- ((schematic_plot / wrapped_pop_var_g_p1) | wrapped_fitness_plot) +
  plot_layout(widths = c(4, 8)) +
  plot_annotation(
    tag_levels = list(c("A", "C", "B")),
    tag_prefix = "",
    tag_suffix = ""
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 25),
    plot.tag.position = c(0.01, 0.98),
    plot.margin = margin(5, 0, 0, 5)
  )

combined_plot

ggsave(combined_plot, filename = "figure1.pdf", width = 8, height = 6, units = "in", dpi = 300, device = cairo_pdf)
ggsave(combined_plot, filename = "figure1.jpg", width = 8, height = 6, units = "in", dpi = 300)
