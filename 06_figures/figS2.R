library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyverse)
library(ggcorrplot)
library(patchwork)
library(Cairo)

# plotting theme:
theme_gge <- function (base_size = 15, base_family = "sans")
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title = element_text(size = 7),
      panel.background = element_rect(fill = "white", colour = NA),
      #panel.border = element_rect(fill = NA, colour = "black"),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", size = 0.25),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      axis.ticks = element_blank(),
      legend.key = element_rect(fill = "white", colour = NA),
      complete = TRUE)
}

# Define kk_pop_ids to remove
excluded_kk_pop_ids <- c("a1.1_H06", "a3_H02", "a3_A03")

# population fitness data
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

# plot population data:
pop_pcplot <- ggplot(PD4, aes(x = env, y = fitness_gain_mu_pop, col = regime, group = kk_pop_id)) + geom_line(alpha= 0.5, size = 0.75) + theme_gge() +
  scale_color_manual(values = c('steelblue','darkorange2')) +
  scale_y_continuous(limits = c(-0.15,0.15), breaks = seq(-0.2,0.15,0.05)) +
  scale_x_discrete(expand = c(0.075, 0)) +  # Reduce white space on either side
  theme(
    #axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.position = c(0.9,0.9),
    legend.title = element_blank()
  ) +
  geom_hline(yintercept = 0, lty = 'dashed',col='gray40',alpha=0.5) +
  xlab('') + ylab('Fitness Gain') #+
#coord_flip()
pop_pcplot

# ggsave(pop_pcplot, filename = 'figureS2a.pdf', width = 8, height = 5, units = 'in', dpi = 300)
# ggsave(pop_pcplot, filename = 'figureS2a.jpg', width = 8, height = 5, units = 'in', dpi = 300)

# # Pivot to wide format: one row per population, one column per environment
# fitness_wide <- PD4 %>%
#   select(kk_pop_id, env, fitness_gain_mu_pop) %>%
#   pivot_wider(names_from = env, values_from = fitness_gain_mu_pop)
# 
# # Drop the kk_pop_id column to get a pure numeric matrix for correlation
# fitness_matrix <- fitness_wide %>% select(-kk_pop_id)
# 
# # Compute correlation matrix
# cor_matrix <- cor(fitness_matrix, use = "pairwise.complete.obs", method = "spearman")
# env_order <- c("YPD", "30°C", "37°C", "pH 7.3", "+NaCl", "-P")
# cor_matrix <- cor_matrix[env_order, env_order]
# 
# # Plot correlation heatmap
# ggcorrplot(cor_matrix, 
#            lab = TRUE, 
#            type = "upper",
#            legend.title = "Spearman\nCorrelation",
#            colors = c("black", "white", "red")) +
#   scale_y_discrete(limits = rev) +
#   theme(
#     axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
#     panel.grid = element_blank(),
#     panel.border = element_blank(),
#     panel.background = element_blank(),
#     axis.ticks = element_blank(),
#     legend.position = c(0.85, 0.75),
#     legend.title = element_text(size = 12),
#     legend.text = element_text(size = 10),
#     legend.key.height = unit(1.2, "cm"),                  # taller color bar
#     legend.key.width = unit(1, "cm")                    # wider color bar
#   )

# Define the environment order
env_order <- c("YPD", "30°C", "37°C", "pH 7.3", "+NaCl", "-P")

# Function to create the plot for a regime
plot_corr_by_regime <- function(data, regime_label) {
  fitness_wide <- data %>%
    filter(regime == regime_label) %>%
    select(kk_pop_id, env, fitness_gain_mu_pop) %>%
    pivot_wider(names_from = env, values_from = fitness_gain_mu_pop)
  
  fitness_matrix <- fitness_wide %>% select(-kk_pop_id)
  
  cor_matrix <- cor(fitness_matrix, use = "pairwise.complete.obs", method = "spearman")
  cor_matrix <- cor_matrix[env_order, env_order]
  
  ggcorrplot(cor_matrix, 
             lab = TRUE, 
             type = "upper",
             legend.title = "Spearman\nCorrelation",
             colors = c("black", "white", "red")) +
    scale_y_discrete(limits = rev) +
    labs(title = paste(regime_label)) +  # add plot title
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.position = c(0.9, 0.75),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.6, "cm"),
      legend.key.width = unit(0.5, "cm")
    )
}

# Generate the plots
plot_asexual <- plot_corr_by_regime(PD4, "asexual")
plot_sexual  <- plot_corr_by_regime(PD4, "sexual")

# Combine side by side
corr_plots <- plot_asexual + plot_spacer() + plot_sexual + plot_layout(widths = c(1,0.01,1))
corr_plots

combined_plot <- wrap_elements(pop_pcplot) + wrap_elements(corr_plots) +
  plot_layout(nrow = 2, heights = c(1.2,1)) +
  plot_annotation(tag_levels = 'A', tag_prefix = '', tag_suffix = '') &
  theme(
    plot.tag = element_text(face = 'bold', size = 25),
    plot.tag.position = c(0.01, 0.98),
    plot.margin = margin(5, 5, 5, 5)
  )
combined_plot

ggsave(combined_plot, filename = 'figureS2.pdf', width = 8, height = 9, units = 'in', dpi = 300)
ggsave(combined_plot, filename = 'figureS2.jpg', width = 8, height = 9, units = 'in', dpi = 300)

### STATS ###

# Spearman tests for each pair of envs
get_correlation_table <- function(data, regime_label) {
  fitness_wide <- data %>%
    filter(regime == regime_label) %>%
    select(kk_pop_id, env, fitness_gain_mu_pop) %>%
    pivot_wider(names_from = env, values_from = fitness_gain_mu_pop)
  
  fitness_matrix <- fitness_wide %>% select(-kk_pop_id)
  envs <- colnames(fitness_matrix)
  
  # If fewer than 2 environments, skip
  if (length(envs) < 2) {
    return(data.frame())
  }
  
  results <- list()
  
  for (i in seq_along(envs)) {
    for (j in seq(i + 1, length(envs))) {
      e1 <- envs[i]
      e2 <- envs[j]
      if (is.na(e1) || is.na(e2)) next  # skip invalid pairs
      
      vals <- na.omit(fitness_matrix[, c(e1, e2)])
      if (nrow(vals) >= 3) {
        test <- cor.test(vals[[1]], vals[[2]], method = "spearman", exact = FALSE)
        results[[length(results) + 1]] <- data.frame(
          regime = regime_label,
          env1 = e1,
          env2 = e2,
          rho = test$estimate,
          p_value = test$p.value
        )
      }
    }
  }
  
  if (length(results) == 0) return(data.frame())
  
  results_df <- bind_rows(results)
  results_df$p_adj_BH <- p.adjust(results_df$p_value, method = "BH")
  return(results_df)
}
cor_table_asexual <- get_correlation_table(PD4, "asexual")
cor_table_sexual  <- get_correlation_table(PD4, "sexual")
cor_table_all <- bind_rows(cor_table_asexual, cor_table_sexual)
print(cor_table_all)
#write.csv(cor_table_all, "figureS2_spearman_corrs_with_BH.csv", row.names = FALSE)

# Permutation test across all environment pairs
spearman_env_permutation_test <- function(fitness_matrix, n_perm = 1000, seed = 42) {
  set.seed(seed)
  
  # Step 1: Compute observed mean Spearman correlation
  observed_cor <- cor(fitness_matrix, method = "spearman", use = "pairwise.complete.obs")
  observed_mean_cor <- mean(observed_cor[upper.tri(observed_cor)], na.rm = TRUE)
  
  # Step 2: Permute columns and compute null distribution
  permuted_mean_cors <- numeric(n_perm)
  for (i in 1:n_perm) {
    permuted_matrix <- apply(fitness_matrix, 2, sample)  # permute within each environment
    cor_perm <- cor(permuted_matrix, method = "spearman", use = "pairwise.complete.obs")
    permuted_mean_cors[i] <- mean(cor_perm[upper.tri(cor_perm)], na.rm = TRUE)
  }
  
  # Step 3: Compute empirical p-value (one-sided, greater)
  p_val <- mean(permuted_mean_cors >= observed_mean_cor)
  
  # Return results
  list(
    observed_mean_correlation = observed_mean_cor,
    permuted_mean_correlations = permuted_mean_cors,
    p_value = p_val
  )
}
# Helper to prepare matrix
prepare_fitness_matrix <- function(data, regime_label) {
  fitness_wide <- data %>%
    filter(regime == regime_label) %>%
    select(kk_pop_id, env, fitness_gain_mu_pop) %>%
    pivot_wider(names_from = env, values_from = fitness_gain_mu_pop)
  
  fitness_matrix <- fitness_wide %>% select(-kk_pop_id)
  return(fitness_matrix)
}
# Prepare matrices
matrix_asexual <- prepare_fitness_matrix(PD4, "asexual")
matrix_sexual  <- prepare_fitness_matrix(PD4, "sexual")
# Run permutation tests
spearman_test_asexual <- spearman_env_permutation_test(matrix_asexual, n_perm = 10000)
spearman_test_sexual  <- spearman_env_permutation_test(matrix_sexual, n_perm = 10000)
# Print results
cat("\nAsexual:\n")
cat("Observed Mean Correlation: ", round(spearman_test_asexual$observed_mean_correlation, 3), "\n")
cat("Permutation p-value: ", spearman_test_asexual$p_value, "\n")
cat("\nSexual:\n")
cat("Observed Mean Correlation: ", round(spearman_test_sexual$observed_mean_correlation, 3), "\n")
cat("Permutation p-value: ", spearman_test_sexual$p_value, "\n")

# Permutation test plots

# Combine into tidy format
df_nulls <- tibble(
  correlation = c(spearman_test_asexual$permuted_mean_correlations,
                  spearman_test_sexual$permuted_mean_correlations),
  regime = rep(c("asexual", "sexual"), each = length(spearman_test_asexual$permuted_mean_correlations))
)
# Add observed values for vertical lines
obs_vals <- tibble(
  regime = c("asexual", "sexual"),
  obs = c(spearman_test_asexual$observed_mean_correlation,
          spearman_test_sexual$observed_mean_correlation)
)
obs_vals <- obs_vals %>%
  mutate(
    label = paste0("rho[plain('obs')] == ", round(obs, 3)),
    y_pos = Inf  # we'll use this to position text at the top
  )
# Factor for plotting order
df_nulls$regime <- factor(df_nulls$regime, levels = c("asexual", "sexual"))
obs_vals$regime <- factor(obs_vals$regime, levels = c("asexual", "sexual"))
null_plot <- ggplot(df_nulls, aes(x = correlation, fill = regime)) +
  geom_histogram(bins = 50, alpha = 0.8, color = "black") +
  geom_vline(data = obs_vals, aes(xintercept = obs), 
             color = "red", linetype = "dashed", linewidth = 1) +
  geom_text(
    data = obs_vals,
    aes(x = obs, y = y_pos, label = label),
    inherit.aes = FALSE,
    color = "red",
    size = 4, parse = TRUE,
    hjust = 1.1,  # push slightly to the right of the line
    vjust = 1.2    # push just below top edge
  ) +
  facet_wrap(~regime, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("asexual" = "steelblue", "sexual" = "darkorange2")) +
  theme_gge() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14),
    legend.position = "none",
    plot.title = element_text(size = 14),
  ) +
  labs(
    x = "Mean Spearman Correlation (Across Environments)",
    y = "Count",
    title = "Null Distributions of Mean Spearman Correlation"
  )
null_plot
#ggsave("null_distribution_spearman.jpg", plot = null_plot, width = 6, height = 6, units = "in", dpi = 300)
### END STATS ###