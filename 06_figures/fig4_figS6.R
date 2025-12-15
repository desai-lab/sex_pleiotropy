library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)
library(scales)
library(Cairo)
library(grid)
library(lme4)
library(lmerTest)
library(lawstat)
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


# import fitness (mean-centered) data for parents from fa03:
comb1 <- read.csv(file = 'parental_mc_fitnesses_fa03_fa08.csv')

# import F1 information for FA06, FA07, and FA08:
fa06_f1 <- read.csv(file = 'fa06_F1_final.csv')
fa06_p  <- read.csv(file = 'fa06_Parents_final.csv')

fa07_f1 <- read.csv(file = 'fa07_F1_final.csv')
fa07_p  <- read.csv(file = 'fa07_Parents_final.csv')

fa08_f1 <- read.csv(file = 'fa08_all_F1_final.csv')
fa08_p  <- read.csv(file = 'fa08_all_Parents_final.csv')

# merge together by relevant columns:
focal_cols <- c('kk_strain_id', 'f1_unique_id', 'env', 'regime', 'fitness_gain_adj')

fa06_f1 <- fa06_f1[,focal_cols]
fa07_f1 <- fa07_f1[,focal_cols]
fa08_f1 <- fa08_f1[,focal_cols]

# add fa label:
fa06_f1$assay <- 'FA06'
fa07_f1$assay <- 'FA07'
fa08_f1$assay <- 'FA08'

# rbind together:
fa_all <- bind_rows(fa06_f1,fa07_f1,fa08_f1)

# now combine parental fitnesses in just the same way:
focal_cols_p <- c('kk_strain_id', 'env', 'regime', 'fitness_gain_adj')
fa06_p <- fa06_p[,focal_cols_p]
fa07_p <- fa07_p[,focal_cols_p]
fa08_p <- fa08_p[,focal_cols_p]

fa06_p$assay <- 'FA06'
fa07_p$assay <- 'FA07'
fa08_p$assay <- 'FA08'

fa_all_p <- bind_rows(fa06_p,fa07_p,fa08_p)

# at this point, the parents are aligned with their F1 progeny.
# now I need to calculate differences between each of these parents and the fitness of them via fa03:

fa_all_p2 <- merge(fa_all_p,comb1[,!names(comb1) %in% c('fitness_gain_adj')], by = c('kk_strain_id','env','regime'))

# calculate differential. This value will be subtracted from both parents and their F1s:
#fa_all_p2$diff <- fa_all_p2$fitness_gain_adj - fa_all_p2$mc_avg

names(fa_all_p2)[which(names(fa_all_p2)=='fitness_gain_adj')] <- 'parental_fitness'

# add back diff by merging diff by merging on kk_strain_id, env, and regime:
fa_all2 <- merge(fa_all, fa_all_p2, by = c('kk_strain_id','env','regime','assay'))

# subtract off parental fitness. Now all F1s should be expressed relative to their parents (who are at zero)
fa_all2$parent_centered_fitness <- fa_all2$fitness_gain_adj - fa_all2$parental_fitness

# now add parental fitness from fa03 to string them back on a line:
fa_all2$re_centered_fitness <- fa_all2$parent_centered_fitness + fa_all2$mc_avg

# now let's try to plot. First I need to sort strain sets by parental fitness rank per env:
order1 <- dplyr::filter(fa_all_p2, env == 'YPD') %>% arrange(desc(regime), desc(mc_avg))
order2 <- unique(paste0(order1$kk_strain_id))
order2 <- order2[order2 %in% fa_all2$kk_strain_id]
fa_all2$kk_strain_id <- factor(fa_all2$kk_strain_id, levels = order2)
fa_all2$regime <- factor(fa_all2$regime, levels = c('sexual','asexual'))

# sort parents:
fa_all_p2 <- dplyr::filter(fa_all_p2, kk_strain_id %in% order2)
fa_all_p2$kk_strain_id <- factor(fa_all_p2$kk_strain_id, levels = order2)

# Define the strings to filter out
exclude_patterns <- c("D04_a1.1_H06", "C04_a3_H02", "C02_a3_A03")

# Create a single regular expression pattern
exclude_regex <- paste(exclude_patterns, collapse = "|")

# Filter out rows in fa_all2 where kk_strain_id contains any of the exclude patterns
fa_all2 <- fa_all2[!grepl(exclude_regex, fa_all2$kk_strain_id), ]
fa_all2$kk_strain_id <- droplevels(fa_all2$kk_strain_id)

# Filter out rows in fa_all_p2 where kk_strain_id contains any of the exclude patterns
fa_all_p2 <- fa_all_p2[!grepl(exclude_regex, fa_all_p2$kk_strain_id), ]
fa_all_p2$kk_strain_id <- droplevels(fa_all_p2$kk_strain_id)

# now plot, for YPD only:
F1_parent_dists1 <- ggplot(dplyr::filter(fa_all2, env=='YPD')) +
  #geom_point(data = dplyr::filter(fa_all_p2, env=='YPD'), aes(x = kk_strain_id, y = fitness_gain - diff), col = 'white',fill='white',size=3) +
  geom_jitter(width = 0.1, alpha = 0.5, aes(x = kk_strain_id, y = re_centered_fitness, col = regime), size = 0.5) +
  geom_boxplot(alpha = 0.8, aes(x = kk_strain_id, y = re_centered_fitness, col = regime),outlier.colour = NA) +
  theme_gge() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = 'none') +
  geom_point(data = dplyr::filter(fa_all_p2, env=='YPD'), aes(x = kk_strain_id, y = mc_avg), col = 'black',fill='white',size=2) +
  scale_color_manual(values = c('darkorange2','steelblue')) +
  facet_grid(~ regime, space = 'free', scales = 'free_x') +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(limits = c(-0.2,0.15), breaks = seq(-0.2,0.15,0.05))

# compute quantiles per regime:
# split up into kk_strain_id:
sets <- split(fa_all2, fa_all2$kk_strain_id)
# x <- sets[[1]]
# the_env <- 'YPD'
compute_quantile <- function(x,the_env){
  x2 <- dplyr::filter(x, env == the_env)
  # calculate quantile of parent in F1 distr:
  ps <- unique(x2$mc_avg)
  fit_sort <- sort(x2$re_centered_fitness)
  pq <- sum(fit_sort<=ps) / length(fit_sort)
  # also compute F1 summary statistics:"
  f1_quantiles <- data.frame(t(quantile(fit_sort)))
  f1_mean = mean(fit_sort, na.rm=T)
  d1 <- data.frame(kk_strain_id = x2$kk_strain_id[1],
                   env = x2$env[1],
                   regime = x2$regime[1],
                   mc_avg = x2$mc_avg[1],
                   pq = pq,
                   f1_mean = f1_mean,
                   f1_q0 = f1_quantiles[1],
                   f1_q0.25 = f1_quantiles[2],
                   f1_q0.50 = f1_quantiles[3],
                   f1_q0.75 = f1_quantiles[4],
                   f1_q100 = f1_quantiles[5]
  )
  return(d1)
}

sets2 <- lapply(sets, function(x) compute_quantile(x, the_env='YPD'))
sets3 <- do.call(rbind, sets2)
sets3$regime <- factor(sets3$regime, levels = c('sexual','asexual'))

mut_dat2 <- read.csv('mut_dat2.csv')
# change case of clone_id:
mut_dat2$clone_id <- sapply(mut_dat2$clone_id, function(x) gsub('_A1.1_','_a1.1_',paste0(x)))
mut_dat2$clone_id <- sapply(mut_dat2$clone_id, function(x) gsub('_A3_','_a3_',paste0(x)))

# combine with F1 summary information:
sets4 <- merge(sets3, mut_dat2[,!names(mut_dat2) %in% 'regime'], by.x = 'kk_strain_id', by.y = 'clone_id')

# Define the annotation text for specific populations
sets3$kk_strain_id <- as.character(sets3$kk_strain_id)

# Define the annotation text for specific populations
sets3$annotation <- NA  # Add a new column for annotations
# Function to check if a string starts with any value from a given vector
startsWithAny <- function(x, patterns) {
  Reduce("|", lapply(patterns, startsWith, x = x))
}

# Assign values for the "sexual" regime
sexual_patterns <- c("D03", "E04", "G04", "A04")
sexual_matches <- startsWithAny(sets3$kk_strain_id, sexual_patterns)

sets3$annotation[sets3$regime == "sexual" & sexual_matches] <- 
  match(substr(sets3$kk_strain_id[sets3$regime == "sexual" & sexual_matches], 1, 3), sexual_patterns)

# Assign values for the "asexual" regime
asexual_patterns <- c("A03", "B05", "D05", "D06")
asexual_matches <- startsWithAny(sets3$kk_strain_id, asexual_patterns)

sets3$annotation[sets3$regime == "asexual" & asexual_matches] <- 
  match(substr(sets3$kk_strain_id[sets3$regime == "asexual" & asexual_matches], 1, 3), asexual_patterns)

sets3$kk_strain_id <- as.factor(sets3$kk_strain_id)

# Add a new column for jittered x-positions to ensure alignment
set.seed(20)  # Set seed for reproducibility
sets3$regime <- factor(sets3$regime, levels = c("asexual", "sexual"))
sets3$jittered_x <- jitter(as.numeric(as.factor(sets3$regime)), amount = 0.2) 

# Compute overall median for each regime
median_values <- sets3 %>%
  group_by(regime) %>%
  summarize(median_fitness = median(f1_mean - mc_avg, na.rm = TRUE))

# now combine mean_diffs1 and mean_diffs2 for Fig5 plot.
mean_diffs1c <- ggplot(sets3, aes(x = jittered_x, y = f1_mean - mc_avg, fill = regime)) + 
  geom_segment(data = median_values, 
               aes(x = ifelse(regime == "asexual", 0.7, 1.7),  # Start at 0.7 for sexual, 1.7 for asexual
                   xend = ifelse(regime == "asexual", 1.3, 2.3),  # End at 1.3 for sexual, 2.3 for asexual
                   y = median_fitness, 
                   yend = median_fitness),
               size = 3, color = 'black') +
  geom_point(shape = 21, size = 2.5, color = "black", stroke = 1) +
  scale_fill_manual(values = c('steelblue', 'darkorange2')) +
  scale_color_manual(values = c('steelblue','darkorange2')) +
  theme_gge() + geom_hline(yintercept = 0, alpha = 0.5) +
  geom_text(aes(label = annotation, color = regime, x = jittered_x), hjust = -0.8, vjust = -0.2, size = 4, na.rm = TRUE) +  # Add annotations
  theme(legend.position = 'none', strip.text = element_blank(), axis.title.y = element_text(size = 14), axis.text.x = element_text(size = 13, color = "black", vjust = -3), axis.line.x = element_blank()) +
  facet_wrap(~ env, ncol = 2) +
  scale_y_continuous(limits = c(-0.04, 0.025),breaks = seq(-0.04,0.02,0.01)) +
  scale_x_continuous(breaks = c(1, 2), labels = c("asexual", "sexual")) +
  labs(y = "Mean F1 − Parent Fitness", x = NULL) +
  annotate("rect", ymin = -Inf, ymax = 0, xmin = -Inf, xmax = Inf, alpha = 0.1, fill = "gray50") +
  annotate("segment", x = 2.35, xend = 2.35, y = -0.036, yend = 0, col = "grey70", linetype = "solid", size = 0.5) +
  annotate("segment", x = 2.25, xend = 2.35, y = -0.036, yend = -0.036, col = "grey70", linetype = "solid", size = 0.5) +
  annotate("segment", x = 2.25, xend = 2.35, y = 0, yend = 0, col = "grey70", linetype = "solid", size = 0.5) +
  annotate("text", x = 2.2, y = -0.0385, label = "78% ≤ 0", hjust = 1, size = 3.5, color = "grey30") +
  annotate("segment", x = 0.6, xend = 0.6, y = 0, yend = 0.02, col = "grey70", linetype = "solid", size = 0.5) +
  annotate("segment", x = 0.6, xend = 0.7, y = 0.02, yend = 0.02, col = "grey70", linetype = "solid", size = 0.5) +
  annotate("segment", x = 0.6, xend = 0.7, y = 0, yend = 0, col = "grey70", linetype = "solid", size = 0.5) +
  annotate("text", x = 1.15, y = 0.0235, label = "60% ≥ 0", hjust = 1, size = 3.5, color = "grey30")

# Define the annotation text for specific populations
sets4$kk_strain_id <- as.character(sets4$kk_strain_id)

# Define the annotation text for specific populations
sets4$annotation <- NA  # Add a new column for annotations
prefix_map <- c("D03" = "1", "E04" = "2", "G04" = "3", "A04" = "4")
sets4$annotation[sets4$regime == "sexual"] <- sapply(
  sets4$kk_strain_id[sets4$regime == "sexual"],
  function(id) {
    prefix <- substr(id, 1, 3)
    prefix_map[prefix]
  }
)
sets4$annotation[sets4$regime == "asexual" & (startsWith(sets4$kk_strain_id, "A03") |
                                                startsWith(sets4$kk_strain_id, "B05") |
                                                startsWith(sets4$kk_strain_id, "D05") |
                                                startsWith(sets4$kk_strain_id, "D06"))] <- c("1", "2", "3", "4")

sets4$kk_strain_id <- as.factor(sets4$kk_strain_id)

### STATS ###

# Create binary outcome: is F1 mean > parent fitness?
sets3$above_zero <- ifelse((sets3$f1_mean - sets3$mc_avg) > 0, 1, 0)
# Create contingency table
contingency <- table(sets3$regime, sets3$above_zero)
# View table
print(contingency)
# Run Fisher's Exact Test
fisher_result <- fisher.test(contingency)
print(fisher_result)

sexual_vals <- dplyr::filter(sets3, regime == "sexual")$f1_mean - dplyr::filter(sets3, regime == "sexual")$mc_avg
asexual_vals <- dplyr::filter(sets3, regime == "asexual")$f1_mean - dplyr::filter(sets3, regime == "asexual")$mc_avg

# Perform Shapiro-Wilk test
shapiro_sexual <- shapiro.test(sexual_vals)
shapiro_asexual <- shapiro.test(asexual_vals)
shapiro_df <- data.frame(
  regime = c("sexual", "asexual"),
  n = c(length(sexual_vals), length(asexual_vals)),
  W = c(shapiro_sexual$statistic, shapiro_asexual$statistic),
  raw_p = c(shapiro_sexual$p.value, shapiro_asexual$p.value)
)
# Add BH-adjusted p-values
shapiro_df$BH_adj_p <- p.adjust(shapiro_df$raw_p, method = "BH")
print(shapiro_df)

# Welch's t-test (default in t.test for two groups)
welch_test <- t.test(sexual_vals, asexual_vals)
print(welch_test)

# Mann-Whitney U test (Wilcoxon rank-sum test)
mw_test <- wilcox.test(sexual_vals, asexual_vals, exact = FALSE)  # exact=FALSE if n > ~50
print(mw_test)

# Perform one-sample t-tests
t_sexual <- t.test(sexual_vals, mu = 0)
t_asexual <- t.test(asexual_vals, mu = 0)
t_test_df <- data.frame(
  regime = c("sexual", "asexual"),
  n = c(length(sexual_vals), length(asexual_vals)),
  t_statistic = c(t_sexual$statistic, t_asexual$statistic),
  df = c(t_sexual$parameter, t_asexual$parameter),
  p_value = c(t_sexual$p.value, t_asexual$p.value)
)
# Apply BH correction (optional but consistent)
t_test_df$BH_adj_p <- p.adjust(t_test_df$p_value, method = "BH")
print(t_test_df)

# Perform Wilcoxon signed-rank tests
w_sexual <- wilcox.test(sexual_vals, mu = 0)
w_asexual <- wilcox.test(asexual_vals, mu = 0)
wilcox_df <- data.frame(
  regime = c("sexual", "asexual"),
  n = c(length(sexual_vals), length(asexual_vals)),
  W_statistic = c(w_sexual$statistic, w_asexual$statistic),
  p_value = c(w_sexual$p.value, w_asexual$p.value)
)
wilcox_df$BH_adj_p <- p.adjust(wilcox_df$p_value, method = "BH")
print(wilcox_df)

# Filter to just YPD environment
fa_ypd <- dplyr::filter(fa_all2, env == 'YPD')
# Mixed model: F1 − Parent Fitness ~ Regime + (random intercept per parent strain)
model_centered <- lmer(parent_centered_fitness ~ regime + (1 | kk_strain_id), data = fa_ypd)
# View results
summary(model_centered)

# # Order strains by average fitness or regime if needed
# fa_ypd$kk_strain_id <- factor(fa_ypd$kk_strain_id, levels = unique(fa_ypd$kk_strain_id[order(fa_ypd$regime)]))

strain_medians <- fa_ypd %>%
  group_by(kk_strain_id) %>%
  summarize(median_fitness = median(parent_centered_fitness, na.rm = TRUE)) %>%
  arrange(desc(median_fitness))  # descending order

fa_ypd$kk_strain_id <- factor(fa_ypd$kk_strain_id, levels = strain_medians$kk_strain_id)

# Mini violins per parent
ggplot(fa_ypd, aes(x = kk_strain_id, y = parent_centered_fitness, fill = regime)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.4) +
  stat_summary(fun = median, geom = "point", shape = 21, fill = "white", color = "black", size = 2) +  # ← median points
  stat_summary(fun = mean, geom = "point", shape = 21, fill = "red", color = "black", size = 2) +
  scale_fill_manual(values = c("darkorange2", "steelblue")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_gge() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "none"
  ) +
  labs(x = "Strain", y = "F1 − Parent Fitness") +
  facet_grid(~ regime, scales = "free_x", space = "free")

### END STATS ###

# Plot with annotations
mean_diffs2b <- ggplot(sets4[sets4$hit_tot > 2, ], aes(x = single_count - multi_count, y = f1_mean - mc_avg, fill = regime, col = regime, group = kk_strain_id)) +
  geom_point(shape = 21, size = 2.5, color = "black", stroke = 1) +
  geom_text(aes(label = annotation), hjust = -0.8, vjust = -0.5, size = 4, na.rm = TRUE) +  # Add annotations
  scale_fill_manual(values = c('darkorange2', 'steelblue')) +
  scale_color_manual(values = c('darkorange2', 'steelblue')) +
  theme_gge() +
  geom_hline(yintercept = 0, alpha = 0.5) +
  annotate("rect", ymin = -Inf, ymax = 0, xmin = -Inf, xmax = Inf, alpha = 0.1, fill = "gray50") +
  theme(legend.position = 'none', axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), axis.line.x = element_blank()) +
  scale_y_continuous(limits = c(-0.04, 0.025), breaks = seq(-0.04, 0.02, 0.01)) +
  scale_x_continuous(breaks = scales::breaks_width(10), expand = expansion(mult = c(0.2, 0.2))) +
  labs(y = "Mean F1 − Parent Fitness", x = "Singleton − Multi-hit Mutation Count")

### STATS ###

# Subset to relevant data (already filtered for hit_tot > 2)
cor_data <- sets4[sets4$hit_tot > 2, ]
# Separate by regime
cor_sexual <- subset(cor_data, regime == "sexual")
cor_asexual <- subset(cor_data, regime == "asexual")
# Spearman tests
spearman_sexual <- cor.test(
  cor_sexual$single_count - cor_sexual$multi_count,
  cor_sexual$f1_mean - cor_sexual$mc_avg,
  method = "spearman"
)
spearman_asexual <- cor.test(
  cor_asexual$single_count - cor_asexual$multi_count,
  cor_asexual$f1_mean - cor_asexual$mc_avg,
  method = "spearman"
)
# View results
spearman_sexual
spearman_asexual

# Create mutation bias column in sets4
sets4$mutation_bias <- sets4$single_count - sets4$multi_count
# Merge into fa_all2 (clone-level data)
fa_all2_mut <- merge(fa_all2, sets4[, c("kk_strain_id", "mutation_bias", "regime")], by = "kk_strain_id")
# Subset to YPD and drop clones with missing mutation data (just to be safe)
fa_ypd_mut <- dplyr::filter(fa_all2_mut, env == "YPD" & !is.na(mutation_bias))
fa_ypd_mut$regime <- fa_ypd_mut$regime.x
# Full model: mutation bias + regime + interaction, random intercept for strain
model_mutbias <- lmer(parent_centered_fitness ~ mutation_bias * regime + (1 | kk_strain_id), data = fa_ypd_mut)
summary(model_mutbias)

# Sexual and asexual clones only in YPD
fa_ypd_sexual <- subset(fa_ypd_mut, regime == "sexual")
fa_ypd_asexual <- subset(fa_ypd_mut, regime == "asexual")
# Sexual model
model_sexual <- lmer(parent_centered_fitness ~ mutation_bias + (1 | kk_strain_id), data = fa_ypd_sexual)
summary(model_sexual)
# Asexual model
model_asexual <- lmer(parent_centered_fitness ~ mutation_bias + (1 | kk_strain_id), data = fa_ypd_asexual)
summary(model_asexual)


### END STATS ###


## Figure 4c (spearman correlations)

fa_all2_h <- fa_all2[fa_all2$env=='YPD',c('f1_unique_id','fitness_gain_adj')]
names(fa_all2_h)[2] <- 's_at_home'
fa_all3 <- merge(fa_all2[fa_all2$env!='YPD',], fa_all2_h, by = 'f1_unique_id')

# now do parents:
fa_all_p2_h <- fa_all_p2[fa_all_p2$env=='YPD',c('kk_strain_id','mc_avg')]
names(fa_all_p2_h)[2] <- 's_at_home'
fa_all_p3 <- merge(fa_all_p2[fa_all_p2$env!='YPD',], fa_all_p2_h, by = 'kk_strain_id')
fa_all_p3 <- dplyr::filter(fa_all_p3,kk_strain_id %in% unique(paste0(fa_all3$kk_strain_id)))

# calculate spearman ranks for asexual and sexual, for each environment:
x <- fa_all3[fa_all3$regime=='sexual',]

cor_by_env <- function(x,regime) {
  # exclude rows where there are any NA values:
  x <- droplevels(x[!(is.na(x$s_at_home) | is.na(x$re_centered_fitness)),])
  x2 <- split(x, list(x$env,x$kk_strain_id))
  # calculate cor for ecah env:
  cor_calc <- function(xx){
    x4 <- as.matrix(xx[,c('s_at_home','re_centered_fitness')])
    the_cor <- cor(x4[,1],x4[,2], method = 'spearman')
    # construct results row:
    data.frame(kk_strain_id = xx$kk_strain_id[1],
               regime = regime,
               env = xx$env[1],
               the_cor = as.vector(the_cor))
  }
  the_cors <- lapply(x2,function(x) cor_calc(x))
  do.call(rbind,the_cors)
}

scors <- cor_by_env(fa_all3[fa_all3$regime=='sexual',], regime = 'sexual')
acors <- cor_by_env(fa_all3[fa_all3$regime=='asexual',], regime = 'asexual')

all_cors <- bind_rows(scors,acors)

# connect the dots:
spearman_plot1 <- ggplot(all_cors, aes(x = env, y = the_cor, group = kk_strain_id, col = regime)) +
  geom_point(aes(fill = regime), shape = 21, size = 2.5, color = "black", stroke = 1) +
  geom_line() +
  #facet_wrap(~ regime) +
  labs( y = "Home-Away Fitness Correlation") +
  theme_gge() + 
  theme(
    legend.position = 'none',
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14, color = "black", vjust = -3),
    axis.line.x = element_blank()
  ) +
  scale_x_discrete(labels = c("SC37C" = "37°C", "SC_pH7.3" = "pH 7.3", "SC_0.2M_NaCl" = "+NaCl")) +
  scale_fill_manual(values = c('steelblue','darkorange2')) +
  scale_color_manual(values = c('steelblue','darkorange2')) +
  scale_y_continuous(limits = c(-0.2, 0.8), breaks = seq(-0.2, 0.8, 0.2)) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  annotate("rect", ymin = -Inf, ymax = 0, xmin = -Inf, xmax = Inf, alpha = 0.1, fill = "gray50")

### STATS

# split data into groups
sexuals  <- all_cors$the_cor[all_cors$regime == "sexual"]
asexuals <- all_cors$the_cor[all_cors$regime == "asexual"]

# Brunner–Munzel test
res_bm <- brunner.munzel.test(sexuals, asexuals, alternative = "two.sided")
print(res_bm)

### END STATS

ypd_f1_both <- mean_diffs1c + plot_spacer() + spearman_plot1 +
  plot_layout(ncol = 3, widths = c(0.6, 0.15, 1)) +
  plot_annotation(tag_levels = list(c("B", "D")), tag_prefix = '', tag_suffix = '') & 
  theme(plot.tag = element_text(face = 'bold', size = 25), plot.tag.position = c(0.01, 1.01), plot.margin = margin(7, 5, 5, 5))

## Figure 4d (fitness distributions)

# generate F1 fitness distributions for 4 asexual and 4 sexual strains. Also produce plot of ALL of the data for inspection.
# pick 4 of each:
# data.frame with all offspring fitnesses is:
YPD_fit <- dplyr::filter(fa_all2, env=='YPD')

asex_set <- grep('D06_|A03_|D05_|B05_',paste0(sets4[sets4$hit_tot>2 & sets4$regime=='asexual','kk_strain_id']),value=T)
# sex_set  <- grep('E04_|G04_|D03_|A04_',paste0(sets4[sets4$hit_tot>2 & sets4$regime=='sexual','kk_strain_id']),value=T)
sex_set <- c(
  grep('D03_', paste0(sets4[sets4$hit_tot > 2 & sets4$regime == 'sexual', 'kk_strain_id']), value = TRUE),
  grep('E04_', paste0(sets4[sets4$hit_tot > 2 & sets4$regime == 'sexual', 'kk_strain_id']), value = TRUE),
  grep('G04_', paste0(sets4[sets4$hit_tot > 2 & sets4$regime == 'sexual', 'kk_strain_id']), value = TRUE),
  grep('A04_', paste0(sets4[sets4$hit_tot > 2 & sets4$regime == 'sexual', 'kk_strain_id']), value = TRUE)
)

YPD_fita <- dplyr::filter(YPD_fit, kk_strain_id %in% asex_set)
YPD_fits <- dplyr::filter(YPD_fit, kk_strain_id %in% sex_set)

# grab matching parental fitnesses:
YPD_parents_a <- dplyr::filter(fa_all_p2, env=='YPD', kk_strain_id %in% asex_set)
YPD_parents_s <- dplyr::filter(fa_all_p2, env=='YPD', kk_strain_id %in% sex_set)

# grab unique values:
YPD_parents_a <- unique(YPD_parents_a[,c('kk_strain_id','mc_avg','regime')])
YPD_parents_s <- unique(YPD_parents_s[,c('kk_strain_id','mc_avg','regime')])

# re-order factors for plotting:
YPD_parents_a$kk_strain_id <- factor(YPD_parents_a$kk_strain_id, levels = asex_set)
YPD_parents_s$kk_strain_id <- factor(YPD_parents_s$kk_strain_id, levels = sex_set)

YPD_fit_both <- bind_rows(YPD_fita,YPD_fits)
YPD_parents <- bind_rows(YPD_parents_a,YPD_parents_s)

# generate new factor for facet_grid to work on:
YPD_fit_both$strain_fact <- match(YPD_fit_both$kk_strain_id,asex_set)
YPD_fit_both$strain_fact[is.na(YPD_fit_both$strain_fact)] <- match(YPD_fit_both$kk_strain_id,sex_set)[!is.na(match(YPD_fit_both$kk_strain_id,sex_set))]
YPD_fit_both$strain_fact <- factor(YPD_fit_both$strain_fact)

# add same factor to parents:
YPD_parents <- merge(YPD_parents, unique(YPD_fit_both[names(YPD_fit_both) %in% c('kk_strain_id','regimne','strain_fact')]), by = 'kk_strain_id')

# Filter data for sexual and asexual rows
sexual_data <- YPD_fit_both[YPD_fit_both$regime == "sexual", ]
asexual_data <- YPD_fit_both[YPD_fit_both$regime == "asexual", ]
sexual_parents <- YPD_parents[YPD_parents$regime == "sexual", ]
asexual_parents <- YPD_parents[YPD_parents$regime == "asexual", ]

# Determine shared x-axis limits based on the data
#x_limits <- range(YPD_fit_both$re_centered_fitness, na.rm = TRUE)
x_limits <- c(-0.1, 0.1)
y_limits <- c(0,30)

# Plot for sexual row
sexual_plot <- ggplot(sexual_data, aes(x = re_centered_fitness)) +
  geom_histogram(bins = 20, fill = 'gray40', alpha = 0.5) +
  facet_grid(. ~ strain_fact) +  # Only one row
  geom_vline(data = sexual_parents, aes(xintercept = mc_avg, col = regime), size = 1.5, alpha = 0.5) +
  theme_gge() + 
  theme(
    legend.position = 'none',
    strip.text.y = element_blank(), strip.text.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    plot.title = element_text(size = 10, hjust = 0.01, vjust = 2, face = "bold"),
    axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)
  ) +
  scale_color_manual(values = c('darkorange2')) +  # Color for sexuals
  scale_x_continuous(limits = x_limits) +  # Set consistent x-axis limits
  scale_y_continuous(limits = y_limits) +
  labs(title = "Sexual", y = "No. of F1 clones", x = "Mean-centered Fitness")  # Keep x-axis label for alignment

# Plot for asexual row
asexual_plot <- ggplot(asexual_data, aes(x = re_centered_fitness)) +
  geom_histogram(bins = 20, fill = 'gray40', alpha = 0.5) +
  facet_grid(. ~ strain_fact) +  # Only one row
  geom_vline(data = asexual_parents, aes(xintercept = mc_avg, col = regime), size = 1.5, alpha = 0.5) +
  theme_gge() + 
  theme(
    legend.position = 'none',
    strip.text.y = element_blank(), strip.text.x = element_text(size = 10, margin = margin(b = 10)),
    panel.spacing = unit(1, "lines"),
    plot.title = element_text(size = 10, hjust = 0.01, vjust = -6, face = "bold"),
    axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)
  ) +
  scale_color_manual(values = c('steelblue')) +  # Color for asexuals
  scale_x_continuous(limits = x_limits) +  # Set consistent x-axis limits
  scale_y_continuous(limits = y_limits) +
  labs(title = "Asexual", y = "No. of F1 clones", x = NULL)  # Remove x-axis label

tagged_asexual_plot <- asexual_plot + 
  labs(tag = "C")

F1_hists_plot <- tagged_asexual_plot / plot_spacer() / sexual_plot +
  plot_layout(guides = "collect", axis_titles = "collect", heights = c(1, 0.005, 1)) +
  labs(x = "Mean-centered Fitness", y = "No. of F1 clones") +
  plot_annotation(
    theme = theme(plot.margin = margin(5, 5, 5, 5))
  ) & 
  theme(plot.tag = element_text(face = 'bold', size = 25), plot.tag.position = c(0.01, 0.98), plot.margin = margin(5, 5, 5, 5))


wrapped_top_plot <- wrap_elements(ypd_f1_both)
wrapped_bottom_plot <- wrap_elements(F1_hists_plot)

combined_plot <- wrapped_top_plot / wrapped_bottom_plot +
  plot_layout(widths = c(1, 1)) +  
  plot_annotation(
    theme = theme(plot.margin = margin(5, 5, 5, 5))
  )

combined_plot

external_image <- image_read_pdf("GGE_backcross_figure.pdf")
external_grob <- grid::rasterGrob(external_image, interpolate = TRUE)
external_plot <- wrap_elements(external_grob)

final_plot <- external_plot / combined_plot +
  plot_layout(widths = c(1, 1), heights = c(0.6,2)) +
  plot_annotation(tag_levels = list(c("A", "")), tag_prefix = '', tag_suffix = '') &
  theme(plot.tag = element_text(face = 'bold', size = 25), plot.tag.position = c(0.015, 0.95), plot.margin = margin(3, 3, 3, 3))

final_plot

# Save the combined plot
ggsave(final_plot, filename = "figure4.pdf", width = 10, height = 15, units = "in", dpi = 300, device = cairo_pdf)
ggsave(final_plot, filename = "figure4.jpg", width = 10, height = 15, units = "in", dpi = 300)


### Figure S6 ###

quantile_flagger_2 <- function(x){
  
  # first split by environment:
  x2 <- reshape2::dcast(x, f1_unique_id ~ env, value.var = 're_centered_fitness')
  
  # remove rows with incomplete data:
  x2 <- x2[complete.cases(x2),]
  
  # scale first three cols:
  x2[,-c(1)] <- apply(x2[,-c(1)],2,scale)
  
  # now calculate average:
  x2$pleio_z <- apply(x2[,-c(1)],1,function(x) mean(x,na.rm=T))
  
  # add back other meta-data:
  x2$regime <- x$regime[1]
  x2$kk_strain_id <- x$kk_strain_id[1]
  
  return(x2)
}

# need to compile strains with pleiotropy information first:
plei_tmp <- rowSums(table(fa_all2$kk_strain_id[fa_all2$env!='YPD'],fa_all2$env[fa_all2$env!='YPD']))
all_pleio_strains <- names(plei_tmp[plei_tmp>0])
fa_plei <- droplevels(dplyr::filter(fa_all2, kk_strain_id %in% all_pleio_strains))

all_F1s <- split(fa_plei,fa_plei$kk_strain_id)
#x <- all_F1s[[1]]

all_F1s_2 <- lapply(all_F1s, quantile_flagger_2)
all_F1s_3 <- do.call(rbind, all_F1s_2)

# Determine global x and y limits
x_min <- min(all_F1s_3$YPD, na.rm = TRUE)
x_max <- max(all_F1s_3$YPD, na.rm = TRUE)
y_min <- min(all_F1s_3$pleio_z, na.rm = TRUE)
y_max <- max(all_F1s_3$pleio_z, na.rm = TRUE)

# Define consistent limits
x_limits <- c(x_min, x_max)
y_limits <- c(y_min, y_max)

# Adding spearman correlations and p values to plot
rho_labels <- all_F1s_3 %>%
  group_by(kk_strain_id, regime) %>%
  summarise(
    test = list(cor.test(YPD, pleio_z, method = "spearman", exact = FALSE)),
    .groups = "drop"
  ) %>%
  mutate(
    rho = sapply(test, function(x) x$estimate),
    p_value = sapply(test, function(x) x$p.value)
  ) %>%
  # Apply Benjamini-Hochberg correction across all p-values
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    label = paste0("ρ = ", round(rho, 2), "\n", 
                   "p ", format.pval(p_adj, digits = 2, eps = 0.001)),
    x = x_max - 2,
    y = y_min + 3
  )

# now plot:
f1z_a <- ggplot(all_F1s_3[all_F1s_3$regime=='asexual',], aes(x = YPD, y = pleio_z)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = 0, alpha = 0.2) +
  geom_vline(xintercept = 0, alpha = 0.2) +
  geom_point(fill='steelblue', shape = 21, size = 2, stroke = 1, col = "black") +
  geom_text(
    data = filter(rho_labels, regime == "asexual"),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    family = "sans",
    size = 3,
    hjust = 0
  ) +
  theme_gge() +
  facet_wrap(~kk_strain_id,ncol=3,nrow=2) +
  labs(y = "Mean Fitness in Away Environments (z-scored)", 
       x = "Fitness in Home Environment (z-scored)",
       title = "Asexual"
  ) +
  xlim(x_limits) + ylim(y_limits) +
  theme(
    legend.position = 'none',
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    strip.text = element_blank(),
    plot.title = element_text(size = 10, hjust = 0.01, vjust = -2, face = "bold")
  )

f1z_s <- ggplot(all_F1s_3[all_F1s_3$regime=='sexual',], aes(x = YPD, y = pleio_z)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = 0, alpha = 0.2) +
  geom_vline(xintercept = 0, alpha = 0.2) +
  geom_point(fill='darkorange2', shape = 21, size = 2, stroke = 1, col = "black") +
  geom_text(
    data = filter(rho_labels, regime == "sexual"),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    family = "sans",
    size = 3,
    hjust = 0
  ) +
  theme_gge() +
  facet_wrap(~kk_strain_id,ncol=3,nrow=2) +
  labs(y = "Mean Fitness in Away Environments (z-scored)", 
       x = "Fitness in Home Environment (z-scored)",
       title = "Sexual"
  ) +  
  xlim(x_limits) + ylim(y_limits) +
  theme(
    legend.position = 'none',
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    strip.text = element_blank(),
    plot.title = element_text(size = 10, hjust = 0.01, vjust = -2, face = "bold")
  )

home_away_corrs <- f1z_a + f1z_s + 
  plot_layout(ncol = 1, heights = c(1, 1), axis_title="collect")

home_away_corrs

ggsave(home_away_corrs, filename = 'figureS6.pdf', width = 8, height = 10, units = 'in', dpi = 300)
ggsave(home_away_corrs, filename = 'figureS6.jpg', width = 8, height = 10, units = 'in', dpi = 300)
