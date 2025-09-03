# ==============================================================================
#
#    R Script for Simulation & Analysis for the F_ws 
#    "A Split-Group Permutation F-Test for Robust ANOVA under Heteroscedasticity"
#    FINAL VERSION: Includes Main Study + P-Value Accuracy Deep Dive
#
# ==============================================================================

# ==============================================================================
#  1. SETUP: CHECK AND INSTALL PACKAGES
# ==============================================================================

required_packages <- c("tidyverse", "future", "furrr")

not_installed <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if(length(not_installed) > 0) {
  cat("Installing missing packages:", paste(not_installed, collapse=", "), "\n")
  install.packages(not_installed)
}

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(furrr))

select <- dplyr::select
cat("Setup complete. Libraries loaded.\n")

# ==============================================================================
#  2. SIMULATION CONTROL PANEL
# ==============================================================================
set.seed(1951)
TEST_MODE <- FALSE # <-- SET TO `TRUE` for a quick test run

if (TEST_MODE) {
  SIM_REPLICATIONS <- 100
  PERMUTATIONS <- 49
  cat("--- RUNNING IN TEST MODE ---\n")
} else {
  SIM_REPLICATIONS <- 2000
  PERMUTATIONS <- 999
  cat("--- RUNNING IN PUBLICATION MODE ---\n")
}

# =====================================
#  3. DEFINE CORE CUSTOM FUNCTIONS
# =====================================

generate_data <- function(k, n_per_group, means, sds) {
  observations <- vector("list", k)
  for (i in 1:k) {
    observations[[i]] <- rnorm(n = n_per_group[i], mean = means[i], sd = sds[i])
  }
  tibble(observation = unlist(observations), group = factor(rep(1:k, times = n_per_group)))
}

# ===================================================================================
#  #  This  uses Base R for memory stability AND fixes the column naming bug.
# ===================================================================================

calculate_fws_statistic <- function(data) {
    
  # Add the subgroup column using base R
  data$subgroup <- ave(1:nrow(data), data$group, FUN = function(x) sample(rep(1:2, length.out = length(x))))
  
  # Calculate subgroup stats using aggregate
  means <- aggregate(observation ~ group + subgroup, data = data, FUN = mean)
  vars <- aggregate(observation ~ group + subgroup, data = data, FUN = var)
  ns <- aggregate(observation ~ group + subgroup, data = data, FUN = length)
  
  # --- THE FIX IS HERE ---
  # Combine into a single data frame by merging correctly by group and subgroup
  # and then rename the data columns properly.
  temp1 <- merge(means, vars, by = c("group", "subgroup"))
  subgroup_stats <- merge(temp1, ns, by = c("group", "subgroup"))
  # The aggregated columns will be named observation, observation.x, observation.y
  names(subgroup_stats) <- c("group", "subgroup", "mean", "var", "n")
  
  # Check the subgroup size constraint
  if (any(subgroup_stats$n <= 3)) {
    return(NA_real_)
  }
  
  # Reshape to wide format using base R
  wide_stats_list <- split(subgroup_stats, subgroup_stats$group)
  process_group <- function(df) {
    data.frame(
      group = df$group[1],
      n_1 = df$n[df$subgroup == 1], var_1 = df$var[df$subgroup == 1], mean_1 = df$mean[df$subgroup == 1],
      n_2 = df$n[df$subgroup == 2], var_2 = df$var[df$subgroup == 2], mean_2 = df$mean[df$subgroup == 2]
    )
  }
  wide_stats <- do.call(rbind, lapply(wide_stats_list, process_group))

  # Calculate F_ws components using base R vectorized operations
  w_i <- ((wide_stats$n_1 - 3) / ((wide_stats$n_1 - 1) * wide_stats$var_1)) + 
         ((wide_stats$n_2 - 3) / ((wide_stats$n_2 - 1) * wide_stats$var_2))
  
  x_i_star <- (((wide_stats$mean_1 * (wide_stats$n_1 - 3) / ((wide_stats$n_1 - 1) * wide_stats$var_1)) + 
                (wide_stats$mean_2 * (wide_stats$n_2 - 3) / ((wide_stats$n_2 - 1) * wide_stats$var_2))) / w_i)
  
  w_total <- sum(w_i,na.rm=TRUE)
  x_w_star <- sum(w_i * x_i_star,na.rm=TRUE) / w_total
  
  numerator <- sum(w_i * (x_i_star - x_w_star)^2)

  return(numerator)
}
# =========================================================
#   MEMORY-EFFICIENT PERMUTATION FUNCTION
# =========================================================

perm_fws_test <- function(data) {
  # Calculate the observed statistic only once
  t_obs <- calculate_fws_statistic(data)
  if (is.na(t_obs)) return(NA_real_)

  # Store the original labels and create a copy of the data that we can modify
  original_labels <- data$group
  perm_data <- data # This is the ONLY full copy we make.

  # Pre-allocate the vector to store results.
  t_perm <- numeric(PERMUTATIONS)

  for (i in 1:PERMUTATIONS) {
    # The EFFICIENT part: Shuffle ONLY the labels vector.
    # Do NOT create a new data.frame.
    perm_data$group <- sample(original_labels)
    
    # Recalculate the statistic on the modified data.
    t_perm[i] <- calculate_fws_statistic(perm_data)
  }
  
  p_value <- (sum(t_perm >= t_obs, na.rm = TRUE) + 1) / (PERMUTATIONS + 1)
  return(p_value)
}

run_scenario <- function(params) {
  cat("Running scenario:", paste(names(params), params, sep="=", collapse=", "), "\n")
  
  if (params$effect_type == "null") {
    means <- rep(0, 4)
    sds <- params$sd_structure[[1]]
  } else if (params$effect_type == "mean_shift") {
    means <- c(0, 0, 0.5, 0.8)
    sds <- params$sd_structure[[1]]
  }
  
  p_values_list <- replicate(SIM_REPLICATIONS, {
    sim_data <- generate_data(k = 4, n_per_group = rep(params$n, 4), means = means, sds = sds)
    c(classic_f = summary(aov(observation ~ group, data = sim_data))[[1]][["Pr(>F)"]][1],
      welch_f = oneway.test(observation ~ group, data = sim_data, var.equal = FALSE)$p.value,
      perm_fws = perm_fws_test(sim_data))
  }, simplify = "matrix")

  power_results <- tibble(
    test = rownames(p_values_list),
    value = rowMeans(p_values_list < 0.05, na.rm = TRUE)
  )
  
  as_tibble(params) %>%
    mutate(sd_structure = paste(params$sd_structure[[1]], collapse = ":")) %>%
    bind_cols(power_results) %>%
    mutate(result_metric = if_else(params$effect_type == "null", "type1_error", "power"),
           n_sims = SIM_REPLICATIONS)
}

# ==========================
#  4. SETUP AND RUN STUDY
# ==========================
# Let's manually set a lower number of cores to avoid maxing out RAM.
# With 32 GB RAM (16 GB available), 3 or 4 is a much safer number.
num_cores <- 2
plan(multisession, workers = num_cores)

cat(paste("Parallel processing enabled on", num_cores, "cores.\n"))

simulation_grid <- expand_grid(
  n = c(20, 50),
  sd_structure = list(c(1, 1, 1, 1), c(1, 2, 4, 8)),
  effect_type = c("null", "mean_shift")
)

cat("Starting Main F_ws simulation study...\n")
start_time <- Sys.time()
main_sim_results <- future_map(
  .x = 1:nrow(simulation_grid),
  ~ run_scenario(params = as.list(simulation_grid[.x, ])),
  .progress = TRUE,
  .options = furrr_options(seed = TRUE, packages = "tidyverse")
)
end_time <- Sys.time()
cat("Main simulation study finished. Total time:", format(end_time - start_time), "\n")

# =================================
#  5. P-VALUE ACCURACY DEEP DIVE
# =================================

cat("\nStarting P-Value Accuracy deep-dive simulation...\n")

# Define the single, challenging null scenario to investigate
null_params_for_accuracy <- list(
  n = 20,
  sd_structure = list(c(1, 2, 4, 8)),
  effect_type = "null"
)

# This block is also parallelized by the plan() we set up earlier
all_p_values_matrix <- replicate(SIM_REPLICATIONS, {
    sim_data <- generate_data(
      k = 4,
      n_per_group = rep(null_params_for_accuracy$n, 4),
      means = rep(0, 4),
      sds = null_params_for_accuracy$sd_structure[[1]]
    )
    c(
      classic_f = summary(aov(observation ~ group, data = sim_data))[[1]][["Pr(>F)"]][1],
      welch_f = oneway.test(observation ~ group, data = sim_data, var.equal = FALSE)$p.value,
      perm_fws = perm_fws_test(sim_data)
    )
}, simplify = "matrix")

cat("P-Value Accuracy simulation finished.\n")

# =================================
#  6. SAVE, PROCESS, AND PLOT
# =================================

# --- Process and Save Main Results ---
final_results_table <- bind_rows(main_sim_results)
print("--- Final Results Table (Main Study) ---")
print(final_results_table)
write_csv(final_results_table, "simulation_results_fws_final.csv")
cat("\nFinal main results saved to simulation_results_fws_final.csv\n")

# --- Prepare data for plotting ---
plot_data <- final_results_table %>%
  mutate(test_label = factor(case_when(
    test == "classic_f" ~ "Classic ANOVA",
    test == "welch_f"   ~ "Welch's ANOVA",
    test == "perm_fws"  ~ "Permutation F_ws (Ours)"
  ), levels = c("Classic ANOVA", "Welch's ANOVA", "Permutation F_ws (Ours)"))) %>%
  mutate(sd_label = if_else(sd_structure == "1:1:1:1", "Homoscedastic", "Heteroscedastic"))

# --- PLOT 1: Type I Error ---
type1_error_plot_data <- plot_data %>% filter(effect_type == "null")
plot1 <- ggplot(type1_error_plot_data, aes(x = test_label, y = value, fill = test_label)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) +
  geom_text(aes(label = sprintf("%.3f", value)), vjust = -0.5, size = 3.5) +
  facet_grid(sd_label ~ paste("n =", n)) +
  labs(title = "Figure 1: Type I Error Control", x = "Test Type", y = "Empirical Type I Error Rate") +
  scale_fill_manual(values = c("gray70", "skyblue3", "darkorange")) +
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.15))

print(plot1)
ggsave("plot_type1_error.png", plot = plot1, width = 10, height = 8, dpi = 300)

# --- PLOT 2: Statistical Power ---
power_plot_data <- plot_data %>%
  filter(effect_type == "mean_shift", !(sd_label == "Heteroscedastic" & test_label == "Classic ANOVA"))
plot2 <- ggplot(power_plot_data, aes(x = factor(n), y = value, group = test_label, color = test_label)) +
  geom_line(size = 1.2) + geom_point(size = 4, aes(shape = test_label)) +
  geom_text(aes(label = sprintf("%.2f", value)), vjust = -1.5, size = 4) +
  facet_wrap(~sd_label) +
  labs(title = "Figure 2: Statistical Power Comparison", x = "Sample Size per Group (n)", y = "Empirical Power", color = "Test", shape = "Test") +
  scale_color_manual(values = c("gray70", "skyblue3", "darkorange")) +
  scale_shape_manual(values = c(15, 17, 19)) +
  theme_bw(base_size = 14) + theme(legend.position = "bottom") + coord_cartesian(ylim = c(0, 1.0))

print(plot2)
ggsave("plot_power_comparison.png", plot = plot2, width = 12, height = 7, dpi = 300)

# --- PLOT 3: P-Value Accuracy ---
p_value_dist_data <- as_tibble(t(all_p_values_matrix)) %>%
  pivot_longer(everything(), names_to = "test", values_to = "p_value") %>%
  mutate(test_label = factor(case_when(
    test == "classic_f" ~ "Classic ANOVA",
    test == "welch_f"   ~ "Welch's ANOVA",
    test == "perm_fws"  ~ "Permutation F_ws (Ours)"
  ), levels = c("Classic ANOVA", "Welch's ANOVA", "Permutation F_ws (Ours)")))

plot3 <- ggplot(p_value_dist_data, aes(sample = p_value)) +
  stat_qq(distribution = qunif, size = 1, alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size=1) +
  facet_wrap(~test_label) +
  labs(title = "Figure 3: P-Value Accuracy under Heteroscedastic Null (n=20)",
       subtitle = "Accurate p-values should follow the red dashed line (Uniform Distribution)",
       x = "Theoretical Quantiles (Uniform)", y = "Sample Quantiles (Observed P-Values)") +
  theme_bw(base_size = 14)

print(plot3)

ggsave("plot_pvalue_accuracy_qq.png", plot = plot3, width = 12, height = 5, dpi = 300)
