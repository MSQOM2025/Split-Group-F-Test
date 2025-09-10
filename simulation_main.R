# ==============================================================================
#
#    R Script for Simulation & Analysis for the F_ws 
#    "AThe Price of Purity: On the Sample-Size Dependency of a Split-Group Permutation Test for the Behrens-Fisher Problem"
#
# ==============================================================================

# ==============================================================================
#  1. SETUP: CHECK AND INSTALL PACKAGES
# ==============================================================================
required_packages <- c("tidyverse", "future", "furrr","qqplotr","ggplot2")
not_installed <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if(length(not_installed) > 0) {
  cat("Installing missing packages:", paste(not_installed, collapse=", "), "\n")
  install.packages(not_installed)
}
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(qqplotr))
suppressPackageStartupMessages(library(ggplot2))

select <- dplyr::select


# ==============================================================================
#  2. SIMULATION CONTROL PANEL
# ==============================================================================
set.seed(123)
TEST_MODE <- FALSE   # SET TO `TRUE` FOR A QUICK TEST RUN
if (TEST_MODE) { SIM_REPLICATIONS <- 100; PERMUTATIONS <- 49 } else { SIM_REPLICATIONS <- 2000; PERMUTATIONS <- 999 }


# ===================================================================================
#  3. DEFINE CORE CUSTOM FUNCTIONS
# ===================================================================================


generate_data <- function(k, n_per_group, means, sds) {
    observations <- vector("list", k); for (i in 1:k) { observations[[i]] <- rnorm(n_per_group[i], means[i], sds[i]) }
    tibble(observation = unlist(observations), group = factor(rep(1:k, times = n_per_group)))
}

calculate_fws_statistic <- function(data) {
  data <- data[!is.na(data$observation), ]; if (nrow(data) < 2) return(NA_real_)
  data$subgroup <- ave(1:nrow(data), data$group, FUN = function(x) sample(rep(1:2, length.out=length(x))))
  means <- aggregate(observation ~ group + subgroup, data=data, FUN=mean); vars  <- aggregate(observation ~ group + subgroup, data=data, FUN=var); ns    <- aggregate(observation ~ group + subgroup, data=data, FUN=length)
  if(any(is.na(vars$observation))) return(NA_real_)
  temp1 <- merge(means, vars, by = c("group", "subgroup")); subgroup_stats <- merge(temp1, ns, by = c("group", "subgroup")); names(subgroup_stats) <- c("group", "subgroup", "mean", "var", "n")
  if (any(subgroup_stats$n <= 1)) return(NA_real_)
  wide_stats_list <- split(subgroup_stats, subgroup_stats$group)
  process_group <- function(df) {
    if (nrow(df) < 2) return(NULL); data.frame(group=df$group[1], n_1=df$n[df$subgroup==1], var_1=df$var[df$subgroup==1], mean_1=df$mean[df$subgroup==1], n_2=df$n[df$subgroup==2], var_2=df$var[df$subgroup==2], mean_2=df$mean[df$subgroup==2])
  }
  wide_stats <- do.call(rbind, lapply(wide_stats_list, process_group)); if (is.null(wide_stats) || nrow(wide_stats) == 0) return(NA_real_)
  w_i <- (1 / wide_stats$var_1) + (1 / wide_stats$var_2)
  x_i_star <- ((wide_stats$mean_1 / wide_stats$var_1) + (wide_stats$mean_2 / wide_stats$var_2)) / w_i
  w_total <- sum(w_i, na.rm=TRUE); if (w_total == 0) return(NA_real_)
  x_w_star <- sum(w_i * x_i_star, na.rm=TRUE) / w_total
  sum(w_i * (x_i_star - x_w_star)^2, na.rm=TRUE)
}


perm_fws_test <- function(data) {
  t_obs <- calculate_fws_statistic(data); if (is.na(t_obs)) return(list(p.value=NA_real_, statistic=NA_real_))
  original_labels <- data$group; perm_data <- data; t_perm <- numeric(PERMUTATIONS)
  for (i in 1:PERMUTATIONS) { perm_data$group <- sample(original_labels); t_perm[i] <- calculate_fws_statistic(perm_data) }
  p_value <- (sum(t_perm >= t_obs, na.rm = TRUE) + 1) / (PERMUTATIONS + 1)
  return(list(p.value=p_value, statistic=t_obs))
}


# ==========================================================
# run_scenario function 
# ==========================================================


run_scenario <- function(params) {
  cat("Running scenario:", paste(names(params), params, sep="=", collapse=", "), "\n")

  if (params$effect_type == "null") { 
    means <- rep(0, 4)
    sds <- params$sd_structure[[1]] 
  } else { 
    means <- c(0, 0, 0.5, 0.8)
    sds <- params$sd_structure[[1]] 
  }
  
  # Pre-allocate result vectors
  results_classic <- numeric(SIM_REPLICATIONS)
  results_welch <- numeric(SIM_REPLICATIONS)
  results_fws <- numeric(SIM_REPLICATIONS)

  for (i in 1:SIM_REPLICATIONS) {
    sim_data <- generate_data(k=4, n_per_group=rep(params$n, 4), means=means, sds=sds)

    # --- Classic ANOVA ---
    p_classic <- try(summary(aov(observation ~ group, data=sim_data))[[1]][["Pr(>F)"]][1], silent=TRUE)
    results_classic[i] <- if(inherits(p_classic, "try-error") || length(p_classic)==0) NA else p_classic
    
    # --- Welch's ANOVA ---
    p_welch <- try(oneway.test(observation ~ group, data=sim_data)$p.value, silent=TRUE)
    results_welch[i] <- if(inherits(p_welch, "try-error") || length(p_welch)==0) NA else p_welch
    
    # --- Perm F_ws Test ---
    p_fws_result <- try(perm_fws_test(sim_data), silent=TRUE)
    p_fws <- if (is.list(p_fws_result)) p_fws_result$p.value else NA
    results_fws[i] <- if(inherits(p_fws_result, "try-error") || is.null(p_fws) || length(p_fws)==0) NA else p_fws
  }

  # --- Combine results ---
  power_results <- tibble(
    test = c("classic_f", "welch_f", "perm_fws"),
    value = c(mean(results_classic < 0.05, na.rm = TRUE),
              mean(results_welch < 0.05, na.rm = TRUE),
              mean(results_fws < 0.05, na.rm = TRUE))
  )

  # Attach parameters and return
  as_tibble(params) %>%
    mutate(sd_structure = paste(params$sd_structure[[1]], collapse = ":")) %>%
    bind_cols(power_results) %>%
    mutate(
      result_metric = if_else(params$effect_type == "null", "type1_error", "power"),
      n_sims = SIM_REPLICATIONS
    )
}

# 4. SETUP AND RUN STUDY
# ==========================
num_cores <- 15
plan(multisession, workers = num_cores)

simulation_grid <- expand_grid(n=c(20,50), sd_structure=list(c(1,1,1,1), c(1,2,4,8)), effect_type=c("null","mean_shift"))

cat("Starting Main F_ws simulation study...\n")
start_time <- Sys.time()
main_sim_results <- future_map(
  .x = 1:nrow(simulation_grid),
  ~ run_scenario(params = as.list(simulation_grid[.x, ])),
  .progress = TRUE,
  .options = furrr_options(seed=TRUE, packages=c("tidyverse"), globals=c("SIM_REPLICATIONS", "PERMUTATIONS", "simulation_grid", "run_scenario", "generate_data", "perm_fws_test", "calculate_fws_statistic"))
)
end_time <- Sys.time()

# 5. P-VALUE ACCURACY 
start_time_pvalue <- Sys.time()

cat("\nStarting P-Value Accuracy deep-dive simulation...\n")
null_params_for_accuracy <- list(n=20, sd_structure=list(c(1,2,4,8)), effect_type="null")

# Pre-allocate result vectors
pvals_classic <- numeric(SIM_REPLICATIONS)
pvals_welch <- numeric(SIM_REPLICATIONS)
pvals_fws <- numeric(SIM_REPLICATIONS)

for (i in 1:SIM_REPLICATIONS) {
    sim_data <- generate_data(k=4, n_per_group=rep(null_params_for_accuracy$n, 4), means=rep(0,4), sds=null_params_for_accuracy$sd_structure[[1]])

    p_classic <- try(summary(aov(observation ~ group, data=sim_data))[[1]][["Pr(>F)"]][1], silent=TRUE)
    pvals_classic[i] <- if(inherits(p_classic, "try-error") || length(p_classic)==0) NA else p_classic
    
    p_welch <- try(oneway.test(observation ~ group, data=sim_data)$p.value, silent=TRUE)
    pvals_welch[i] <- if(inherits(p_welch, "try-error") || length(p_welch)==0) NA else p_welch
    
    p_fws_result <- try(perm_fws_test(sim_data), silent=TRUE)
    p_fws <- if(is.list(p_fws_result)) p_fws_result$p.value else NA
    pvals_fws[i] <- if(inherits(p_fws_result, "try-error") || is.null(p_fws) || length(p_fws)==0) NA else p_fws
}

# Combine the vectors into a matrix for easy processing
all_p_values_matrix <- rbind(
  classic_f = pvals_classic,
  welch_f = pvals_welch,
  perm_fws = pvals_fws
)

cat("P-Value Accuracy simulation finished.\n")
end_time_pvalue <- Sys.time()

# 6. SAVE, PROCESS, AND PLOT
final_results_table <- bind_rows(main_sim_results)
write_csv(final_results_table, "simulation_results_fws_final.csv")
print("--- Final Results Table (Main Study) ---")
print(final_results_table)
cat("\nFinal main results saved to simulation_results_fws_final.csv\n")

# Prepare data for plotting
plot_data <- final_results_table %>%
  mutate(test_label = factor(case_when(
    test == "classic_f" ~ "Classic ANOVA",
    test == "welch_f"   ~ "Welch's ANOVA",
    test == "perm_fws"  ~ "Permutation F_ws (Ours)"
  ), levels = c("Classic ANOVA", "Welch's ANOVA", "Permutation F_ws (Ours)"))) %>%
  mutate(sd_label = if_else(sd_structure == "1:1:1:1", "Homoscedastic", "Heteroscedastic"))

type1_error_plot_data <- plot_data %>% filter(effect_type == "null")
plot1 <- ggplot(type1_error_plot_data, aes(x = test_label, y = value, fill = test_label)) +
  geom_bar(stat = "identity", color = "black") +
  geom_hline(yintercept=0.05, linetype="dashed", color="red", linewidth=1) +
  geom_text(aes(label = sprintf("%.3f", value)), vjust=-0.5) +
  facet_grid(sd_label ~ paste("n =", n)) +
  labs(title="Figure 1: Type I Error Control", x="Test", y="Empirical Type I Error") +
  scale_fill_manual(values=c("Classic ANOVA"="gray70", "Welch's ANOVA"="skyblue3", "Permutation F_ws (Ours)"="darkorange")) +
  theme_bw(base_size = 14) + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none") +
  coord_cartesian(ylim=c(0,0.20))
ggsave("plot_type1_error.png", plot=plot1, width=10, height=8, dpi=300)
cat("Saved plot_type1_error.png\n")

# --- PLOT 2: Statistical Power ---
plot2_data <- plot_data %>%
  filter(effect_type == "mean_shift",
         !(sd_label == "Heteroscedastic" & test_label %in% c("Classic ANOVA", "Permutation F_ws (Ours)")))
plot2 <- ggplot(plot2_data, aes(x=factor(n), y=value, group=test_label, color=test_label)) +
  geom_line(linewidth=1.2) + geom_point(size=4, aes(shape=test_label)) +
  facet_wrap(~sd_label) +
  labs(title="Figure 2: Statistical Power Comparison", x="Sample Size (n)", y="Power", color="Test", shape="Test") +
  scale_color_manual(values=c("Classic ANOVA"="gray70", "Welch's ANOVA"="skyblue3", "Permutation F_ws (Ours)"="darkorange"), drop=FALSE) +
  scale_shape_manual(values=c("Classic ANOVA"=15, "Welch's ANOVA"=17, "Permutation F_ws (Ours)"=19), drop=FALSE) +
  theme_bw(base_size = 14) + theme(legend.position="bottom") + coord_cartesian(ylim=c(0,1))
ggsave("plot_power_comparison.png", plot=plot2, width=12, height=7, dpi=300)
cat("Saved plot_power_comparison.png\n")

# --- PLOT 3: P-Value Accuracy ---
p_value_dist_data <- as_tibble(t(all_p_values_matrix)) %>%
  pivot_longer(everything(), names_to="test", values_to="p_value") %>%
  mutate(test_label = factor(case_when(
    test == "classic_f" ~ "Classic ANOVA",
    test == "welch_f"   ~ "Welch's ANOVA",
    test == "perm_fws"  ~ "Permutation F_ws (Ours)"
  ), levels=c("Classic ANOVA", "Welch's ANOVA", "Permutation F_ws (Ours)"))) %>%
  filter(!is.na(p_value))

plot3 <- ggplot(p_value_dist_data, aes(sample = p_value)) +
  ggplot2::stat_qq(distribution = stats::qunif, alpha = 0.6) +
  ggplot2::stat_qq_line(distribution = stats::qunif, color = "red", linetype = "dashed", linewidth=1) +
  facet_wrap(~test_label) +
  labs(title = "Figure 3: P-Value Accuracy under Heteroscedastic Null (n=20)",
       subtitle = "Accurate p-values should follow the red dashed line (Uniform Distribution)",
       x = "Theoretical Quantiles (Uniform)",
       y = "Sample Quantiles (Observed P-Values)") +
  theme_bw(base_size = 14)

ggsave("plot_pvalue_accuracy_qq.png", plot=plot3, width=10, height=8, dpi=300)
cat("Saved plot_pvalue_accuracy_qq.png\n")


cat("\nAll plots have been generated.\n")


cat("Main simulation study finished. Total time:", format(end_time - start_time), "\n")
cat("P-Value Accuracy deep-dive finished. Total time:", format(end_time_pvalue - start_time_pvalue), "\n")
