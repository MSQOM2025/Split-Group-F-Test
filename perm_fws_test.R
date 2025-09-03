# ==============================================================================
#
#    Permutation Split-Group F-Test (Perm F_ws)
#    A standalone R function for the method described in:
#    M.S (2025). "A Split-Group Permutation F-Test for Robust ANOVA under Heteroscedasticity"
#
# ==============================================================================

# ===================================================================================
#  CORE STATISTIC FUNCTION
#  This version uses Base R for memory stability.
# ===================================================================================

# An internal helper function to compute the bias-corrected, weighted
# between-group sum of squares statistic.
# INPUT: A data.frame or tibble with two columns: 'observation' and 'group'.
# OUTPUT: The numeric value of the F_ws statistic.
calculate_fws_statistic <- function(data) {
  data$subgroup <- ave(1:nrow(data), data$group, FUN = function(x) sample(rep(1:2, length.out = length(x))))
  
  means <- aggregate(observation ~ group + subgroup, data = data, FUN = mean)
  vars <- aggregate(observation ~ group + subgroup, data = data, FUN = var)
  ns <- aggregate(observation ~ group + subgroup, data = data, FUN = length)
  
  temp1 <- merge(means, vars, by = c("group", "subgroup"))
  subgroup_stats <- merge(temp1, ns, by = c("group", "subgroup"))
  names(subgroup_stats) <- c("group", "subgroup", "mean", "var", "n")
  
  if (any(subgroup_stats$n <= 3)) {
    return(NA_real_)
  }
  
  wide_stats_list <- split(subgroup_stats, subgroup_stats$group)
  process_group <- function(df) {
    if(nrow(df) < 2) return(NULL)
    data.frame(
      group = df$group[1],
      n_1 = df$n[df$subgroup == 1], var_1 = df$var[df$subgroup == 1], mean_1 = df$mean[df$subgroup == 1],
      n_2 = df$n[df$subgroup == 2], var_2 = df$var[df$subgroup == 2], mean_2 = df$mean[df$subgroup == 2]
    )
  }
  wide_stats <- do.call(rbind, lapply(wide_stats_list, process_group))

  # Return NA if splitting resulted in no valid groups
  if (is.null(wide_stats) || nrow(wide_stats) == 0) return(NA_real_)

  w_i <- ((wide_stats$n_1 - 3) / ((wide_stats$n_1 - 1) * wide_stats$var_1)) + 
         ((wide_stats$n_2 - 3) / ((wide_stats$n_2 - 1) * wide_stats$var_2))
  
  x_i_star <- (((wide_stats$mean_1 * (wide_stats$n_1 - 3) / ((wide_stats$n_1 - 1) * wide_stats$var_1)) + 
                (wide_stats$mean_2 * (wide_stats$n_2 - 3) / ((wide_stats$n_2 - 1) * wide_stats$var_2))) / w_i)
  
  # --- THE ROBUST VERSION OF THE FINAL CALCULATIONS ---
  w_total <- sum(w_i, na.rm = TRUE)
  # Prevent division by zero if all weights were NA
  if (w_total == 0) return(NA_real_)
  x_w_star <- sum(w_i * x_i_star, na.rm = TRUE) / w_total
  
  numerator <- sum(w_i * (x_i_star - x_w_star)^2, na.rm = TRUE)

  return(numerator)
}```

# ==============================================================================
#  MAIN PUBLIC FUNCTION
# ==============================================================================

# Performs a robust, permutation-based F-test for comparing k >= 2 group
# means under variance heterogeneity (the Behrens-Fisher problem).
perm_fws_test <- function(data, perms = 999) {
  
  # --- Input validation ---
  if (!is.data.frame(data) || ncol(data) != 2) {
    stop("Input 'data' must be a data.frame or tibble with exactly two columns.")
  }
  names(data) <- c("observation", "group")
  if (!is.numeric(data$observation)) {
    stop("The first column ('observation') must be numeric.")
  }
  if (!is.factor(data$group)) {
    data$group <- as.factor(data$group)
    warning("The group column was not a factor and has been converted.")
  }

  # --- Test Execution ---
  t_obs <- calculate_fws_statistic(data)
  if (is.na(t_obs)) {
    warning("Could not calculate observed statistic, possibly due to small subgroup size.")
    return(list(p.value = NA_real_, statistic = NA_real_, permutations = perms))
  }

  # Memory-efficient permutation loop
  original_labels <- data$group
  perm_data <- data 
  t_perm <- numeric(perms)

  for (i in 1:perms) {
    perm_data$group <- sample(original_labels)
    t_perm[i] <- calculate_fws_statistic(perm_data)
  }
  
  p_value <- (sum(t_perm >= t_obs, na.rm = TRUE) + 1) / (perms + 1)
  
  # --- Return a structured result ---
  result <- list(
    p.value = p_value,
    statistic = t_obs,
    permutations = perms
  )
  
  return(result)
}

# ==============================================================================
#  EXAMPLE USAGE
# ==============================================================================

# This code will only run if the script is sourced directly, not when the
# function is just loaded.
if (sys.nframe() == 0) {

  cat("--- Running example on the CO2 dataset ---\n")
  
  # Load the CO2 dataset from base R
  data(CO2)
  
  # Create the interaction groups (k=4, n=21 per group)
  # This dataset is known to be heteroscedastic.
  co2_data <- data.frame(
    observation = CO2$uptake,
    group = interaction(CO2$Type, CO2$Treatment, sep = "-")
  )
  
  cat("\n--- Bartlett's test for Heteroscedasticity ---\n")
  print(bartlett.test(observation ~ group, data = co2_data))
  
  cat("\n--- Running Permutation F_ws Test ---\n")
  
  # Set a seed for a reproducible example
  set.seed(123)
  
  # Run the test
  test_result <- perm_fws_test(co2_data, perms = 999)
  
  # Print the results
  cat("Observed Statistic:", test_result$statistic, "\n")
  cat("P-value:", test_result$p.value, "\n")
  

}
