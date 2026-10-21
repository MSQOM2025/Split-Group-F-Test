# ==============================================================================
# Script: perm_fws_engine.R
# Purpose: High-performance, memory-efficient implementation of the 
#          Split-Group Permutation Test (Perm F_ws) in pure Base R.
# ==============================================================================

#' Calculate the Split-Group Dispersion Statistic (T)
#'
#' @param x A numeric vector of observations.
#' @param group_sizes An integer vector containing sample sizes (n_1, ..., n_k).
#' @param n1_vec An integer vector containing subgroup 1 sizes (n_{11}, ..., n_{k1}).
#' @param n2_vec An integer vector containing subgroup 2 sizes (n_{12}, ..., n_{k2}).
#' @param split_idx A list of index pairs for slicing subgroups.
#' @return The scalar test statistic T.
calc_fws_stat_fast <- function(x, split_idx, k) {
  x_bar_star <- numeric(k)
  w <- numeric(k)
  
  for (i in 1:k) {
    # Extract observations for each subgroup using pre-calculated index slices
    sub1 <- x[split_idx[[i]]$idx1]
    sub2 <- x[split_idx[[i]]$idx2]
    
    n1 <- split_idx[[i]]$n1
    n2 <- split_idx[[i]]$n2
    
    # Subgroup means
    m1 <- sum(sub1) / n1
    m2 <- sum(sub2) / n2
    
    # Subgroup unbiased sample variances (Sum of squared deviations / (n - 1))
    v1 <- sum((sub1 - m1)^2) / (n1 - 1)
    v2 <- sum((sub2 - m2)^2) / (n2 - 1)
    
    # Split-group point estimator of group mean: Eq. (4)
    x_bar_star[i] <- 0.5 * (m1 + m2)
    
    # Unbiased plug-in estimator of Var(x_bar_star): Eq. (7)
    # V_i^* = 0.25 * (S_{i1}^2 / n_{i1} + S_{i2}^2 / n_{i2})
    v_star <- 0.25 * (v1 / n1 + v2 / n2)
    
    # Guard against degenerate variance estimates in extreme edge-cases
    if (v_star <= .Machine$double.eps) {
      return(NA_real_)
    }
    
    # Group precision weight: Eq. (8)
    w[i] <- 1 / v_star
  }
  
  # Weighted grand mean: \bar{X}_w^*
  sum_w <- sum(w)
  x_bar_w <- sum(w * x_bar_star) / sum_w
  
  # Final quadratic form statistic: Eq. (9)
  # T = sum( w_i * (x_bar_star_i - x_bar_w)^2 )
  t_stat <- sum(w * (x_bar_star - x_bar_w)^2)
  return(t_stat)
}


#' Split-Group Permutation Test for One-Way Heteroscedastic ANOVA
#'
#' @param x A numeric vector of observations of length N.
#' @param group A factor or vector indicating group membership (length N).
#' @param B Number of random permutations (default = 999).
#' @return A list containing the observed statistic, p-value, and permuted null distribution.
perm_fws_test <- function(x, group, B = 999) {
  if (length(x) != length(group)) {
    stop("Lengths of observations vector 'x' and 'group' vector must match.")
  }
  
  # Ensure group is factor and identify structure
  group <- as.factor(group)
  group_levels <- levels(group)
  k <- length(group_levels)
  
  if (k < 2) {
    stop("At least two groups are required.")
  }
  
  # Calculate group sizes
  group_sizes <- as.integer(table(group))
  n1_vec <- group_sizes %/% 2
  n2_vec <- group_sizes - n1_vec
  
  if (any(n1_vec <= 1 | n2_vec <= 1)) {
    stop("Sample size per group too small; each split must contain at least 2 observations.")
  }
  
  # Pre-calculate contiguous slice indices for maximum performance
  split_idx <- vector("list", k)
  current_pos <- 1L
  for (i in 1:k) {
    n1 <- n1_vec[i]
    n2 <- n2_vec[i]
    idx1 <- current_pos:(current_pos + n1 - 1L)
    idx2 <- (current_pos + n1):(current_pos + n1 + n2 - 1L)
    split_idx[[i]] <- list(idx1 = idx1, idx2 = idx2, n1 = n1, n2 = n2)
    current_pos <- current_pos + n1 + n2
  }
  
  # Sort input vector according to group order for the observed statistic
  order_idx <- order(group)
  x_sorted <- x[order_idx]
  
  # Calculate observed test statistic
  t_obs <- calc_fws_stat_fast(x_sorted, split_idx, k)
  if (is.na(t_obs)) {
    return(list(statistic = NA_real_, p_value = NA_real_, B = B))
  }
  
  # Permutation loop: Shuffling pooled vector directly
  n_total <- length(x)
  t_perm <- numeric(B)
  
  for (b in 1:B) {
    # Generate random permutation of observations
    x_perm <- x[sample.int(n_total)]
    t_perm[b] <- calc_fws_stat_fast(x_perm, split_idx, k)
  }
  
  # Drop degenerate permutations if any occurred
  valid_perm <- !is.na(t_perm)
  n_valid <- sum(valid_perm)
  
  if (n_valid < (0.5 * B)) {
    warning("More than 50% of permutations produced invalid variance estimates.")
    return(list(statistic = t_obs, p_value = NA_real_, B = B))
  }
  
  # Calculate permutation p-value: Eq. (10)
  p_value <- (sum(t_perm[valid_perm] >= t_obs) + 1) / (n_valid + 1)
  
  return(list(
    statistic = t_obs,
    p_value = p_value,
    perm_distribution = t_perm[valid_perm],
    df = k - 1,
    B = n_valid
  ))
}
# ==============================================================================
# Script: simulation_table1.R
# Purpose: Reproduce the COMPLETE Table 1 of the manuscript (M = 2000 replications)
#          including both Type I Error (Null) and Empirical Power (Alternative).
# ==============================================================================

# 1. Source the core testing engine
source("perm_fws.R")

set.seed(110)

M_REPS <- 2000     # Number of Monte Carlo replications
B_PERM <- 999      # Number of permutations per test
ALPHA  <- 0.05

# ------------------------------------------------------------------------------
# Core Simulation Engine Function
# ------------------------------------------------------------------------------
run_monte_carlo_scenario <- function(n, sds, effect_delta = 0) {
  k <- 4
  group <- rep(1:k, each = n)
  means <- c(0, 0, 0.5 * effect_delta, 0.8 * effect_delta)
  
  classic_p <- numeric(M_REPS)
  welch_p   <- numeric(M_REPS)
  perm_p    <- numeric(M_REPS)
  
  cat(sprintf("\nRunning: n = %d | sds = (%s) | Delta = %.1f ...\n", 
              n, paste(sds, collapse=","), effect_delta))
  
  pb <- txtProgressBar(min = 0, max = M_REPS, style = 3)
  
  for (m in 1:M_REPS) {
    # Generate data
    x <- unlist(lapply(1:k, function(i) rnorm(n, mean = means[i], sd = sds[i])))
    
    # 1. Classic ANOVA
    classic_p[m] <- summary(aov(x ~ factor(group)))[[1]][["Pr(>F)"]][1]
    
    # 2. Welch ANOVA
    welch_p[m] <- oneway.test(x ~ group, var.equal = FALSE)$p.value
    
    # 3. Split-Group Permutation Test (Our method)
    perm_p[m] <- perm_fws_test(x, group, B = B_PERM)$p_value
    
    if (m %% 100 == 0) setTxtProgressBar(pb, m)
  }
  close(pb)
  
  # Calculate Rejection Rates
  res <- c(
    Classic = mean(classic_p < ALPHA, na.rm = TRUE),
    Welch   = mean(welch_p < ALPHA, na.rm = TRUE),
    PermFws = mean(perm_p < ALPHA, na.rm = TRUE)
  )
  return(res)
}

# ==============================================================================
# PART 1: TYPE I ERROR SCENARIOS (Delta = 0.0)
# ==============================================================================
cat("\n=======================================================\n")
cat("      PART 1: EVALUATING EMPIRICAL TYPE I ERROR        \n")
cat("=======================================================\n")

# Scenario 1: n = 20, Homoscedastic Null
res_n20_homo_null <- run_monte_carlo_scenario(n = 20, sds = c(1, 1, 1, 1), effect_delta = 0)

# Scenario 2: n = 20, Heteroscedastic Null (Severe)
res_n20_hetero_null <- run_monte_carlo_scenario(n = 20, sds = c(1, 2, 4, 8), effect_delta = 0)

# Scenario 3: n = 50, Homoscedastic Null
res_n50_homo_null <- run_monte_carlo_scenario(n = 50, sds = c(1, 1, 1, 1), effect_delta = 0)

# Scenario 4: n = 50, Heteroscedastic Null (Severe)
res_n50_hetero_null <- run_monte_carlo_scenario(n = 50, sds = c(1, 2, 4, 8), effect_delta = 0)


# ==============================================================================
# PART 2: STATISTICAL POWER SCENARIOS (Delta = 1.0)
# ==============================================================================
cat("\n=======================================================\n")
cat("      PART 2: EVALUATING STATISTICAL POWER             \n")
cat("=======================================================\n")

# Scenario 5: n = 20, Homoscedastic Alternative
power_n20_homo <- run_monte_carlo_scenario(n = 20, sds = c(1, 1, 1, 1), effect_delta = 1.0)

# Scenario 6: n = 20, Heteroscedastic Alternative (Severe)
power_n20_hetero <- run_monte_carlo_scenario(n = 20, sds = c(1, 2, 4, 8), effect_delta = 1.0)

# Scenario 7: n = 50, Homoscedastic Alternative
power_n50_homo <- run_monte_carlo_scenario(n = 50, sds = c(1, 1, 1, 1), effect_delta = 1.0)

# Scenario 8: n = 50, Heteroscedastic Alternative (Severe)
power_n50_hetero <- run_monte_carlo_scenario(n = 50, sds = c(1, 2, 4, 8), effect_delta = 1.0)


# ==============================================================================
# CONSOLIDATED SUMMARY: FULL TABLE 1 REPRODUCTION
# ==============================================================================
cat("\n=================================================================================\n")
cat("                  COMPLETE REPRODUCTION OF TABLE 1 (M = 2000)                    \n")
cat("=================================================================================\n")
cat(sprintf("%-10s | %-12s | %-15s | %-18s | %-15s\n", 
            "Sample (n)", "Variances", "Test", "Type I Error (Null)", "Power (Delta=1)"))
cat("---------------------------------------------------------------------------------\n")

# n = 20, Homo
cat(sprintf("%-10s | %-12s | %-15s | %-18.3f | %-15.3f\n", "n = 20", "1:1:1:1", "Classic ANOVA", res_n20_homo_null[1], power_n20_homo[1]))
cat(sprintf("%-10s | %-12s | %-15s | %-18.3f | %-15.3f\n", "", "", "Welch ANOVA", res_n20_homo_null[2], power_n20_homo[2]))
cat(sprintf("%-10s | %-12s | %-15s | %-18.3f | %-15.3f\n", "", "", "Perm F_ws", res_n20_homo_null[3], power_n20_homo[3]))
cat("---------------------------------------------------------------------------------\n")

# n = 20, Hetero
cat(sprintf("%-10s | %-12s | %-15s | %-18.3f | %-15.3f\n", "n = 20", "1:4:16:64", "Classic ANOVA", res_n20_hetero_null[1], power_n20_hetero[1]))
cat(sprintf("%-10s | %-12s | %-15s | %-18.3f | %-15.3f\n", "", "", "Welch ANOVA", res_n20_hetero_null[2], power_n20_hetero[2]))
cat(sprintf("%-10s | %-12s | %-15s | %-18.3f | %-15.3f\n", "", "", "Perm F_ws", res_n20_hetero_null[3], power_n20_hetero[3]))
cat("---------------------------------------------------------------------------------\n")

# n = 50, Homo
cat(sprintf("%-10s | %-12s | %-15s | %-18.3f | %-15.3f\n", "n = 50", "1:1:1:1", "Classic ANOVA", res_n50_homo_null[1], power_n50_homo[1]))
cat(sprintf("%-10s | %-12s | %-15s | %-18.3f | %-15.3f\n", "", "", "Welch ANOVA", res_n50_homo_null[2], power_n50_homo[2]))
cat(sprintf("%-10s | %-12s | %-15s | %-18.3f | %-15.3f\n", "", "", "Perm F_ws", res_n50_homo_null[3], power_n50_homo[3]))
cat("---------------------------------------------------------------------------------\n")

# n = 50, Hetero
cat(sprintf("%-10s | %-12s | %-15s | %-18.3f | %-15.3f\n", "n = 50", "1:4:16:64", "Classic ANOVA", res_n50_hetero_null[1], power_n50_hetero[1]))
cat(sprintf("%-10s | %-12s | %-15s | %-18.3f | %-15.3f\n", "", "", "Welch ANOVA", res_n50_hetero_null[2], power_n50_hetero[2]))
cat(sprintf("%-10s | %-12s | %-15s | %-18.3f | %-15.3f\n", "", "", "Perm F_ws", res_n50_hetero_null[3], power_n50_hetero[3]))
cat("=================================================================================\n")
# ==============================================================================
# Script: real_data_benchmarks.R
# Purpose: Reproduce Table 2 of the manuscript across the three empirical 
#          benchmarks: InsectSprays (n=12), CO2 (n=21), and iris (n=50).
# ==============================================================================

# Ensure the standalone engine is sourced
source("perm_fws.R")

set.seed(110)
B_PERM <- 9999  # High-precision permutation count for empirical data

cat("======================================================================\n")
cat("      REPRODUCING TABLE 2: EMPIRICAL BENCHMARK COMPARISONS            \n")
cat("======================================================================\n\n")

# ------------------------------------------------------------------------------
# 1. InsectSprays Dataset (k = 6, n_i = 12) - Small-Sample Fragility
# ------------------------------------------------------------------------------
data("InsectSprays")
bart_sprays <- bartlett.test(count ~ spray, data = InsectSprays)
p_classic_sprays <- summary(aov(count ~ spray, data = InsectSprays))[[1]][["Pr(>F)"]][1]
p_welch_sprays <- oneway.test(count ~ spray, data = InsectSprays, var.equal = FALSE)$p.value
res_perm_sprays <- perm_fws_test(InsectSprays$count, InsectSprays$spray, B = B_PERM)

cat("1. InsectSprays (n = 12 per group):\n")
cat("   Bartlett p-value: ", format.pval(bart_sprays$p.value, digits = 4), "\n")
cat("   Classic ANOVA p:  ", format.pval(p_classic_sprays, digits = 4), "\n")
cat("   Welch ANOVA p:    ", format.pval(p_welch_sprays, digits = 4), "\n")
cat("   Perm F_ws T_obs:  ", round(res_perm_sprays$statistic, 4), "\n")
cat("   Perm F_ws p-value:", format.pval(res_perm_sprays$p_value, digits = 4), "\n\n")

# ------------------------------------------------------------------------------
# 2. CO2 Interaction Dataset (k = 4, n_i = 21) - Intermediate Transition
# ------------------------------------------------------------------------------
data("CO2")
CO2$Group <- interaction(CO2$Type, CO2$Treatment)
bart_co2 <- bartlett.test(uptake ~ Group, data = CO2)
p_classic_co2 <- summary(aov(uptake ~ Group, data = CO2))[[1]][["Pr(>F)"]][1]
p_welch_co2 <- oneway.test(uptake ~ Group, data = CO2, var.equal = FALSE)$p.value
res_perm_co2 <- perm_fws_test(CO2$uptake, CO2$Group, B = B_PERM)

cat("2. CO2 Interaction (n = 21 per group):\n")
cat("   Bartlett p-value: ", format.pval(bart_co2$p.value, digits = 4), "\n")
cat("   Classic ANOVA p:  ", format.pval(p_classic_co2, digits = 4), "\n")
cat("   Welch ANOVA p:    ", format.pval(p_welch_co2, digits = 4), "\n")
cat("   Perm F_ws T_obs:  ", round(res_perm_co2$statistic, 4), "\n")
cat("   Perm F_ws p-value:", format.pval(res_perm_co2$p_value, digits = 4), "\n\n")

# ------------------------------------------------------------------------------
# 3. Iris Dataset (k = 3, n_i = 50) - Mature, Fully Calibrated Regime
# ------------------------------------------------------------------------------
data("iris")
bart_iris <- bartlett.test(Sepal.Length ~ Species, data = iris)
p_classic_iris <- summary(aov(Sepal.Length ~ Species, data = iris))[[1]][["Pr(>F)"]][1]
p_welch_iris <- oneway.test(Sepal.Length ~ Species, data = iris, var.equal = FALSE)$p.value
res_perm_iris <- perm_fws_test(iris$Sepal.Length, iris$Species, B = B_PERM)

cat("3. Iris Sepal Length (n = 50 per group):\n")
cat("   Bartlett p-value: ", format.pval(bart_iris$p.value, digits = 4), "\n")
cat("   Classic ANOVA p:  ", format.pval(p_classic_iris, digits = 4), "\n")
cat("   Welch ANOVA p:    ", format.pval(p_welch_iris, digits = 4), "\n")
cat("   Perm F_ws T_obs:  ", round(res_perm_iris$statistic, 4), "\n")
cat("   Perm F_ws p-value:", format.pval(res_perm_iris$p_value, digits = 4), "\n")
cat("======================================================================\n")
