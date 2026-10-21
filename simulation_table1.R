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