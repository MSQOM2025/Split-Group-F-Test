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