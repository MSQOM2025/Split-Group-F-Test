# ==============================================================================
#
#    R Script for Power Surface Simulation of the Perm F_ws Test
#  
#
# ==============================================================================

# 1. SETUP
# --- Check and load required packages ---
required_packages <- c("tidyverse", "future", "furrr", "viridis")
not_installed <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if(length(not_installed) > 0) { install.packages(not_installed) }

library(tidyverse)
library(future)
library(furrr)
library(viridis) # For the nice color scale in the plot
select <- dplyr::select
set.seed(123)
# --- Source the standalone function file ---
# This loads the correct, final `perm_fws_test` function.
if (!file.exists("perm_fws_test.R")) {
  stop("Critical file missing: 'perm_fws_test.R'. Please place it in the same directory.")
}
source("perm_fws_test.R")


# 2. CONTROL PANEL
SIM_REPLICATIONS <- 500  # Number of Monte Carlo runs for each grid point
PERMS_FOR_SURFACE <- 199 # Use fewer permutations for a faster surface plot


# 3. DATA GENERATION FUNCTION
generate_data <- function(n, effect_size) {
    k <- 4
    sds <- c(1, 2, 4, 8) # Heteroscedastic case
    means <- c(0, 0, 0.5, 1.0) * effect_size 
    
    observations <- vector("list", k)
    for (i in 1:k) { observations[[i]] <- rnorm(n, means[i], sds[i]) }
    tibble(observation = unlist(observations), group = factor(rep(1:k, each = n)))
}


# 4. SIMULATION EXECUTION
# --- Define the simulation grid ---
simulation_grid <- expand_grid(
  n = c(20, 30, 40, 50, 60, 80, 100),
  effect_size = seq(0, 1.5, by = 0.25)
)

# --- Set up parallel plan ---
num_cores <- availableCores() - 1
plan(multisession, workers = num_cores)
cat(paste("Parallel processing enabled on", num_cores, "cores.\n"))

# --- Run the simulation in parallel ---
cat("Starting power surface simulation...\n")
start_time <- Sys.time()

surface_results <- future_map_dfr(
  .x = 1:nrow(simulation_grid),
  ~ {
    # This block of code is sent to a parallel worker.
    # It receives the index `.x`.
    params <- as.list(simulation_grid[.x, ])
    
    p_values <- replicate(SIM_REPLICATIONS, {
      sim_data <- generate_data(n = params$n, effect_size = params$effect_size)
      
      # Use the robust, sourced function and extract p-value
      test_result <- perm_fws_test(sim_data, perms = PERMS_FOR_SURFACE) 
      test_result$p.value
    })
    
    # Return a single-row tibble with the results for this grid point.
    tibble(
      n = params$n,
      effect_size = params$effect_size,
      power = mean(p_values < 0.05, na.rm = TRUE)
    )
  },
  .progress = TRUE,
  # --- THE DEFINITIVE FIX ---
  .options = furrr_options(
    seed = TRUE, 
    # Global variables/functions the workers need to know about
    globals = c("SIM_REPLICATIONS", "PERMS_FOR_SURFACE", 
                "simulation_grid", "generate_data", 
                "perm_fws_test", "calculate_fws_unb_statistic"),
    # Packages the workers need to load in their sessions
    packages = c("tidyverse")
  )
)

end_time <- Sys.time()
cat("Power surface simulation finished. Total time:", format(end_time - start_time), "\n")


# 5. SAVE AND PLOT RESULTS
write_csv(surface_results, "power_surface_data.csv")
cat("Power surface data saved to 'power_surface_data.csv'\n")

# --- Create the Power Surface Plot (Heatmap) ---
power_surface_plot <- ggplot(surface_results, aes(x = n, y = effect_size, fill = power)) +
  geom_tile(color = "white") +
  # Use the viridis color scale for a professional look
  scale_fill_viridis_c(name = "Power", labels = scales::percent, limits=c(0,1)) +
  geom_text(aes(label = sprintf("%.2f", power)), color="white", size=3.5) +
  labs(
    title = "Figure 5: Power Surface of the Perm F-ws Test",
    subtitle = "Performance across sample size and effect size under heteroscedasticity",
    x = "Sample Size per Group (n)",
    y = "Effect Size"
  ) +
  # Ensure all axis ticks are shown
  scale_x_continuous(breaks = unique(surface_results$n)) +
  scale_y_continuous(breaks = unique(surface_results$effect_size)) +
  theme_minimal(base_size = 14)

print(power_surface_plot)
ggsave("plot_power_surface.png", plot = power_surface_plot, width = 10, height = 8, dpi = 300)

cat("\nPower surface plot saved to 'plot_power_surface.png'\n")
