# --- Addendum Script to Visualize Weight Instability ---
library(tidyverse)

# A simplified function to get a vector of weights
get_weights <- function(data) {
  subgroup_stats <- data %>%
    group_by(group) %>%
    mutate(subgroup = sample(rep(1:2, length.out = n()))) %>%
    group_by(group, subgroup) %>%
    summarise(n = n(), var = var(observation), .groups = 'drop')
  
  wide_stats <- subgroup_stats %>%
    pivot_wider(id_cols = group, names_from = subgroup, values_from = c(n, var))
  
  if (any(is.na(wide_stats))) return(NA_real_)
  
  # Return the standard inverse variance weights for simplicity and demonstration
  w_i <- (1 / wide_stats$var_1) + (1 / wide_stats$var_2)
  return(w_i)
}

# --- Simulate weights for n=20 and n=50 ---
set.seed(123)
n20_weights <- replicate(500, {
  sim_data <- tibble(
      observation = c(rnorm(20, 0, 1), rnorm(20, 0, 2), rnorm(20, 0, 4), rnorm(20, 0, 8)),
      group = factor(rep(1:4, each = 20))
  )
  get_weights(sim_data)
})

n50_weights <- replicate(500, {
  sim_data <- tibble(
      observation = c(rnorm(50, 0, 1), rnorm(50, 0, 2), rnorm(50, 0, 4), rnorm(50, 0, 8)),
      group = factor(rep(1:4, each = 50))
  )
  get_weights(sim_data)
})

# --- Prepare data for plotting ---
plot_data <- bind_rows(
  tibble(weight = as.numeric(n20_weights), n = "n=20"),
  tibble(weight = as.numeric(n50_weights), n = "n=50")
) %>% filter(!is.na(weight))

# --- CORRECTED CODE to Generate Weight Instability Plot ---
# (with 'linewidth' aesthetic to remove the warning)

weight_plot <- ggplot(plot_data, aes(x = weight, fill = n)) +
  geom_density(alpha = 0.6) +
  scale_x_log10(labels = scales::label_number_auto()) +
  labs(
    title = "Figure 4: Distribution of Simulated Inverse-Variance Weights",
    subtitle = "Instability of weights for small sample sizes",
    x = "Weight Value (w_i, log scale)",
    y = "Density",
    fill = "Sample Size"
  ) +
  theme_bw(base_size = 14) +
  # --- THE FIX IS HERE: `size` changed to `linewidth` ---
  geom_vline(xintercept = median(plot_data$weight[plot_data$n=="n=50"]), 
             linetype="dashed", color="blue", linewidth=1) +
  geom_vline(xintercept = median(plot_data$weight[plot_data$n=="n=20"]), 
             linetype="dashed", color="red", linewidth=1)
  
print(weight_plot)
ggsave("plot_weight_instability.png", plot = weight_plot)