# Code for A Split-Group Permutation F-Test for Robust ANOVA under Heteroscedasticity

This repository contains the R scripts necessary to reproduce all simulations, analyses, and figures presented in the manuscript, A Split-Group Permutation F-Test for Robust ANOVA under Heteroscedasticity by [Your Name], (Year).

## Description

This research introduces and evaluates the Permutation Split-Group F-Test (`Perm F_ws`), a novel non-parametric test for comparing k=2 group means under variance heterogeneity (the Behrens-Fisher problem). The core finding is that the test's validity is dependent on sample size; it is anti-conservative in small samples but becomes a valid, powerful, and theoretically exact alternative to Welch's ANOVA in moderate-to-large samples.

This repository provides the code to replicate the simulation studies and the real-data applications that support these conclusions.

## Repository Contents

This repository contains three primary R scripts

1.  `simulation_main.R` This is the main simulation script. It performs the large-scale Monte Carlo study that compares the Type I error rate and statistical power of the `Perm F_ws` test against the classic ANOVA F-test and Welch's F-test across various scenarios. It also performs the deep dive analysis to generate the data for the p-value accuracy (Q-Q) plots. Warning This is a computationally intensive script that may take a long time to run.

2.  `simulation_addendum_weights.R` This is a smaller, targeted simulation script that generates the data for the Weight Instability Plot (Figure 4 in the manuscript). It visualizes the distribution of the inverse-variance weights for small (`n=20`) vs. moderate (`n=50`) sample sizes to provide evidence for the source of the small-sample instability.

3.  `perm_fws_test.R` This is a clean, well-commented, standalone implementation of our proposed statistical test. It contains the core `perm_fws_test()` function and a working example. This script is provided to allow other researchers to easily apply our method to their own data.

## System Requirements

The code was developed and tested using R version 4.4.x. The following R packages are required and will be automatically checked for (and prompted for installation if missing) by the main simulation script

   `tidyverse` For data manipulation (`dplyr`) and plotting (`ggplot2`).
   `future` For setting up the parallel processing backend.
   `furrr` For executing the parallel simulation loops.

## How to Run

1.  Main Simulation
       Open the `simulation_main.R` script.
       At the top of the script, in the SIMULATION CONTROL PANEL, you can set `TEST_MODE - TRUE` for a very fast (~5-10 minute) test run or `TEST_MODE - FALSE` for the full, publication-quality run (which may take over 24 hours).
       The script will automatically configure parallel processing to use a safe number of CPU cores.
       Upon completion, it will save the full numerical results to `simulation_results_fws_final.csv` and generate all manuscript figures (`plot_type1_error.png`, etc.) in the working directory.

2.  Standalone Function
       To use the test on your own data, you can simply `source(perm_fws_test.R)` and then call the `perm_fws_test(your_data_frame)` function as shown in the example at the bottom of the script file.

## Citation

If you use the code or methods from this repository in your own work, please cite the original manuscript

[Full citation to be added here upon publication]
