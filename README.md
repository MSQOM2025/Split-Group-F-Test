# The Price of Purity: On the Sample-Size Dependency of a Split-Group Permutation Test for the Behrens-Fisher Problem

[![R](https://img.shields.io/badge/Language-R_%3E%3D_4.0.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository provides the official, fully reproducible **Base R** implementation and simulation pipeline for the manuscript:

> **The Price of Purity: On the Sample-Size Dependency of a Split-Group Permutation Test for the Behrens-Fisher Problem**  
> *Mohsen Salehi (Department of Statistics, University of Qom)*

---

## Overview

Comparing group means across $k \ge 2$ populations under variance heterogeneity (the Behrens-Fisher problem) remains a foundational challenge in statistical inference. While Welch's ANOVA provides an established parametric approximation based on Satterthwaite degrees of freedom, non-parametric permutation tests traditionally face exchangeability obstacles under heteroscedasticity.

This repository implements the **Split-Group Permutation Test ($\text{Perm } F_{ws}$)**, a resampling method that uses intra-group data-splitting to construct an asymptotically studentized dispersion statistic:

$$T = \sum_{i=1}^k w_i (\bar{X}_i^* - \bar{X}_w^*)^2$$

### Key Methodological Findings
1. **The Price of Purity (Small-Sample Behavior):** In small-sample regimes ($n = 20$), halving the group sample sizes leaves subgroup variances with limited degrees of freedom ($\nu = 9$). The inverse-gamma distribution of the precision weights exhibits heavy polynomial tails, resulting in a mild Type~I error elevation ($6.4\%$ under a severe $64:1$ variance ratio).
2. **Asymptotic Stabilization ($n \ge 40$ to $50$):** As sample sizes increase, precision weights contract around their population parameters. At $n = 50$, the empirical Type~I error settles at an optimal **$5.8\%$**, achieving nominal calibration while delivering empirical power (**$13.6\%$**) that matches and slightly exceeds Welch's ANOVA (**$13.4\%$**).
3. **Pure Base R Performance:** The core engine is implemented without external dependencies, evaluating $9{,}999$ permutations on typical experimental samples in approximately **$0.05$ seconds**.

---

## Repository Contents

```text
├── perm_fws.R               # Standalone, dependency-free Base R test implementation
├── simulation_table1.R      # Monte Carlo pipeline reproducing Table 1 (M = 2000 replications)
├── real_data_benchmarks.R   # Complete reproduction of Table 2 (InsectSprays, CO2, iris)
└── README.md                # Project documentation and quick-start guide
