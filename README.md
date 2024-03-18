# Hardy-Weinberg Equilibrium Package

This R package provides tools for conducting exact tests under the Hardy-Weinberg Equilibrium (HWE) with adjustments for inbreeding, represented by the parameter theta. This is based on the work of [Wigginton](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1199378/) and [Weir](https://academic.oup.com/genetics/article/180/3/1609/6063905). 

## Installation
You can install the latest version of the package directly from GitHub:
```
# install.packages("devtools") # Uncomment if devtools is not installed
devtools::install_github("frhl/hardyr")
```

## Usage
The `hardyr` package allows for precise hypothesis testing concerning Hardy-Weinberg equilibrium (HWE), even with large sample sizes. Below are examples demonstrating how to perform exact tests for P-values and conduct power calculations across a range of population structures and minor allele frequencies (MAFs).

### Performing Exact HWE Testing
You can perform an exact test for HWE using the `hwe_exact_test` function. This is particularly useful for analyzing large datasets. In the example below, we specify a population size (N), the count of the minor allele (K), and the count of homozygotes (M). The theta parameter adjusts the test for inbreeding effects, with the option to compute it dynamically using `calc_theta_from_f(pA, f)`.
```
library(hardyr)

# Example data
N <- 100000
K <- 1000  # Minor allele count
M <- 24  # Homozygotes

# Perform the HWE exact test
hwe_exact_test(N, K, nM, theta = 4, alternative = "less")
```
To account for inbreeding when testing for HWE, you can dynamically calculate the theta parameter as shown:
```
theta <- calc_theta_from_f(pA, f)  # pA is the allele frequency, and f is the inbreeding coefficient
```
### Power Calculation Across MAF Spectrum
For understanding the statistical power of HWE tests across different population sizes and MAFs, we can simulate scenarios using a combination of hardyr and ggplot2 for visualization. The following script demonstrates power calculation for multiple populations with varying F across a range of minor allele frequencies:
```
library(hardyr)
library(ggplot2)

# Define sequences for population sizes and inbreeding coefficients
N_seq <- c(10000, 50000)
f_seq <- c(0, 0.02, 0.05)

# Perform power calculations
out <- do.call(rbind, lapply(f_seq, function(f) {
  do.call(rbind, lapply(N_seq, function(N) {
    nA_seq <- seq(0, (N*2)*0.1, by=250)
    power_vec <- vector()
    do.call(rbind, lapply(nA_seq, function(nA) {
      pA <- nA / (N * 2)
      theta <- calc_theta_from_f(pA, f)
      power <- hwe_exact_power(N, K=nA, theta=theta, alternative="less", sig.level = 0.05)
      power_vec <- c(power_vec, power)
      data.table(N, pA, nA, f, power)
    }))
  }))
}))

# Prepare data for plotting
out$N_label <- factor(sprintf("%sK", out$N / 1000))
out$F_label <- factor(sprintf("F=%.2f", out$f))

# Plot the results
ggplot(out, aes(x=pA, y=power)) +
  geom_point() +
  geom_line() +
  facet_grid(F_label ~ N_label, labeller = label_both) +
  theme_minimal() +
  labs(x = "Minor Allele Frequency (pA)", y = "Power", title = "Power Calculation Across MAF Spectrum")
````
This script first calculates the power for detecting deviations from HWE across different population structures and sizes and frequencies of the minor allele. It then plots these calculations, helping to visualize the power of the HWE test under various conditions.

<img src="img/sim_power01.png" width=50% height=50%>

## License
This package is licensed under the MIT License - see the LICENSE file for details.


