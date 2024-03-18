# Hardy-Weinberg Equilibrium Package

This R package provides tools for conducting exact tests under the Hardy-Weinberg Equilibrium (HWE) with adjustments for inbreeding, represented by the parameter theta. This is based on the work of [Wigginton](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1199378/) and [Weir](https://academic.oup.com/genetics/article/180/3/1609/6063905). 

## Installation
You can install the latest version of the package directly from GitHub:
```
# install.packages("devtools") # Uncomment if devtools is not installed
devtools::install_github("frhl/hardyr")
```

## Details
The `hardyr` package allows for precise hypothesis testing concerning Hardy-Weinberg equilibrium (HWE), even with large sample sizes. Below are examples demonstrating how to perform exact tests for P-values and conduct power calculations across a range of population structures and minor allele frequencies (MAFs).

## License
This package is licensed under the MIT License - see the LICENSE file for details.


