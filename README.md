# Gene-level Hardy-Weinberg Equilibrium Packag

This R package provides tools for conducting exact tests under the Hardy-Weinberg Equilibrium (HWE) with adjustments for inbreeding, represented by the parameter theta. The implementation is the work of [Wigginton](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1199378/) and [Weir](https://academic.oup.com/genetics/article/180/3/1609/6063905). 

## Features
* Exact Testing for HWE: Perform precise statistical tests to assess if a population is in Hardy-Weinberg equilibrium.
* Inbreeding Coefficient Adjustment: Incorporate the inbreeding coefficient through the theta parameter to adjust the equilibrium expectations.
* Logarithmic Calculations: Utilize logarithmic transformations for numerical stability and handling large datasets.

## Installation
You can install the latest version of the package directly from GitHub:
```
# install.packages("devtools") # Uncomment if devtools is not installed
devtools::install_github("frhl/hardyr")
```

## License
This package is licensed under the MIT License - see the LICENSE file for details.


