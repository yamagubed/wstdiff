# wstdiff: Welch-Satterthwaite Approximation for t-Distribution Differences

## Overview

The `wstdiff` package implements the Welch-Satterthwaite approximation for 
differences of non-standardized t-distributed random variables in both 
univariate and multivariate settings.

## Installation

```r
# Install from GitHub (once available)
# devtools::install_github("yourusername/wstdiff")

# Or install locally
devtools::install_local("path/to/wstdiff")
```

## Usage

### Univariate Case

```r
library(wstdiff)

# Basic example
result <- ws_tdiff_univariate(
  mu1 = 0, sigma1 = 1, nu1 = 10,
  mu2 = 0, sigma2 = 1.5, nu2 = 15
)
print(result)

# Distribution functions
dtdiff(0, result)           # Density
ptdiff(0, result)           # CDF
qtdiff(c(0.025, 0.975), result)  # Quantiles
samples <- rtdiff(1000, result)  # Random generation
```

### Multivariate Case (Independent Components)

```r
result <- ws_tdiff_multivariate_independent(
  mu1 = c(0, 1), 
  sigma1 = c(1, 1.5), 
  nu1 = c(10, 12),
  mu2 = c(0, 0), 
  sigma2 = c(1.2, 1), 
  nu2 = c(15, 20)
)
```

### Multivariate Case (General Covariance)

```r
Sigma1 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
Sigma2 <- matrix(c(1.5, 0.5, 0.5, 1.2), 2, 2)

result <- ws_tdiff_multivariate_general(
  mu1 = c(0, 1), 
  Sigma1 = Sigma1, 
  nu1 = 10,
  mu2 = c(0, 0), 
  Sigma2 = Sigma2, 
  nu2 = 15
)
```

## Reference

Yamaguchi, Y., Homma, G., Maruo, K., & Takeda, K. 
Welch-Satterthwaite Approximation for Difference of Non-Standardized 
t-Distributed Variables. (unpublished).

## License

MIT License
