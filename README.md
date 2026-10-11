# F2Adonis2

F2Adonis2 is a tool for multivariate hypothesis testing, developed to
support clinical metabolomics and high-dimensional omics research. It
addresses the multivariate Behrens-Fisher Problem in dissimilarity-based
analyses (such as PERMANOVA), providing a statistically sound approach
for datasets where group dispersions are heterogeneous. It is only
available on GitHub.

## Motivation

In high-dimensional ecological and metabolomics studies ($P \gg N$),
standard PERMANOVA ($F_1$) methods often yield inflated Type I error
rates when group dispersions are unequal (heteroskedasticity) and group
sizes are unbalanced.

This package implements the methodological solutions proposed by
**Anderson et al. (2017)**. Rather than pooling variance across groups
(which leads to bias), the $F_2$ statistic uses a weighted combination
of group-specific dispersion estimates to maintain level accuracy under
heteroskedastic conditions.

## Installation

You must install the `remotes` package first.

``` r
install.packages("remotes")
```

Then, use `install_github` to install this package.

``` r
remotes::install_github("LakviruPerera/F2Adonis2")
```

## Quick Start

`F2Adonis2` acts as a drop-in replacement for `vegan::adonis2`. It
follows the standard R formula interface and accepts either a data
matrix (distances are computed with `vegdist`) or a precomputed
dissimilarity (`dist`) object on the left-hand side.

``` r
library(F2Adonis2)

# Example data included with the package
data("pkd_metadata", package = "F2Adonis2")
data("pkd_metabolites", package = "F2Adonis2")

pkd_meta <- data.frame(Group = factor(pkd_metadata$Group))

# Run the F2 test (use more permutations, e.g. 999, for final results)
fit_f2 <- F2_adonis2(
  as.matrix(pkd_metabolites) ~ Group,
  data = pkd_meta,
  method = "euclidean",
  permutations = 99
)

# View results
print(fit_f2)
```

## Support and contributing

Questions and bug reports are welcome on the
[issue tracker](https://github.com/LakviruPerera/F2Adonis2/issues).
See [CONTRIBUTING.md](CONTRIBUTING.md) for how to report a bug or propose a change.

## References

Anderson, M.J., Walsh, D.C.I., Robert Clarke, K., Gorley, R.N. and
Guerra-Castro, E. (2017), Some solutions to the multivariate
Behrens–Fisher problem for dissimilarity-based analyses. Aust. N. Z. J.
Stat., 59: 57-79. <https://doi.org/10.1111/anzs.12176>

Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin P, O’Hara
R, Solymos P, Stevens M, Szoecs E, Wagner H, Barbour M, Bedward M,
Bolker B, Borcard D, Borman T, Carvalho G, Chirico M, De Caceres M,
Durand S, Evangelista H, FitzJohn R, Friendly M, Furneaux B, Hannigan G,
Hill M, Lahti L, Martino C, McGlinn D, Ouellette M, Ribeiro Cunha E,
Smith T, Stier A, Ter Braak C, Weedon J (2026). vegan: Community Ecology
Package. R package version 2.8-0, <https://vegandevs.github.io/vegan/>.
