# F2Adonis2 (development version)

* Added `NEWS.md`, `CONTRIBUTING.md`, and a support section in the README.
* Added a GitHub Actions workflow that runs `R CMD check` on Linux, macOS, and Windows.

# F2Adonis2 0.1.0

* First tagged release.
* `F2_adonis2()` runs PERMANOVA with the F2 statistic (Anderson et al. 2017), for comparing groups when dispersions differ and group sizes are unequal.
* p-values by permutation, separate-sample bootstrap, or bias-adjusted bootstrap.
* Plot of the permutation null distribution.
* Includes the Norway benthic macrofauna and PKD urinary metabolomics datasets, and two vignettes.
* The dispersion calculation is much faster than in earlier development versions, with matching results.
* Groups with fewer than two observations now give a clear error.
