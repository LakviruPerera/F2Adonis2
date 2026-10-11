test_that("singleton groups stop with a clear error", {
  set.seed(1)
  Y <- matrix(rnorm(7 * 5), nrow = 7)
  grp <- factor(c("A", "A", "A", "B", "B", "B", "C"))
  expect_error(
    F2_adonis2(Y ~ grp, data = data.frame(grp = grp),
               permutations = 19, plot.dist = FALSE),
    "at least 2 observations"
  )
})

test_that("a non-factor predictor is rejected", {
  set.seed(1)
  Y <- matrix(rnorm(12 * 4), nrow = 12)
  grp <- rep(c("A", "B", "C"), each = 4)
  expect_error(
    F2_adonis2(Y ~ grp, data = data.frame(grp = grp),
               permutations = 19, plot.dist = FALSE),
    "must be a factor"
  )
})

test_that("V_k matches a brute-force loop calculation", {
  set.seed(2)
  Y <- matrix(rnorm(8 * 3), nrow = 8)
  g <- factor(c(rep("A", 3), rep("B", 5)))
  dsq <- as.matrix(dist(Y))^2

  brute <- sapply(levels(g), function(l) {
    idx <- which(g == l)
    s <- 0
    for (i in seq_along(idx)) for (j in seq_along(idx)) {
      if (i < j) s <- s + dsq[idx[i], idx[j]]
    }
    s / (length(idx) * (length(idx) - 1))
  })

  res <- calculate_F2_and_Vi(dsq, g, SS.A = 1, n = nrow(dsq))
  expect_equal(unname(res$Vi), unname(brute))
})

test_that("F2 reproduces the known value on the Norway data", {
  data("norway_pa", package = "F2Adonis2", envir = environment())
  data("norway_metadata", package = "F2Adonis2", envir = environment())
  fit <- F2_adonis2(norway_pa ~ Area, data = norway_metadata,
                    method = "jaccard", permutations = 99,
                    plot.dist = FALSE)
  expect_equal(fit["Model", "F2"], 13.5099, tolerance = 1e-4)
})

test_that("a dist object gives the same F2 as the raw data", {
  set.seed(3)
  Y <- matrix(rnorm(12 * 4), nrow = 12)
  d <- data.frame(grp = factor(rep(c("A", "B", "C"), each = 4)))
  fit_raw  <- F2_adonis2(Y ~ grp, data = d, method = "euclidean",
                         permutations = 19, plot.dist = FALSE)
  fit_dist <- F2_adonis2(dist(Y) ~ grp, data = d,
                         permutations = 19, plot.dist = FALSE)
  expect_equal(fit_raw["Model", "F2"], fit_dist["Model", "F2"])
})

test_that("bootstrap and bias adjustment add p-value columns", {
  set.seed(4)
  Y <- matrix(rnorm(12 * 6), nrow = 12)
  d <- data.frame(grp = factor(rep(c("A", "B", "C"), each = 4)))
  fit <- F2_adonis2(Y ~ grp, data = d, method = "euclidean",
                    permutations = 19, bootstrap = TRUE,
                    bias.adjust = TRUE, plot.dist = FALSE)
  expect_true("Pr(>F2.boot)" %in% names(fit))
  expect_true("Pr(>F2.ba)" %in% names(fit))
  expect_true(fit["Model", "Pr(>F2.boot)"] > 0 &&
                fit["Model", "Pr(>F2.boot)"] <= 1)
})

test_that("the same seed gives the same p-value", {
  Y <- matrix(rnorm(12 * 4), nrow = 12)
  d <- data.frame(grp = factor(rep(c("A", "B", "C"), each = 4)))
  set.seed(10)
  p1 <- F2_adonis2(Y ~ grp, data = d, permutations = 49,
                   plot.dist = FALSE)["Model", "Pr(>F2)"]
  set.seed(10)
  p2 <- F2_adonis2(Y ~ grp, data = d, permutations = 49,
                   plot.dist = FALSE)["Model", "Pr(>F2)"]
  expect_equal(p1, p2)
})
