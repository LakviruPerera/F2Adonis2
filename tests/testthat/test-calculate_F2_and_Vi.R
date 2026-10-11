test_that("calculate_F2_and_Vi matches the known Norway reference value", {
  data("norway_pa", package = "F2Adonis2", envir = environment())
  data("norway_metadata", package = "F2Adonis2", envir = environment())

  grp <- factor(norway_metadata$Area)
  dmat_sq <- as.matrix(vegan::vegdist(norway_pa, method = "jaccard"))^2
  n <- nrow(dmat_sq)
  G <- F2Adonis2:::compute_G(dmat_sq)
  H <- F2Adonis2:::get_H(model.matrix(~ grp))
  SS_A <- sum(diag(H %*% G))

  res <- F2Adonis2:::calculate_F2_and_Vi(dmat_sq, grp, SS_A, n)

  expect_equal(res$F2, 13.509911, tolerance = 1e-5)
})
