#' Permutational Multivariate Analysis of Variance Using the F2 Statistic
#'
#' Analysis of variance using distance matrices for partitioning variation
#' among sources. Unlike standard PERMANOVA, this function uses the F2 statistic
#' (Anderson 2017) which accounts for heterogeneous dispersions (Behrens-Fisher problem).
#'
#' @param formula Model formula. The LHS must be a community data matrix or a
#'     dissimilarity matrix. The RHS must define a single factor predictor.
#' @param data A data frame containing the independent variables.
#' @param permutations Number of permutations for the test, or a permutation matrix.
#' @param method Pairwise distance method used in \code{vegdist}. Common options include
#'     "jaccard" (for binary data), "bray" (for abundance data), or "euclidean".
#' @param sqrt.dist Logical; take square root of dissimilarities to "euclidify" them.
#' @param add Add a constant to dissimilarities to avoid negative eigenvalues.
#'     Options are "lingoes" or "cailliez".
#' @param by Validation restricted to \code{NULL}. F2 currently only supports global tests.
#' @param parallel Number of parallel processes to use.
#' @param na.action Handling of missing values on the RHS.
#' @param strata Groups within which to constrain permutations.
#' @param bootstrap Logical; perform separate-sample residual bootstrapping in
#'     dissimilarity space (Anderson 2017).
#' @param bias.adjust Logical; apply empirical bias correction to bootstrap variances.
#' @param plot.dist Logical; if TRUE (default), plot the permutation distribution of the F2 statistic.
#' @param ... Other arguments passed to \code{vegdist}.
#'
#' @details
#' \code{F2_adonis2} implements the multivariate F2 statistic and the
#' matrix-based bootstrap procedure derived in Anderson et al. (2017).
#'
#' @importFrom vegan vegdist
#' @importFrom stats formula terms model.matrix as.dist na.fail
#'
#' @references
#' Anderson, M. J., Walsh, D. C., Clarke, K. R., Gorley, R. N., & Guerra-Castro, E. (2017).
#' Some solutions to the multivariate Behrens–Fisher problem for dissimilarity-based
#' analyses. \emph{Australian & New Zealand Journal of Statistics}.
#'
#' @author Lakviru Perera (based on logic by Marti J. Anderson)
#'
#' @export
F2_adonis2 <-
  function(formula, data, permutations = 999, method = "euclidean",
           sqrt.dist = FALSE, add = FALSE, by = NULL,
           parallel = getOption("mc.cores"), na.action = na.fail,
           strata = NULL, bootstrap = FALSE, bias.adjust = FALSE,
           plot.dist = TRUE, ...) # Default changed to TRUE
  {
    if (missing(data)) data <- parent.frame()
    else data <- eval(match.call()$data, parent.frame(), enclos = environment(formula))

    if (!is.null(by)) {
      stop("\nError: F2 PERMANOVA currently only supports by = NULL.", call. = FALSE)
    }

    formula_obj <- formula(terms(formula, data = data))
    term_labels <- attr(terms(formula_obj), "term.labels")

    if (length(term_labels) > 1) {
      stop("F2_adonis2 supports only ONE predictor variable.", call. = FALSE)
    }

    predictor_val <- tryCatch({
      eval(parse(text = term_labels[1]), envir = data)
    }, error = function(e) data[[term_labels[1]]])

    if (!is.factor(predictor_val)) {
      stop("The predictor variable must be a factor.", call. = FALSE)
    }

    lhs <- eval(formula_obj[[2]], envir = parent.frame(), enclos = environment(formula_obj))
    original_data <- if (!inherits(lhs, "dist") && !isSymmetric(unname(as.matrix(lhs)))) as.matrix(lhs) else NULL
    original_method <- method

    if (!inherits(lhs, "dist")) lhs <- vegdist(as.matrix(lhs), method = method, ...)

    d <- vegan:::ordiParseFormula(formula = formula_obj, data = data, na.action = na.action, subset = NULL, X = lhs)
    lhs_val <- d$X

    if (sqrt.dist) lhs_val <- sqrt(lhs_val)
    if (is.logical(add) && add) add <- "lingoes"
    if (is.character(add)) {
      add <- match.arg(add, c("lingoes", "cailliez"))
      if (add == "lingoes") {
        ac <- vegan:::addLingoes(as.matrix(lhs_val))
        lhs_val <- sqrt(lhs_val^2 + 2 * ac)
      } else if (add == "cailliez") {
        ac <- vegan:::addCailliez(as.matrix(lhs_val))
        lhs_val <- lhs_val + ac
      }
    }

    sol <- F2_adonis0(lhs_val, d$Y, d$Z, original_data, original_method)
    sol$formula <- match.call()
    sol$terms <- d$terms
    sol$terminfo <- vegan:::ordiTerminfo(d, data)

    perm <- vegan:::getPermuteMatrix(permutations, NROW(lhs_val), strata = strata)
    out <- anova(sol, permutations = perm, by = by, parallel = parallel, bootstrap = bootstrap, bias.adjust = bias.adjust)

    att <- attributes(out)
    att$heading[1] <- "Permutation test for F2-based PERMANOVA"
    attributes(out) <- att

    # Trigger plot automatically by default
    if (plot.dist) plot(out)

    return(out)
  }

F2_adonis0 <- function(lhs, X, Z, original_data = NULL, original_method = "euclidean") {
  G <- vegan:::initDBRDA(lhs)
  sol <- list(Ybar = G, tot.chi = sum(diag(G)), method = "F2_adonis",
              dmat_squared = as.matrix(lhs)^2, original_data = original_data,
              original_method = original_method, X_original = X)

  if (!is.null(Z) && ncol(Z) > 0) {
    Z  <- scale(Z, scale = FALSE)
    QZ <- qr(Z)
    G  <- qr.resid(QZ, G); G <- qr.resid(QZ, t(G))
    sol$pCCA <- list(rank = QZ$rank, tot.chi = sum(diag(G)), QR = QZ)
    X  <- cbind(Z, X)
  }

  X  <- scale(X, scale = FALSE); qrhs <- qr(X)
  Gfit <- qr.fitted(qrhs, G); Gfit <- qr.fitted(qrhs, t(Gfit))
  Gres <- qr.resid(qrhs, G); Gres <- qr.resid(qrhs, t(Gres))

  sol$CCA <- if(qrhs$rank > 0) list(rank = qrhs$rank, tot.chi = sum(diag(Gfit))) else NULL
  sol$CA  <- list(rank = nrow(lhs) - qrhs$rank - 1, tot.chi = sum(diag(Gres)))
  class(sol) <- c("F2_adonis2", "dbrda", "rda", "cca")
  return(sol)
}

#' @export
anova.F2_adonis2 <- function(object, permutations, bootstrap = FALSE, bias.adjust = FALSE, ...) {
  X_orig <- object$X_original
  groups <- if (ncol(X_orig) == 1) as.factor(X_orig[, 1]) else interaction(as.data.frame(X_orig), drop = TRUE)
  G_obs <- object$Ybar; dmat_sq <- object$dmat_squared; n <- nrow(dmat_sq)

  H_obs <- get_H(model.matrix(~ groups)); SS_A_obs <- sum(diag(H_obs %*% G_obs))
  res_obs <- calculate_F2_and_Vi(dmat_sq, groups, SS_A_obs, n)
  F2_obs <- res_obs$F2; Vi_obs <- res_obs$Vi

  n_perms <- nrow(permutations); F2_p <- numeric(n_perms)
  for (i in 1:n_perms) {
    grp_p <- groups[permutations[i, ]]; H_p <- get_H(model.matrix(~ grp_p))
    F2_p[i] <- calculate_F2_and_Vi(dmat_sq, grp_p, sum(diag(H_p %*% G_obs)), n)$F2
  }
  p_perm <- (sum(F2_p >= F2_obs, na.rm = TRUE) + 1) / (n_perms + 1)

  p_boot <- NA; p_ba <- NA
  if (bootstrap) {
    IH <- diag(n) - H_obs
    R_mat <- IH %*% G_obs %*% IH
    grp_levels <- levels(groups); n_per_group <- table(groups)
    H_res_i_list <- lapply(grp_levels, function(g) IH %*% diag(as.numeric(groups == g)) %*% IH)

    F2_b <- numeric(n_perms); Vi_bm <- matrix(NA, n_perms, length(grp_levels))

    for (i in 1:n_perms) {
      beta <- integer(n)
      for (g in grp_levels) { idx_g <- which(groups == g); beta[idx_g] <- sample(idx_g, replace = TRUE) }
      R_boot <- R_mat[beta, beta]; ss_ab  <- sum(diag(H_obs %*% R_boot))
      Vi_boot <- sapply(1:length(grp_levels), function(k) sum(diag(H_res_i_list[[k]] %*% R_boot)) / (n_per_group[k] - 1))
      den_b <- sum((1 - as.numeric(n_per_group) / n) * Vi_boot)
      F2_b[i] <- ss_ab / den_b; Vi_bm[i,] <- Vi_boot
    }
    p_boot <- (sum(F2_b >= F2_obs, na.rm = TRUE) + 1) / (n_perms + 1)

    if (bias.adjust) {
      bias <- colMeans(Vi_bm, na.rm = TRUE) - Vi_obs
      f2_ba <- sapply(1:n_perms, function(i) {
        den <- sum((1 - as.numeric(n_per_group)/n) * (Vi_bm[i,] - bias))
        if(den > 0) sum(diag(H_obs %*% R_mat[sample(1:n, replace=TRUE), sample(1:n, replace=TRUE)])) / den else NA
      })
      p_ba <- (sum(f2_ba >= F2_obs, na.rm = TRUE) + 1) / (n_perms + 1)
    }
  }

  out <- data.frame(Df = object$CCA$rank, SumOfSqs = object$CCA$tot.chi, F2 = F2_obs, "Pr(>F2)" = p_perm, check.names = FALSE)
  if (!is.na(p_boot)) out$`Pr(>F2.boot)` <- p_boot
  if (!is.na(p_ba)) out$`Pr(>F2.ba)` <- p_ba
  out <- rbind(Model = out, Residual = c(object$CA$rank, object$CA$tot.chi, rep(NA, ncol(out)-2)), Total = c(n-1, object$tot.chi, rep(NA, ncol(out)-2)))

  # Attach distribution attribute
  attr(out, "F.perm") <- F2_p

  class(out) <- c("F2_adonis2", "anova.cca", "anova", "data.frame")
  return(out)
}

#' @export
plot.F2_adonis2 <- function(x, ...) {
  perm_vals <- attr(x, "F.perm")
  if (is.null(perm_vals)) stop("No permutation data found.")

  obs_f2 <- x["Model", "F2"]
  p_val <- x["Model", "Pr(>F2)"]

  # Standard histogram of the distribution
  hist(perm_vals,
       main = "F2 Permutation Distribution",
       xlab = "F2 Statistic",
       col = "skyblue",
       border = "white",
       xlim = range(c(perm_vals, obs_f2)))

  # Vertical line for observed F2
  abline(v = obs_f2, col = "red", lwd = 2, lty = 2)

  # Get the max x-value of the plot to place the legend accurately
  plot_max <- max(c(perm_vals, obs_f2))

  # Legend pushed to the absolute right edge
  legend(x = plot_max,
         y = max(hist(perm_vals, plot=FALSE)$counts), # Puts it at the top
         xjust = 1,  # Right-justifies the box so it stays inside the margin
         legend = c(paste("Observed F2 =", round(obs_f2, 3)),
                    paste("p-value =", round(p_val, 4))),
         col = c("red", NA), lwd = c(2, NA), bty = "n",
         cex = 0.8) # Slightly smaller text to ensure it fits
}

compute_G <- function(dmat_sq) {
  n <- nrow(dmat_sq); A <- -0.5 * dmat_sq
  C <- diag(n) - (1/n) * matrix(1, n, n)
  return(C %*% A %*% C)
}

get_H <- function(rhs) {
  q <- qr(rhs); Q <- qr.Q(q)[, 1:q$rank, drop = FALSE]
  return(tcrossprod(Q))
}

calculate_F2_and_Vi <- function(dmat_sq, grp_ID, SS.A, n) {
  n_i <- table(grp_ID); VL <- rep(0, length(n_i)); names(VL) <- names(n_i)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (grp_ID[i] == grp_ID[j]) { g <- as.character(grp_ID[i]); VL[g] <- VL[g] + (dmat_sq[i, j] / (n_i[g] * (n_i[g] - 1))) }
    }
  }
  denom <- sum((1 - (as.numeric(n_i) / n)) * VL)
  return(list(F2 = SS.A / denom, Vi = VL))
}

#' @export
nobs.F2_adonis2 <- function(object, ...) nrow(object$dmat_squared)
