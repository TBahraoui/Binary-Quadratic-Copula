#' @export
bootstrap_ci <- function(X, fit, n_bootstrap = 200, level = 0.95,
                         n_restarts = 3, verbose = TRUE) {
  X <- as.matrix(X)
  n <- nrow(X); d <- ncol(X)
  Omega_boot <- array(NA_real_, dim = c(d, d, n_bootstrap))
  for (b in seq_len(n_bootstrap)) {
    if (verbose && b %% 50 == 0)
      cat(sprintf('Bootstrap %d/%d\n', b, n_bootstrap))
    idx <- sample.int(n, replace = TRUE)
    Xb <- X[idx, , drop = FALSE]
    if (all(colMeans(Xb) %in% c(0, 1))) next
    fit_b <- tryCatch(
      estimate(Xb, n_restarts = n_restarts, verbose = FALSE),
      error = function(e) NULL)
    if (!is.null(fit_b)) Omega_boot[, , b] <- fit_b$Omega
  }
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  ci <- array(NA_real_, dim = c(d, d, 2))
  for (j in seq_len(d - 1)) for (k in (j + 1):d) {
    vals <- Omega_boot[j, k, ]
    if (sum(!is.na(vals)) < 10) next
    ci[j, k, ] <- stats::quantile(vals, probs = probs, na.rm = TRUE)
    ci[k, j, ] <- ci[j, k, ]
  }
  list(ci = ci, Omega = fit$Omega, level = level)
}
