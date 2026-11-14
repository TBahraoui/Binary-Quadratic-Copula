#' @export
cross_validate <- function(X, K = 10, lambda = 0.5,
                           n_restarts = 5, n_repeats = 1, seed = 42) {
  X <- as.matrix(X)
  n <- nrow(X)
  set.seed(seed)
  cv_pl <- numeric(n_repeats)
  fold_pl_all <- vector('list', n_repeats)
  for (r in seq_len(n_repeats)) {
    folds <- sample(rep(seq_len(K), length.out = n))
    fold_pl <- numeric(K)
    for (k in seq_len(K)) {
      X_tr <- X[folds != k, , drop = FALSE]
      X_te <- X[folds == k, , drop = FALSE]
      fit <- tryCatch(
        estimate(X_tr, n_restarts = n_restarts, lambda = lambda,
                 verbose = FALSE),
        error = function(e) NULL)
      if (is.null(fit)) { fold_pl[k] <- NA; next }
      thr <- stats::qnorm(1 - fit$p_hat)
      fold_pl[k] <- composite_pl_gauss(X_te, thr, fit$Rho)
    }
    fold_pl_all[[r]] <- fold_pl
    cv_pl[r] <- mean(fold_pl, na.rm = TRUE)
  }
  list(cv_pl = cv_pl, fold_pl = fold_pl_all)
}
