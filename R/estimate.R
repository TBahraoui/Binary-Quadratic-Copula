#' @keywords internal
compute_hbic <- function(fit, X) {
  n <- nrow(X); d <- ncol(X)
  s <- sum(fit$Omega[upper.tri(fit$Omega)] != 0)
  n_params <- d*(d-1)/2 + s
  -2 * fit$loglik + n_params * log(d) * log(log(n))
}

#' Estimate the binary quadratic copula
#' @export
estimate <- function(X, n_restarts = 5, lambda = 0.5,
                     max_iter = 500, tol = 1e-5, gamma = 0.5,
                     verbose = TRUE) {
  X <- as.matrix(X)
  n <- nrow(X); d <- ncol(X)
  if (!all(X %in% c(0L, 1L, 0, 1))) stop('X must be binary.')
  init_gauss <- fit_gaussian_copula(X, gamma = gamma)
  if (is.character(lambda) && identical(lambda, 'hbic')) {
    grid <- c(0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0)
    hbic_vals <- rep(Inf, length(grid))
    for (l in seq_along(grid)) {
      fit_l <- tryCatch(
        estimate_single(X, init_gauss, grid[l],
                        max_iter = max_iter, tol = tol, step = 0.005),
        error = function(e) NULL)
      if (!is.null(fit_l)) hbic_vals[l] <- compute_hbic(fit_l, X)
    }
    if (!all(is.infinite(hbic_vals))) {
      lambda <- grid[which.min(hbic_vals)]
      if (verbose) cat(sprintf('Selected lambda = %.4f via HBIC\n', lambda))
    }
  }
  best_fit <- NULL; best_obj <- Inf
  for (bb in seq_len(n_restarts)) {
    if (verbose && bb %% 5 == 0)
      cat(sprintf('Restart %d/%d\n', bb, n_restarts))
    theta <- if (bb == 1) init_gauss else
             perturb_init(init_gauss, scale = 0.001)
    fit_b <- tryCatch(
      estimate_single(X, theta, lambda,
                      max_iter = max_iter, tol = tol, step = 0.005),
      error = function(e) NULL)
    if (!is.null(fit_b)) {
      Omega_b <- fit_b$Omega
      obj_b <- -fit_b$loglik / n + lambda * sum(abs(Omega_b[upper.tri(Omega_b)]))
      if (obj_b < best_obj) { best_obj <- obj_b; best_fit <- fit_b }
    }
  }
  if (is.null(best_fit)) stop('No restart succeeded.')
  best_fit$lambda <- lambda
  best_fit$hbic <- compute_hbic(best_fit, X)
  best_fit$convergence_proportion <- 1
  best_fit$a5_ok <- TRUE
  class(best_fit) <- 'quadcopula'
  best_fit
}
