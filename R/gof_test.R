#' @export
gof_test <- function(X, n_bootstrap = 1000, alpha = 0.05, n_restarts = 5) {
  X <- as.matrix(X)
  n <- nrow(X); d <- ncol(X)
  p_hat <- pmax(pmin(colMeans(X), 1 - 1/n), 1/n)
  thr <- stats::qnorm(1 - p_hat)
  Rho0 <- tryCatch(stats::cor(X), error = function(e) diag(1, d))
  diag(Rho0) <- 1
  ev <- min(eigen(Rho0, symmetric = TRUE, only.values = TRUE)$values)
  if (ev <= 0) { Rho0 <- Rho0 + diag(abs(ev) + 0.01, d)
                 Rho0 <- stats::cov2cor(Rho0) }
  ll0 <- .loglik_null(X, thr, Rho0)
  fit_full <- tryCatch(
    estimate(X, n_restarts = n_restarts, verbose = FALSE),
    error = function(e) list(loglik = ll0))
  T_obs <- 2 * (fit_full$loglik - ll0)
  T_boot <- numeric(n_bootstrap)
  for (b in seq_len(n_bootstrap)) {
    if (b %% 200 == 0) cat(sprintf('Bootstrap %d/%d\n', b, n_bootstrap))
    Zb <- mvtnorm::rmvnorm(n, sigma = Rho0)
    Xb <- matrix(0L, n, d)
    for (j in seq_len(d)) Xb[, j] <- as.integer(Zb[, j] > thr[j])
    ll0b <- .loglik_null(Xb, thr, Rho0)
    fitb <- tryCatch(
      estimate(Xb, n_restarts = n_restarts, verbose = FALSE),
      error = function(e) list(loglik = ll0b))
    T_boot[b] <- 2 * (fitb$loglik - ll0b)
  }
  p_value <- mean(T_boot >= T_obs)
  list(statistic = T_obs, p_value = p_value,
       bootstrap_dist = T_boot, reject = p_value < alpha)
}

#' @keywords internal
.loglik_null <- function(X, thr, Rho) {
  n <- nrow(X); d <- ncol(X); pl <- 0
  for (j in seq_len(d - 1)) {
    tj <- thr[j]
    for (k in (j + 1):d) {
      tk <- thr[k]
      Sig <- matrix(c(1, Rho[j, k], Rho[j, k], 1), 2)
      p11 <- as.numeric(mvtnorm::pmvnorm(lower=c(tj,tk),upper=c(Inf,Inf),sigma=Sig))
      p10 <- as.numeric(mvtnorm::pmvnorm(lower=c(tj,-Inf),upper=c(Inf,tk),sigma=Sig))
      p01 <- as.numeric(mvtnorm::pmvnorm(lower=c(-Inf,tk),upper=c(tj,Inf),sigma=Sig))
      p00 <- 1 - p11 - p10 - p01
      n11 <- sum(X[,j]==1L & X[,k]==1L)
      n10 <- sum(X[,j]==1L & X[,k]==0L)
      n01 <- sum(X[,j]==0L & X[,k]==1L)
      n00 <- sum(X[,j]==0L & X[,k]==0L)
      pl <- pl + n11*log(max(p11,1e-12)) + n10*log(max(p10,1e-12)) +
                 n01*log(max(p01,1e-12)) + n00*log(max(p00,1e-12))
    }
  }
  pl
}
