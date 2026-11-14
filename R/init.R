#' @keywords internal
fit_gaussian_copula <- function(X, gamma = 0.5) {
  n <- nrow(X); d <- ncol(X)
  p_hat <- pmax(pmin(colMeans(X), 1 - 1/n), 1/n)
  Rho <- tryCatch(stats::cor(X), error = function(e) diag(1, d))
  diag(Rho) <- 1
  ev <- tryCatch(min(Re(eigen(Rho, symmetric = TRUE,
                              only.values = TRUE)$values)),
                 error = function(e) NA)
  if (!is.na(ev) && ev <= 0) {
    Rho <- (Rho + t(Rho)) / 2
    e <- eigen(Rho, symmetric = TRUE)
    e$values <- pmax(Re(e$values), 0.01)
    Rho <- stats::cov2cor(e$vectors %*% diag(e$values) %*% t(e$vectors))
  }
  list(Rho = Rho, Omega = matrix(0, d, d), p_hat = p_hat, gamma = gamma)
}

#' @keywords internal
perturb_init <- function(theta_init, scale = 0.001) {
  d <- ncol(theta_init$Rho)
  Rho_p <- theta_init$Rho
  for (j in seq_len(d - 1)) for (k in (j + 1):d) {
    if (abs(Rho_p[j, k]) > 0.05) {
      Rho_p[j, k] <- Rho_p[j, k] + stats::rnorm(1, 0, scale)
      Rho_p[k, j] <- Rho_p[j, k]
    }
  }
  Rho_p <- (Rho_p + t(Rho_p)) / 2; diag(Rho_p) <- 1
  ev <- min(eigen(Rho_p, symmetric = TRUE, only.values = TRUE)$values)
  if (ev <= 0.02) Rho_p <- theta_init$Rho
  list(Rho = Rho_p, Omega = theta_init$Omega,
       p_hat = theta_init$p_hat, gamma = theta_init$gamma)
}
