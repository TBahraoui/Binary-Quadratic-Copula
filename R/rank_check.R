#' @export
rank_check <- function(fit, tol = 1e-8, n_states = 500, n_mc = 1000) {
  d <- ncol(fit$Omega)
  supp <- which(abs(fit$Omega) > tol & upper.tri(fit$Omega), arr.ind = TRUE)
  if (nrow(supp) == 0) {
    return(list(full_rank = TRUE, rank = 0L, cols = 0L,
                singular_values = numeric(0)))
  }
  states <- matrix(sample(0:1, n_states * d, replace = TRUE), n_states, d)
  Z <- mvtnorm::rmvnorm(n_mc, sigma = fit$Rho)
  thr <- stats::qnorm(1 - fit$p_hat)
  J <- matrix(0, n_states, nrow(supp))
  eps <- 1e-3
  for (s in seq_len(nrow(supp))) {
    j <- supp[s, 1]; k <- supp[s, 2]
    Op <- fit$Omega; Om <- fit$Omega
    Op[j, k] <- Op[j, k] + eps; Op[k, j] <- Op[j, k]
    Om[j, k] <- Om[j, k] - eps; Om[k, j] <- Om[j, k]
    p_plus  <- .pmf_quad(Z, Op, thr, states)
    p_minus <- .pmf_quad(Z, Om, thr, states)
    J[, s] <- (p_plus - p_minus) / (2 * eps)
  }
  sv <- svd(J, nu = 0, nv = 0)$d
  full_rank <- length(sv) >= nrow(supp) && min(sv) > tol
  list(full_rank = full_rank, rank = sum(sv > tol),
       cols = nrow(supp), singular_values = sv)
}

#' @keywords internal
.pmf_quad <- function(Z, Omega, thr, states) {
  n <- nrow(Z); d <- ncol(Z)
  Y <- matrix(0, n, d)
  for (j in seq_len(d)) Y[, j] <- Omega[j, j] * Z[, j]^2
  for (j in seq_len(d - 1)) for (k in (j + 1):d) {
    if (abs(Omega[j, k]) > 1e-12) {
      cross <- 2 * Omega[j, k] * Z[, j] * Z[, k]
      Y[, j] <- Y[, j] + cross; Y[, k] <- Y[, k] + cross
    }
  }
  X <- 1L * (Y > matrix(thr, n, d, byrow = TRUE))
  vapply(seq_len(nrow(states)),
         function(s) mean(rowSums(X == matrix(states[s, ], n, d,
                                              byrow = TRUE)) == d),
         numeric(1))
}
