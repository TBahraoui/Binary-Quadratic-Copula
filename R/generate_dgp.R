#' Generate data from the binary quadratic copula (Fix B).
#' @export
generate_dgp <- function(n, d, dgp_type = 1, rho = 0.5, omega = 0.3,
                         c_val = 0.5, gamma = 0.5,
                         p_j = NULL, a_j = NULL, b_j = NULL) {
  stopifnot(dgp_type %in% 1:5, d >= 2, n >= 1)
  if (is.null(p_j)) p_j <- rep(0.5, d)
  if (is.null(a_j)) a_j <- rep(0, d)
  if (is.null(b_j)) b_j <- rep(1, d)
  if (dgp_type == 4) {
    Rho0 <- diag(d)
  } else if (dgp_type == 3) {
    Rho0 <- rho^abs(outer(seq_len(d), seq_len(d), '-'))
  } else {
    Rho0 <- matrix(rho, d, d); diag(Rho0) <- 1
  }
  Rho0 <- (Rho0 + t(Rho0)) / 2
  Omega <- diag(c_val, d)
  if (dgp_type == 1) {
    nb <- 3; bs <- max(d %/% nb, 1)
    for (bb in 0:(nb - 1)) {
      s <- bb * bs + 1; e <- min((bb + 1) * bs, d)
      if (s < e) for (j in s:(e - 1)) for (k in (j + 1):e) {
        v <- omega * stats::runif(1, 0.5, 1.5)
        Omega[j, k] <- Omega[k, j] <- v
      }
    }
  } else if (dgp_type == 2) {
    hubs <- sample(seq_len(d), min(3, d))
    for (h in hubs) {
      conn <- sample(setdiff(seq_len(d), hubs), min(10, max(1, d - 3)))
      for (j in conn) {
        v <- omega * stats::runif(1, 0.5, 1.5)
        Omega[h, j] <- Omega[j, h] <- v
      }
    }
  } else if (dgp_type == 3) {
    for (j in seq_len(d - 1)) for (k in (j + 1):min(j + 3, d)) {
      v <- omega * stats::runif(1, 0.5, 1.5)
      Omega[j, k] <- Omega[k, j] <- v
    }
  } else if (dgp_type == 4) {
    for (j in seq_len(d - 1)) for (k in (j + 1):d) {
      v <- omega * stats::runif(1, -0.5, 0.5) / d
      Omega[j, k] <- Omega[k, j] <- v
    }
  }
  Off <- Omega; diag(Off) <- 0
  R <- Rho0 + gamma * Off
  diag(R) <- 1
  R <- (R + t(R)) / 2
  ev <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
  shrink <- 1
  while (ev <= 0.02) {
    shrink <- shrink * 0.8
    R <- Rho0 + shrink * gamma * Off
    diag(R) <- 1
    R <- (R + t(R)) / 2
    ev <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
  }
  Rho_eff <- stats::cov2cor(R)
  Z <- mvtnorm::rmvnorm(n, sigma = Rho_eff)
  Y <- matrix(0, n, d)
  for (j in seq_len(d)) Y[, j] <- a_j[j] + b_j[j] * Z[, j]
  tau_true <- stats::qnorm(1 - p_j)
  X <- matrix(0L, n, d)
  for (j in seq_len(d)) X[, j] <- as.integer(Y[, j] > tau_true[j])
  colnames(X) <- paste0('X', seq_len(d))
  colnames(Z) <- paste0('Z', seq_len(d))
  list(X = X, Y = Y, Z = Z,
       Rho_true = Rho_eff, Omega_true = Omega,
       a_true = a_j, b_true = b_j, c_true = rep(c_val, d),
       tau_true = tau_true, p_true = p_j, gamma = gamma)
}
