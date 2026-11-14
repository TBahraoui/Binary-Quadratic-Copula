#' Single-start coordinate descent on Rho (Fix B).
#'
#' Soft-threshold on rho_prop with magnitude step * lambda * gamma.
#' NO extra multiplier. Convergence only on max|Rho - Rho_old|.
#' @keywords internal
estimate_single <- function(X, theta_init, lambda,
                            max_iter = 500, tol = 1e-5, step = 0.005) {
  n <- nrow(X); d <- ncol(X)
  Rho   <- theta_init$Rho
  p_hat <- theta_init$p_hat
  gamma <- theta_init$gamma
  Rho0  <- diag(1, d)
  thr   <- stats::qnorm(1 - p_hat)
  off_idx <- which(upper.tri(Rho), arr.ind = TRUE)
  iter <- 0L
  for (iter in seq_len(max_iter)) {
    Rho_old <- Rho
    for (s in seq_len(nrow(off_idx))) {
      j <- off_idx[s, 1]; k <- off_idx[s, 2]
      g <- .score_rho(X, thr, Rho[j, k], j, k) / n
      rho_prop <- Rho[j, k] + step * g
      rho_thr <- sign(rho_prop) * max(abs(rho_prop) - step * lambda * gamma, 0)
      v <- max(min(rho_thr, 0.99), -0.99)
      Rho[j, k] <- v; Rho[k, j] <- v
    }
    ev <- min(eigen(Rho, symmetric = TRUE, only.values = TRUE)$values)
    if (ev <= 0.02) Rho <- project_feasible_rho(Rho)
    if (max(abs(Rho - Rho_old)) < tol) break
  }
  Omega_hat <- (Rho - Rho0) / gamma
  pl <- composite_pl_gauss(X, thr, Rho)
  list(Rho = Rho, Omega = Omega_hat, p_hat = p_hat, gamma = gamma,
       loglik = pl, converged = (iter < max_iter), iterations = iter)
}

#' @keywords internal
.score_rho <- function(X, thr, rho_jk, j, k) {
  n <- nrow(X); tj <- thr[j]; tk <- thr[k]
  rho_eff <- max(min(rho_jk, 0.99), -0.99)
  Sig <- matrix(c(1, rho_eff, rho_eff, 1), 2)
  p11 <- as.numeric(mvtnorm::pmvnorm(lower=c(tj,tk),upper=c(Inf,Inf),sigma=Sig))
  p10 <- as.numeric(mvtnorm::pmvnorm(lower=c(tj,-Inf),upper=c(Inf,tk),sigma=Sig))
  p01 <- as.numeric(mvtnorm::pmvnorm(lower=c(-Inf,tk),upper=c(tj,Inf),sigma=Sig))
  p00 <- 1 - p11 - p10 - p01
  dp11 <-  .phi2(tj,  tk, rho_eff)
  dp10 <- -.phi2(tj, -tk, rho_eff)
  dp01 <- -.phi2(-tj, tk, rho_eff)
  dp00 <-  .phi2(-tj,-tk, rho_eff)
  n11 <- sum(X[,j]==1L & X[,k]==1L)
  n10 <- sum(X[,j]==1L & X[,k]==0L)
  n01 <- sum(X[,j]==0L & X[,k]==1L)
  n00 <- sum(X[,j]==0L & X[,k]==0L)
  n11*dp11/max(p11,1e-12) + n10*dp10/max(p10,1e-12) +
  n01*dp01/max(p01,1e-12) + n00*dp00/max(p00,1e-12)
}

#' @keywords internal
.phi2 <- function(u, v, r) {
  exp(-(u^2+v^2-2*r*u*v)/(2*(1-r^2))) / (2*pi*sqrt(1-r^2))
}

#' @keywords internal
composite_pl_gauss <- function(X, thr, Rho) {
  d <- ncol(X); pl <- 0
  for (j in seq_len(d-1)) {
    tj <- thr[j]
    for (k in (j+1):d) {
      tk <- thr[k]
      re <- max(min(Rho[j,k],0.99),-0.99)
      Sig <- matrix(c(1,re,re,1),2)
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

#' @keywords internal
project_feasible_rho <- function(R) {
  R <- (R+t(R))/2
  e <- eigen(R, symmetric=TRUE)
  e$values <- pmax(Re(e$values), 0.02)
  R2 <- e$vectors %*% diag(e$values) %*% t(e$vectors)
  stats::cov2cor((R2+t(R2))/2)
}

#' @keywords internal
soft_threshold <- function(x, lambda) sign(x)*pmax(abs(x)-lambda, 0)

#' @keywords internal
is_positive_definite <- function(M) {
  ev <- tryCatch(Re(eigen((M+t(M))/2, symmetric=TRUE, only.values=TRUE)$values),
                 error=function(e) NA)
  !is.na(ev[1]) && min(ev) > 1e-8
}

#' @keywords internal
project_feasible <- function(Omega, Rho) {
  dm <- (solve(Rho)-Omega + t(solve(Rho)-Omega))/2
  e <- eigen(dm, symmetric=TRUE)
  e$values <- pmax(Re(e$values), 1e-6)
  dp <- e$vectors %*% diag(e$values) %*% t(e$vectors)
  (solve(Rho)-dp + t(solve(Rho)-dp))/2
}
