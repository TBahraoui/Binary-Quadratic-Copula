#' @export
check_existence <- function(fit, tol = 1e-10) {
  Omega <- fit$Omega; Rho <- fit$Rho
  Rho_inv <- tryCatch(solve(Rho), error = function(e) NULL)
  if (is.null(Rho_inv)) { warning('Rho singular'); return(invisible(FALSE)) }
  diff_mat <- (Rho_inv - Omega + t(Rho_inv - Omega)) / 2
  ev <- Re(eigen(diff_mat, symmetric = TRUE, only.values = TRUE)$values)
  if (min(ev) <= tol) { warning('A5 violated'); return(invisible(FALSE)) }
  if (any(diag(Omega) <= 0)) { warning('Non-positive diag Omega');
    return(invisible(FALSE)) }
  invisible(TRUE)
}
