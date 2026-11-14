#' Marginal CDF F_j(y) for the binary quadratic copula (Prop. 3.5)
#'
#' F_j(y) = Phi((-b + sqrt(Delta)) / (2c)) - Phi((-b - sqrt(Delta)) / (2c))
#' Delta(y) = b^2 - 4c(a - y)
#'
#' @param y Evaluation point.
#' @param a Intercept.
#' @param b Linear coefficient.
#' @param c Quadratic curvature (must be > 0).
#' @return F_j(y), the marginal CDF value.
#' @keywords internal
F_marginal <- function(y, a, b, c) {
  y_min <- a - b^2 / (4 * c)
  if (length(y) > 1) return(vapply(y, F_marginal,
                                    numeric(1), a = a, b = b, c = c))
  if (y <= y_min) return(0)
  Delta <- b^2 - 4 * c * (a - y)
  if (Delta < 0) return(0)
  sq <- sqrt(Delta)
  stats::pnorm((-b + sq) / (2 * c)) - stats::pnorm((-b - sq) / (2 * c))
}

#' Invert F_j to get tau_j = F_j^{-1}(1 - p)
#'
#' @param p Marginal probability P(X_j = 1).
#' @param a Intercept.
#' @param b Linear coefficient.
#' @param c Quadratic curvature.
#' @return tau_j (scalar).
#' @keywords internal
 tau_from_p <- function(p, a, b, c) {
  y_min <- a - b^2 / (4 * c)
  target <- 1 - p
  lo <- y_min + 1e-10
  hi <- y_min + 1
  while (F_marginal(hi, a, b, c) < target) {
    hi <- y_min + 2 * (hi - y_min)
    if (hi > y_min + 1e12) break
  }
  stats::uniroot(
    function(y) F_marginal(y, a, b, c) - target,
    lower = lo, upper = hi, tol = 1e-10
  )$root
}
