#' @export
imhof_cdf <- function(x, lambda, n_grid = 2000, t_max = 50) {
  lambda <- lambda[abs(lambda) > 1e-12]
  if (length(lambda) == 0) return(as.numeric(x >= 0))
  x <- as.numeric(x); out <- numeric(length(x))
  tg <- seq(1e-6, t_max, length.out = n_grid)
  dg <- c(diff(tg), tg[length(tg)] - tg[length(tg) - 1])
  for (i in seq_along(x)) {
    integrand <- vapply(tg, function(t) {
      phi <- exp(-0.5 * sum(log(1 + 4 * t^2 * lambda^2)) -
                 1i * sum(atan(2 * t * lambda)) +
                 1i * t * sum(lambda))
      Im(phi * exp(-1i * t * x[i])) / t
    }, numeric(1))
    out[i] <- 0.5 - (1 / pi) * sum(integrand * dg)
  }
  out
}
