#' @export
binarize <- function(Y, method = c('median', 'quantile', 'mixture'),
                     threshold = 0.5) {
  method <- match.arg(method)
  Y <- as.matrix(Y)
  n <- nrow(Y); d <- ncol(Y)
  X <- matrix(0L, n, d)
  for (j in seq_len(d)) {
    y <- Y[, j]
    if (method == 'median') {
      X[, j] <- as.integer(y > stats::median(y))
    } else if (method == 'quantile') {
      X[, j] <- as.integer(y > stats::quantile(y, threshold))
    } else {
      X[, j] <- .binarize_mixture(y)
    }
  }
  colnames(X) <- colnames(Y)
  X
}

#' @keywords internal
.binarize_mixture <- function(y) {
  if (length(unique(y)) < 3) return(as.integer(y > stats::median(y)))
  mu <- stats::quantile(y, c(0.25, 0.75), names = FALSE)
  sdv <- rep(stats::sd(y), 2)
  pi_v <- c(0.5, 0.5)
  for (iter in 1:100) {
    d1 <- pi_v[1] * stats::dnorm(y, mu[1], sdv[1])
    d2 <- pi_v[2] * stats::dnorm(y, mu[2], sdv[2])
    w1 <- d1 / (d1 + d2 + 1e-300)
    pi_v <- c(mean(w1), mean(1 - w1))
    mu <- c(stats::weighted.mean(y, w1),
            stats::weighted.mean(y, 1 - w1))
    sdv <- c(sqrt(stats::weighted.mean((y - mu[1])^2, w1)),
             sqrt(stats::weighted.mean((y - mu[2])^2, 1 - w1)))
    if (abs(diff(pi_v))[1] < 1e-6) break
  }
  cutoff <- tryCatch({
    stats::uniroot(
      function(c0) pi_v[1] * stats::pnorm(c0, mu[1], sdv[1]) -
                    (1 - pi_v[1]) * (1 - stats::pnorm(c0, mu[2], sdv[2])),
      c(min(y), max(y)))$root
  }, error = function(e) stats::median(y))
  as.integer(y > cutoff)
}
