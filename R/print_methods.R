#' @export
print.quadcopula <- function(x, ...) {
  cat('Binary Quadratic Copula Fit\n===========================\n')
  d <- ncol(x$Omega)
  cat(sprintf('Dimension: %d x %d\n', d, d))
  n_int <- sum(x$Omega[upper.tri(x$Omega)] != 0)
  cat(sprintf('Non-zero interactions: %d\n', n_int))
  cat(sprintf('Sparsity: %.1f%%\n',
              100 * (1 - n_int / max(1, d*(d-1)/2))))
  cat(sprintf('Log-likelihood: %.2f\n', x$loglik))
  cat(sprintf('HBIC: %.2f\n', x$hbic))
  cat(sprintf('Convergence proportion: %.2f\n', x$convergence_proportion))
  cat(sprintf('Assumption A5: %s\n',
              if (x$a5_ok) 'SATISFIED' else 'VIOLATED'))
  invisible(x)
}

#' @export
summary.quadcopula <- function(object, threshold = 0.1, ...) {
  cat('\nTop interactions (|Omega| > ', threshold, '):\n', sep = '')
  d <- ncol(object$Omega)
  rows <- list()
  for (j in seq_len(d - 1)) for (k in (j + 1):d) {
    if (abs(object$Omega[j, k]) > threshold) {
      nm1 <- if (is.null(rownames(object$Omega))) j else rownames(object$Omega)[j]
      nm2 <- if (is.null(rownames(object$Omega))) k else rownames(object$Omega)[k]
      rows[[length(rows) + 1L]] <- data.frame(
        Var1 = nm1, Var2 = nm2,
        Omega = round(object$Omega[j, k], 4),
        Type = ifelse(object$Omega[j, k] > 0,
                      'Synergistic', 'Antagonistic'),
        stringsAsFactors = FALSE)
    }
  }
  if (length(rows) == 0) {
    cat('No interactions above threshold.\n')
  } else {
    out <- do.call(rbind, rows)
    print(out[order(-abs(out$Omega)), ], row.names = FALSE)
  }
  invisible(object)
}
