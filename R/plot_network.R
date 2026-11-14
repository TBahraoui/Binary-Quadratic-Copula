#' @export
plot_network <- function(fit, threshold = 0.1,
                         main = 'Quadratic Copula Network',
                         node_color = NULL, ...) {
  if (!requireNamespace('igraph', quietly = TRUE))
    stop("Package 'igraph' is required.")
  d <- ncol(fit$Omega)
  edges <- which(abs(fit$Omega) > threshold & upper.tri(fit$Omega),
                 arr.ind = TRUE)
  if (nrow(edges) == 0) {
    warning('No edges above threshold.'); return(invisible(NULL))
  }
  vnames <- if (!is.null(rownames(fit$Omega))) rownames(fit$Omega) else
            as.character(seq_len(d))
  edge_df <- data.frame(
    from   = vnames[edges[, 1]],
    to     = vnames[edges[, 2]],
    weight = fit$Omega[edges],
    stringsAsFactors = FALSE)
  g <- igraph::graph_from_data_frame(edge_df, directed = FALSE)
  igraph::V(g)$color <- if (is.null(node_color)) 'steelblue' else node_color
  igraph::V(g)$size  <- 8
  igraph::E(g)$color <- ifelse(edge_df$weight > 0,
                               grDevices::rgb(0.2, 0.6, 0.2, 0.6),
                               grDevices::rgb(0.8, 0.2, 0.2, 0.6))
  igraph::E(g)$width <- pmax(abs(edge_df$weight) * 5, 0.5)
  plot(g, main = main, vertex.label.cex = 0.7, ...)
  invisible(g)
}
