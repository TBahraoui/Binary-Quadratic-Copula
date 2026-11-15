# Single-cell RNA-seq network example
# Input: normalized numeric matrix saved as RDS, cells in rows and genes in columns.
# Usage: Rscript examples/scrna_network.R normalized_expression.rds

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript examples/scrna_network.R normalized_expression.rds")
}
if (!requireNamespace("quadcopula", quietly = TRUE)) {
  stop("Install the quadcopula package before running this example.")
}
if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Install the optional igraph package: install.packages('igraph')")
}

expr <- as.matrix(readRDS(args[[1]]))
if (!is.numeric(expr)) stop("Expression matrix must contain numeric values.")
if (nrow(expr) < 2L || ncol(expr) < 2L) stop("Matrix must have at least two cells and two genes.")
if (anyNA(expr)) stop("Expression matrix contains missing values; clean them before fitting.")

# Select a small variable-gene panel for a tractable, illustrative fit.
gene_variance <- apply(expr, 2, stats::var)
variable_genes <- names(sort(gene_variance[is.finite(gene_variance)], decreasing = TRUE))
selected_genes <- head(variable_genes, 30L)
if (length(selected_genes) < 2L) stop("At least two variable genes are required.")
sc_expr <- expr[, selected_genes, drop = FALSE]

# A high quantile makes the binary state represent relatively high expression.
X <- quadcopula::binarize(sc_expr, method = "quantile", threshold = 0.75)
keep <- apply(X, 2, function(x) length(unique(x)) > 1L)
X <- X[, keep, drop = FALSE]
if (ncol(X) < 2L) stop("Fewer than two genes vary after binarization; adjust preprocessing/threshold.")

fit <- quadcopula::estimate(X, lambda = 0.5, verbose = TRUE)
png("scrna-binary-quadratic-network.png", width = 1000, height = 800)
on.exit(dev.off(), add = TRUE)
quadcopula::plot_network(
  fit,
  threshold = 0.05,
  main = "Single-cell RNA-seq binary quadratic associations"
)
