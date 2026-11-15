# Arabidopsis thaliana isoprenoid-pathway network example
# Input: CSV with samples in rows, gene IDs in columns.
# Usage: Rscript examples/arabidopsis_isoprenoid_network.R expression.csv
# Supply identifiers appropriate to your annotation in pathway_genes below.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript examples/arabidopsis_isoprenoid_network.R expression.csv")
}
if (!requireNamespace("quadcopula", quietly = TRUE)) {
  stop("Install the quadcopula package before running this example.")
}
if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Install the optional igraph package: install.packages('igraph')")
}

# Replace with Arabidopsis gene IDs that match the expression CSV column names.
pathway_genes <- c("REPLACE_WITH_ARABIDOPSIS_GENE_ID")
expr <- as.matrix(read.csv(args[[1]], row.names = 1, check.names = FALSE))
if (!is.numeric(expr)) stop("Expression matrix must contain numeric values.")

matched_genes <- intersect(pathway_genes, colnames(expr))
if (length(matched_genes) < 2L || identical(matched_genes, "REPLACE_WITH_ARABIDOPSIS_GENE_ID")) {
  stop("Set pathway_genes to at least two valid gene IDs present in the input columns.")
}
pathway_expr <- expr[, matched_genes, drop = FALSE]
keep <- apply(pathway_expr, 2, stats::sd, na.rm = TRUE) > 0
pathway_expr <- pathway_expr[, keep, drop = FALSE]
if (ncol(pathway_expr) < 2L) stop("At least two non-constant pathway genes are required.")
if (anyNA(pathway_expr)) stop("Expression matrix contains missing values; clean them before fitting.")

X <- quadcopula::binarize(pathway_expr, method = "quantile", threshold = 0.5)
fit <- quadcopula::estimate(X, lambda = 0.5, verbose = TRUE)
png("arabidopsis-isoprenoid-network.png", width = 1000, height = 800)
on.exit(dev.off(), add = TRUE)
quadcopula::plot_network(
  fit,
  threshold = 0.05,
  main = "Arabidopsis thaliana isoprenoid-pathway associations"
)
