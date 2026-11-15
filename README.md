# Binary-Quadratic-Copula

R code for Binary Quadratic Copula models, estimation, inference, and simulation examples.

## Example: visualize a binary quadratic copula model

The graphic below illustrates a quadratic interaction structure (`Omega`) for a simulated 9-variable model. Nodes are binary variables; edges represent pairwise quadratic interactions. This is a schematic example, not a fit to observed data.

![Example binary quadratic copula interaction network](examples/binary-quadratic-copula-network.svg)

```r
# Install igraph if needed: install.packages("igraph")
set.seed(42)
sim <- quadcopula::generate_dgp(n = 500, d = 9, dgp_type = 1)
png("binary-quadratic-copula-network.png", width = 900, height = 700)
quadcopula::plot_network(
  list(Omega = sim$Omega_true), threshold = 0.05,
  main = "Simulated binary quadratic copula interactions"
)
dev.off()
```

## Applied-data workflow examples

The scripts below show how to prepare expression matrices for two use cases: Arabidopsis thaliana isoprenoid-pathway genes and single-cell RNA-seq. They expect your own expression data; no Arabidopsis or scRNA-seq dataset is bundled in this repository. The resulting networks are exploratory estimated pairwise associations, not causal pathway edges. Gene identifiers must match the column names in your input matrix.

### Arabidopsis thaliana isoprenoid pathway

Prepare a CSV with samples in rows, genes in columns, and gene IDs as column names. Add the pathway gene identifiers appropriate for your annotation to the script, then run:

```sh
Rscript examples/arabidopsis_isoprenoid_network.R expression.csv
```

Edit `pathway_genes` in the script to list the Arabidopsis genes in your pathway data. The script keeps matching, variable genes, binarizes their expression, estimates the model, and saves an interaction-network PNG.

### Single-cell RNA-seq

Provide a normalized expression matrix saved as an RDS file, with cells in rows and genes in columns. The script chooses a small set of variable genes to keep the example computationally manageable, binarizes expression above each gene's 75th percentile, estimates the model, and saves a network plot:

```sh
Rscript examples/scrna_network.R normalized_expression.rds
```

Both workflows require the package and optional `igraph` dependency. See the scripts for input checks and output filenames.
