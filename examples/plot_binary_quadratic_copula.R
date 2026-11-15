# Example plot: simulated model's true quadratic interaction network.
# Requires the quadcopula package and the optional igraph package.

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Install the optional igraph package first: install.packages('igraph')")
}

set.seed(42)
sim <- quadcopula::generate_dgp(n = 500, d = 9, dgp_type = 1)

png("binary-quadratic-copula-network.png", width = 900, height = 700)
on.exit(dev.off(), add = TRUE)
quadcopula::plot_network(
  list(Omega = sim$Omega_true),
  threshold = 0.05,
  main = "Simulated binary quadratic copula interactions"
)
