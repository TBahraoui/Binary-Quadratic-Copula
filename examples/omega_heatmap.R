# Heatmap example for the binary quadratic copula
# Run from the repository root with the package installed.
# Output: binary-quadratic-copula-omega-heatmap.png

set.seed(42)
sim <- quadcopula::generate_dgp(n = 500, d = 9, dgp_type = 1)
fit <- quadcopula::estimate(
  sim$X,
  lambda = 0.1,
  gamma = 0.5,
  n_restarts = 1,
  max_iter = 200,
  verbose = FALSE
)

# Use a shared symmetric color scale so truth and estimate are comparable.
limit <- max(abs(c(sim$Omega_true, fit$Omega)), na.rm = TRUE)
if (!is.finite(limit) || limit == 0) limit <- 1
colors <- grDevices::colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(101)
variable_names <- colnames(sim$X)
if (is.null(variable_names)) variable_names <- paste0("X", seq_len(ncol(sim$X)))

plot_omega <- function(omega, title, labels, zlim, colors) {
  d <- nrow(omega)
  graphics::image(
    x = seq_len(d), y = seq_len(d), z = t(omega[d:1, , drop = FALSE]),
    col = colors, zlim = c(-zlim, zlim), axes = FALSE,
    xlab = "Variable", ylab = "Variable", main = title,
    useRaster = TRUE
  )
  graphics::axis(1, at = seq_len(d), labels = labels, las = 2, cex.axis = 0.8)
  graphics::axis(2, at = seq_len(d), labels = rev(labels), las = 2, cex.axis = 0.8)
  graphics::box()
}

# True and estimated off-diagonal interaction matrices.
grDevices::png(
  "binary-quadratic-copula-omega-heatmap.png",
  width = 1400, height = 700, res = 140
)
old_par <- graphics::par(mfrow = c(1, 2), mar = c(7, 7, 4, 2), oma = c(0, 0, 2, 0))
plot_omega(sim$Omega_true, "True interaction matrix", variable_names, limit, colors)
plot_omega(fit$Omega, "Estimated interaction matrix", variable_names, limit, colors)
graphics::mtext(
  "Binary quadratic copula: true and estimated Ω",
  outer = TRUE, cex = 1.2, font = 2
)
graphics::par(old_par)
grDevices::dev.off()
