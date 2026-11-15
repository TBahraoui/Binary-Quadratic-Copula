# Binary quadratic copula interaction heatmap

The SVG is a schematic illustration of a symmetric interaction matrix `Omega`. The runnable R script creates a PNG with side-by-side heatmaps of the simulated true interaction matrix and the matrix estimated from the simulated binary observations.

From the repository root, with the package installed, run:

```sh
Rscript examples/omega_heatmap.R
```

The script uses `generate_dgp()` and `estimate()` and saves `binary-quadratic-copula-omega-heatmap.png` in the current working directory. Color is scaled symmetrically around zero: red indicates positive values, blue negative values, and near-white values near zero. The diagonal is zero. Since the model simulation uses random draws, estimated values and the resulting heatmap may vary by R/package version.
