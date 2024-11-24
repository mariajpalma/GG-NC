# Instalar remotes si no está instalado
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cran.rstudio.com/")
}

# Paquetes y versiones específicas
packages <- list(
  "rgl" = "1.2.8",
  "ucie" = "1.0.2",
  "chameleon" = "0.2-3",
  "ComplexHeatmap" = "2.14.0",
  "aricode" = "1.0.2",
  "ggplot2" = "3.4.3",
  "dplyr" = "1.1.3",
  "pracma" = "2.4.2",
  "doParallel" = "1.0.17",
  "iterators" = "1.0.14",
  "foreach" = "1.5.2",
  "igraph" = "1.5.1"
)

# Instalar versiones específicas de CRAN
for (pkg in names(packages)) {
  remotes::install_version(pkg, version = packages[[pkg]], repos = "https://cran.rstudio.com/")
}

# Paquetes de Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cran.rstudio.com/")
}
BiocManager::install("ComplexHeatmap", version = "2.14.0", update = FALSE)
