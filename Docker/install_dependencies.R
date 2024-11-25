# Instalar remotes si no está instalado
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cran.rstudio.com/")
}

# Paquetes y versiones específicas
packages <- list(
  "rgl" = "1.2.8",
  "ucie" = "1.0.2",
  "chameleon" = "0.2-3",
  "aricode" = "1.0.2",
  "ggplot2" = "3.4.3",
  "dplyr" = "1.1.3",
  "pracma" = "2.4.2",
  "doParallel" = "1.0.17",
  "iterators" = "1.0.14",
  "foreach" = "1.5.2",
  "igraph" = "1.5.1"
)

for (pkg in names(packages)) {
  remotes::install_version(pkg, version = packages[[pkg]], repos = "https://cran.rstudio.com/")
}

# Instalar ComplexHeatmap

# Install ComplexHeatmap version 2.14.0 from tarball
tarball_url <- "https://mghp.osn.xsede.org/bir190004-bucket01/archive.bioconductor.org/packages/3.16/bioc/src/contrib/ComplexHeatmap_2.14.0.tar.gz"
destfile <- tempfile(fileext = ".tar.gz")

# Download the tarball
download.file(tarball_url, destfile)

# Install the package from the tarball
install.packages(destfile, repos = NULL, type = "source")

# Clean up temporary tarball file
unlink(destfile)
