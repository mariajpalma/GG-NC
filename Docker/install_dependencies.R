# Instalar remotes si no está instalado
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cran.rstudio.com/")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cran.rstudio.com/")
}

# Paquetes y versiones específicas
packages <- list(
  "rgl" = "1.3.18",
  "ucie" = "1.0.2",
  "chameleon" = "0.2-3",
  "aricode" = "1.0.3",
  "ggplot2" = "3.5.2",
  "dplyr" = "1.1.4",
  "pracma" = "2.4.4",
  "doParallel" = "1.0.17",
  "iterators" = "1.0.14",
  "foreach" = "1.5.2",
  "igraph" = "2.1.4",
  "viridis" = "0.6.5",
  "viridisLite" = "0.4.2",
  "circlize" = "0.4.16",
  "data.table" = "1.18.2.1",
  "pryr" = "0.1.6"
  #"ComplexHeatmap" = "2.22.0"
)

for (pkg in names(packages)) {
  remotes::install_version(pkg, version = packages[[pkg]], repos = "https://cran.rstudio.com/")
}

# Instalar ComplexHeatmap desde Bioconductor
BiocManager::install("ComplexHeatmap", ask = FALSE, update = FALSE)

# Instalar ComplexHeatmap

# Install ComplexHeatmap version 2.14.0 from tarball
# <- "https://mghp.osn.xsede.org/bir190004-bucket01/archive.bioconductor.org/packages/3.16/bioc/src/contrib/ComplexHeatmap_2.14.0.tar.gz"
#destfile <- tempfile(fileext = ".tar.gz")

# Download the tarball
#download.file(tarball_url, destfile)

# Install the package from the tarball
#install.packages(destfile, repos = NULL, type = "source")

# Clean up temporary tarball file
#unlink(destfile)
