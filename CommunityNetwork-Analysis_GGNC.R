#PipelineV251023
###############################################################################
################################# LIBRARIES ###################################
###############################################################################

library(igraph)
library(foreach)
library(doParallel)
library(pracma)
library(dplyr)
library(ggplot2)
library(aricode)
library(ComplexHeatmap)
library(chameleon)
library(ucie)
library(rgl)

###############################################################################
################################# FUNCTIONS ###################################
###############################################################################
# This function processes IBD (Identity by Descent) data to identify
# related individuals, hubs, and unrelated individuals based on a given threshold.
# It then constructs a network graph and returns the igraph object along with
# colors for genetic regions.

ibd_network <- function (path, fname1, fname2, max, prune=FALSE) {
  #' Create an IBD Network from IBD Data and Additional Information
  #'
  #' This function reads IBD (Identity-by-Descent) data from a file, processes it to identify related individuals,
  #' filters out hubs (highly connected individuals), and constructs a network using the remaining individuals. 
  #' It can also prune the resulting network based on the input parameters.
  #'
  #' @param path A character string specifying the path to the directory containing the input files.
  #' @param fname1 A character string specifying the name of the file containing the IBD data.
  #' @param fname2 A character string specifying the name of the file containing additional information (e.g., genetic regions).
  #' @param max A numeric value specifying the IBD threshold for selecting related individuals.
  #' @param prune A logical value indicating whether to prune the network (default is FALSE).
  #'
  #' @details
  #' The function performs the following steps:
  #' \itemize{
  #'   \item Reads the IBD data file and additional information file.
  #'   \item Identifies related individuals based on the specified IBD threshold (`max`).
  #'   \item Identifies highly connected individuals (hubs) and removes them.
  #'   \item Constructs an igraph network from the remaining unrelated individuals.
  #'   \item Optionally prunes the network if `prune = TRUE`.
  #' }
  #' The function also writes the list of hub individuals and unrelated individuals to separate text files.
  #'
  #' @return A list containing:
  #' \itemize{
  #'   \item \code{ibd_graph}: An igraph object representing the IBD network.
  #'   \item \code{spop_color}: A data frame mapping genetic regions to their corresponding colors.
  #' }
  #'
  #' @export
  #'
  #' @examples
  #' \dontrun{
  #' path <- "data/"
  #' fname1 <- "ibd_data.txt"
  #' fname2 <- "info_data.csv"
  #' ibd_network(path, fname1, fname2, max = 7, prune = TRUE)
  #' }
  #'
  # Check if files exist
  if (!file.exists(file.path(path, fname1))) {
    stop("The IBD file does not exist at the specified path.")
  }
  if (!file.exists(file.path(path, fname2))) {
    stop("The additional information file does not exist at the specified path.")
  }
  # Read the IBD file
  ibd <- read.table(file.path(path, fname1), header = T)
  # Select related individuals based on IBD threshold
  related <- ibd[which(ibd$IBD_CM_SUM >= max),]
  # Extract IDs of related individuals
  relatives_1 <- related$ID1
  relatives_2 <- related$ID2
  relatives <- c(relatives_1, relatives_2)
  # Calculate frequency of each individual
  freq <- table(relatives)
  # Identify hubs (individuals with IBD connections above threshold)
  hubs <- names(freq[freq > 1])
  # Write hubs to a text file for information
  directorio = file.path(path, "Removed_samples")
  if (!file.exists(directorio)) {
    # If it does not exist, create the directory
    dir.create(directorio)
    cat("'Removed_samples' folder created successfully.\n")
  }
  write.table(
    hubs,
    file = file.path(directorio, "Hubs_IBD.txt"),
    sep = ",",
    row.names = F,
    col.names = F
  )
  # Remove hubs from IBD data
  ibd <- ibd[!(ibd$ID1 %in% hubs | ibd$ID2 %in% hubs), ]
  # Rest of the related individuals
  related_remaining <- ibd[ibd$IBD_CM_SUM >= max, ]
  Nohub <- unique(related_remaining$ID1)
  # Write unrelated individuals to a text file
  write.table(
    Nohub,
    file = file.path(directorio, "NoHubs_IBD.txt"),
    sep = ",",
    row.names = F,
    col.names = F
  )
  # Remove unrelated individuals from IBD data
  ibd <- ibd[!(ibd$ID1 %in% Nohub),]
  # Convert IBD data to data frame format
  unrelated <- ibd
  print(dim(unrelated))
  unrelated_df <- unrelated[, c("ID1", "ID2", "IBD_CM_SUM")]
  colnames(unrelated_df) <- c("ID1", "ID2", "weight")
  # Extract unique IDs
  id1 <- cbind(unrelated$ID1, unrelated$ID2)
  id <- c(id1[, 1], id1[, 2])
  id <- unique(id)
  # Read additional information from second file
  info <-
    read.csv(file.path(path, fname2), sep = "\t", head = TRUE)
  info = info[, c(1, 5, 8)]
  # Filter information based on IDs
  filtrado <- info[info[, 1] %in% id, ]
  id_order <- match(filtrado[, 1], id)
  id_ordenado <- filtrado[order(id_order),]
  colnames(id_ordenado) <- c("name", "genetic_region", "color")
  # Extract unique genetic regions and corresponding colors
  info_nodup <-
    info[!duplicated(info$genetic_region, fromLast = FALSE), ]
  spop_color <- cbind(info_nodup$genetic_region, info_nodup$color)
  colnames(spop_color) <- c("genetic_region", "color")
  # Create an igraph from the unrelated data
  ibd_graph <-
    graph_from_data_frame(unrelated_df, directed = F, vertices = id_ordenado)
  # Check if it is an igraph object
  if (is.igraph(ibd_graph)) {
    print("The igraph object has been created.")
  } else {
    stop("The igraph object has NOT been created.")
  }
  # Set colors of vertices in the graph
  V(ibd_graph)$color <- V(ibd_graph)$color 
  V(ibd_graph)$spop <- V(ibd_graph)$genetic_region
  # Set widths of edges based on weight
  E(ibd_graph)$width <-
    (E(ibd_graph)$weight / max(E(ibd_graph)$weight)) * 50
  # Prune the network if requested
  if (prune)
  {
    ibd_graph <- prune_network(ibd_graph, path)
  }
  # Check if the graph has vertices
  if (vcount(ibd_graph) == 0) {
    stop("The graph has no vertices.")
  }
  # Check if the graph has edges
  if (ecount(ibd_graph) == 0) {
    stop("The graph has no edges.")
  }
  # Return the graph and the colors by superpopulation per genetic region
  return(list(ibd_graph, spop_color))
}

prune_network <- function (net, path, random_seed=NULL){
  #' Prune an IBD Network by Removing Low-Degree and Small Community Nodes
  #'
  #' This function prunes an IBD (Identity-by-Descent) network by iteratively removing vertices (nodes) with a degree of 1 and isolated nodes. 
  #' It then performs community detection using the Louvain method and removes communities with fewer than 25 members.
  #' This function is usefult to get the networks shown in the paper.
  #' @param net An igraph object representing the IBD network to be pruned.
  #' @param path A character string specifying the path to the directory where files with the removed samples (nodes) will be saved.
  #'
  #' @details
  #' The function performs the following steps:
  #' \itemize{
  #'   \item Iteratively removes vertices with a degree of 1 and isolated vertices (degree 0).
  #'   \item Detects communities using the Louvain community detection algorithm.
  #'   \item Identifies and removes communities with fewer than 25 members.
  #'   \item Writes the lists of removed vertices (degree 0, degree 1, and small communities) to text files in the specified path.
  #' }
  #'
  #' @return An igraph object representing the pruned IBD network.
  #'
  #' @export
  #'
  #' @examples
  #' \dontrun{
  #' path <- "data/"
  #' pruned_net <- prune_network(ibd_net, path, random_seed)
  #' }
  #'
  # While there are nodes with a degree of 1
  degree1 <- c()
  degree0 <- c()
  while (sum(degree(net) == 1)> 0){
    # Remove vertices with degree of 1
    degree1 <- c(degree1, names(which(degree(net) == 1)))
    net <- delete_vertices(net, names(which(degree(net) == 1)))
    # Also remove any isolated vertices
    degree0 <- c(degree0, names(which(degree(net) == 0)))
    net <- delete_vertices(net, names(which(degree(net) == 0)))
  }
  # Perform community detection using the Louvain method
  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }
  clusters <- cluster_louvain(net, resolution = 0.01)
  # Identify communities with fewer than 25 members
  to_remove <- names(which(table(clusters$membership)< 25))
  # Initialize an empty list to store the individuals to be removed
  ind_list <- c()
  # Iterate over the communities that need to be removed
  for (i in to_remove){
    # Get the indices of nodes in the current community
    comm <- which(clusters$membership==i)
    # Add the names of the nodes in the community to the removal list
    ind_list <- c(ind_list,clusters$names[comm])
    # Remove the nodes of the current community from the network
    net <- delete_vertices(net, clusters$names[comm])
  }
  write.table(
    degree0,
    file = file.path(path, "/Removed_samples/degree0.txt"),
    sep = ",",
    row.names = F,
    col.names = F
  )
  write.table(
    degree1,
    file = file.path(path, "/Removed_samples/degree1.txt"),
    sep = ",",
    row.names = F,
    col.names = F
  )
  
  write.table(
    ind_list,
    file = file.path(path, "/Removed_samples/communities.txt"),
    sep = ",",
    row.names = F,
    col.names = F
  )
  
  # Return the pruned network
  return(net)
}



# This function reads PCA data and additional information about individuals (such as population labels) from files.
# It then calculates a correlation matrix, creates an adjacency matrix, and constructs a network graph based on
# a specified threshold. Population labels and colors are assigned to vertices in the graph, and the function returns
# the igraph object along with information about genetic regions and their colors.
pca_network <- function (path, fname1, fname2, max) {
  #' Create a PCA Network from PCA Data and Additional Information
  #'
  #' This function constructs a network from Principal Component Analysis (PCA) data, where the nodes represent individuals, and the edges represent correlations between individuals' PCA scores. It also assigns colors to nodes based on genetic regions.
  #'
  #' @param path A character string specifying the path to the directory containing the input files.
  #' @param fname1 A character string specifying the name of the file containing the PCA data.
  #' @param fname2 A character string specifying the name of the file containing additional information (e.g., genetic regions).
  #' @param max A numeric value specifying the threshold for edge weights (correlation) below which edges will be removed.
  #'
  #' @details
  #' The function performs the following steps:
  #' \itemize{
  #'   \item Reads the PCA data and additional information file.
  #'   \item Constructs a correlation matrix from the PCA scores.
  #'   \item Creates an undirected network graph where edges are weighted by the correlation between individuals' PCA scores.
  #'   \item Removes edges with weights below the specified threshold (`max`).
  #'   \item Assigns colors to nodes based on the individuals' genetic regions as specified in the additional information file.
  #' }
  #'
  #' @return A list containing:
  #' \itemize{
  #'   \item \code{pca_graph}: An igraph object representing the PCA network.
  #'   \item \code{spop_color}: A data frame mapping genetic regions to their corresponding colors.
  #' }
  #'
  #' @export
  #'
  #' @examples
  #' \dontrun{
  #' path <- "data/"
  #' fname1 <- "pca_data.txt"
  #' fname2 <- "info_data.csv"
  #' pca_network(path, fname1, fname2, max = 0.7)
  #' }
  #'
  # Check if files exist
  if (!file.exists(file.path(path, fname1))) {
    stop("The PCA file does not exist at the specified path.")
  }
  if (!file.exists(file.path(path, fname2))) {
    stop("The additional information file does not exist at the specified path.")
  }
  # Read the PCA data
  PCS <- read.table(file.path(path, fname1), header = T)
  # Read additional information
  info <- read.csv(file.path(path, fname2), sep = "\t", head = TRUE)
  # Remove individuals from PC
  PCS_noidv = PCS
  rownames(PCS_noidv) <- PCS_noidv[, 1]
  PCS_noidv <- PCS_noidv[, -1]
  # Calculate correlation matrix
  cor_mat <- cor(t(PCS_noidv))
  # Create an adjacency matrix
  net_pca <-
    graph.adjacency(
      adjmatrix = cor_mat,
      weighted = T,
      diag = F,
      mode = "upper"
    )
  #An undirected graph will be created, only the upper right triangle (including the diagonal)
  # is used (for the edge weights).
  # Remove edges with weight below the threshold
  pca_graph <- delete.edges(net_pca, which(E(net_pca)$weight < max))
  # Assign colors based on genetic regions
  indv_spop_color <- info[, c(1, 5, 8)]
  vert <- matrix(0, nrow = length(V(net_pca)$name), ncol = 1)
  vert[, 1] <- V(net_pca)$name
  # Filter and order individuals
  filtrado <- indv_spop_color[indv_spop_color[, 1] %in% vert[, 1],]
  id_order <- match(filtrado[, 1], vert[, 1])
  id_ordenado <- filtrado[order(id_order), ]
  # Colors by genetic region
  info_nodup <-
    info[!duplicated(info$genetic_region, fromLast = FALSE), ]
  spop_color <- cbind(info_nodup$genetic_region, info_nodup$color)
  colnames(spop_color) <- c("genetic_region", "color")
  # Add genetic region and color to the igraph object
  V(pca_graph)$spop <- id_ordenado[, 2]
  V(pca_graph)$color <- id_ordenado[, 3]
  # Return the graph and the colors by superpopulation per genetic region
  # Check if it is an igraph object
  if (is.igraph(pca_graph)) {
    print("The igraph object has been created.")
  } else {
    stop("The igraph object has NOT been created.")
  }
  # Check if the graph has vertices
  if (vcount(pca_graph) == 0) {
    stop("The graph has no vertices.")
  }
  # Check if the graph has edges
  if (ecount(pca_graph) == 0) {
    stop("The graph has no edges.")
  }
  return(list(pca_graph, spop_color))
}


# The function ReadGRMBin was taken from https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM
# R script to read the GRM binary file
ReadGRMBin = function(prefix, AllN = F, size = 4) {
  # Helper function to compute the sum of integers from 1 to i
  sum_i = function(i) {
    return(sum(1:i))
  }
  # Define file names based on the given prefix
  BinFileName = paste(prefix, ".grm.bin", sep = "")
  NFileName = paste(prefix, ".grm.N.bin", sep = "")
  IDFileName = paste(prefix, ".grm.id", sep = "")
  # Read IDs from the ID file
  id = read.table(IDFileName)
  n = dim(id)[1]
  # Open binary GRM file and read data
  BinFile = file(BinFileName, "rb")
  grm = readBin(BinFile,
                n = n * (n + 1) / 2,
                what = numeric(0),
                size = size)
  # Open binary N file and read data
  NFile = file(NFileName, "rb")
  if (AllN == T) {
    N = readBin(NFile,
                n = n * (n + 1) / 2,
                what = numeric(0),
                size = size)
  }
  else
    N = readBin(NFile,
                n = 1,
                what = numeric(0),
                size = size)
  # Compute indices for diagonal elements
  i = sapply(1:n, sum_i)
  # Return the GRM data as a list
  return(list(
    diag = grm[i],
    off = grm[-i],
    id = id,
    N = N
  ))
}



# Function to create a genetic relationship matrix (GRM) network
grm_network <- function (path, fname1, fname2, max) {
  #' Create a Genetic Relationship Matrix (GRM) Network
  #'
  #' This function constructs a network from a Genetic Relationship Matrix (GRM) binary file. It generates a symmetric adjacency matrix from the GRM data, where nodes represent individuals, and edges represent genetic relationships. The edges are weighted by the genetic relationship values, and edges below a certain threshold are removed. The function also assigns colors to nodes based on the genetic regions specified in an additional information file.
  #'
  #' @param path A character string specifying the path to the directory containing the GRM binary files and additional information file.
  #' @param fname1 A character string specifying the base name of the GRM binary files (excluding extensions).
  #' @param fname2 A character string specifying the name of the file containing additional information (e.g., genetic regions).
  #' @param max A numeric value specifying the threshold for edge weights (genetic relationships) below which edges will be removed.
  #'
  #' @details
  #' The function performs the following steps:
  #' \itemize{
  #'   \item Checks if the necessary GRM binary files (.grm.bin, .grm.N.bin, .grm.id) exist in the specified path.
  #'   \item Reads the GRM binary files and constructs a symmetric adjacency matrix from the diagonal and upper triangle elements of the GRM.
  #'   \item Creates an undirected igraph network, where edges are weighted by the genetic relationship values.
  #'   \item Removes edges with weights below the specified threshold (`max`).
  #'   \item Reads additional information and assigns colors to nodes based on genetic regions.
  #'   \item Returns the resulting igraph network and a data frame mapping genetic regions to their corresponding colors.
  #' }
  #'
  #' @return A list containing:
  #' \itemize{
  #'   \item \code{net_grm}: An igraph object representing the GRM network.
  #'   \item \code{spop_color}: A data frame mapping genetic regions to their corresponding colors.
  #' }
  #'
  #' @export
  #'
  #' @examples
  #' \dontrun{
  #' path <- "data/"
  #' fname1 <- "grm_data"
  #' fname2 <- "info_data.csv"
  #' grm_network(path, fname1, fname2, max = 0.2)
  #' }
  #'
  
  # Check if files exist
  datafile = file.path(path, fname1)
  if (!file.exists(paste(datafile, ".grm.bin", sep = ""))) {
    stop(
      "The GRM file ",
      paste(datafile, ".grm.bin", sep = ""),
      "does not exist at the specified path."
    )
  }
  if (!file.exists(paste(datafile, ".grm.N.bin", sep = ""))) {
    stop(
      "The GRM file ",
      paste(datafile, ".grm.N.bin", sep = ""),
      "does not exist at the specified path."
    )
  }
  if (!file.exists(paste(datafile, ".grm.id", sep = ""))) {
    stop(
      "The GRM file ",
      paste(datafile, ".grm.id", sep = ""),
      "does not exist at the specified path."
    )
  }
  # Read the GRM binary files
  GRM = ReadGRMBin(datafile)
  # Extract the diagonal and upper triangle of the matrix
  diagonal <- GRM[[1]]  # Diagonal elements of the GRM
  triangulo_superior <-
    GRM[[2]]  # Upper triangle elements of the GRM
  id <- GRM[[3]] # Read the GRM binary files
  id <- id[, 2] # Extract individual IDs
  
  # Calculate the number of rows and columns in the matrix
  n <- length(diagonal)
  matriz_simetrica <- matrix(0, n, n)
  
  # Llenar la matriz con los valores de la diagonal
  matriz_simetrica[upper.tri(matriz_simetrica)] <-
    triangulo_superior
  matriz_simetrica <- matriz_simetrica + t(matriz_simetrica)
  # Fill the symmetric matrix with the GRM values
  rownames(matriz_simetrica) <- GRM$id[, 2]
  colnames(matriz_simetrica) <- GRM$id[, 2]
  # Create an igraph from the adjacency matrix
  net_grm <-
    graph.adjacency(
      adjmatrix = matriz_simetrica,
      weighted = T,
      diag = F,
      mode = "lower"
    )
  # Remove edges with weights below the threshold
  net_grm <- delete.edges(net_grm, which(E(net_grm)$weight < max))
  # Read additional information for individuals
  infofile <- file.path(path, fname2)
  info <- read.csv(file.path(path, fname2), sep = "\t", head = TRUE)
  indv_spop_color <- info[, c(1, 5, 8)]
  vert <- names(V(net_grm))
  # Match and order the individual information
  filtrado <- indv_spop_color[indv_spop_color[, 1] %in% vert,]
  id_order <- match(filtrado[, 1], vert)
  id_ordenado <- filtrado[order(id_order), ]
  # Create a color mapping for superpopulations
  info_nodup <-
    info[!duplicated(info$genetic_region, fromLast = FALSE), ]
  spop_color <- cbind(info_nodup$genetic_region, info_nodup$color)
  colnames(spop_color) <- c("genetic_region", "color")
  # Add superpopulation names and colors as vertex attributes
  V(net_grm)$spop <- id_ordenado[, 2]
  V(net_grm)$color <- id_ordenado[, 3]
  # Check if it is an igraph object
  if (is.igraph(net_grm)) {
    print("The igraph object has been created.")
  } else {
    stop("The igraph object has NOT been created.")
  }
  # Check if the graph has vertices
  if (vcount(net_grm) == 0) {
    stop("The graph has no vertices.")
  }
  # Check if the graph has edges
  if (ecount(net_grm) == 0) {
    stop("The graph has no edges.")
  }
  # Returns the igraph object and the the genetic regions and their corresponding colors.
  return(list(net_grm, spop_color))
}


# This function makes an igraph object depending on the kind and values of the input data.
input <- function(kind, path, data, info, max, prune) {
  #' Create a Network Based on Input Data Type
  #'
  #' This function creates an igraph object from different types of genetic data: Identical by Descent (IBD), Principal Component Analysis (PCA), or Genetic Relationship Matrix (GRM). It selects the appropriate network creation function depending on the specified data type and parameters.
  #'
  #' @param kind A character string specifying the type of data to create the network from. Acceptable values are \code{"IBD"}, \code{"PCA"}, and \code{"GRM"}.
  #' @param path A character string specifying the path to the directory containing the input files.
  #' @param data A character string specifying the file name for the data (depending on the \code{kind} of network).
  #' @param info A character string specifying the file name for the additional information (e.g., genetic regions).
  #' @param max A numeric value specifying the threshold for edge weights (depending on the type of data) below which edges will be removed.
  #' @param prune A logical value indicating whether or not to prune the network (used only for IBD networks).
  #'
  #' @details
  #' The function performs the following steps:
  #' \itemize{
  #'   \item Determines the type of network to create based on the \code{kind} parameter.
  #'   \item Calls the appropriate network creation function: \code{ibd_network}, \code{pca_network}, or \code{grm_network}.
  #'   \item If the type of data is not one of the accepted options (\code{"IBD"}, \code{"PCA"}, \code{"GRM"}), it displays an error message.
  #' }
  #'
  #' @return A list containing the network (igraph object) and a data frame mapping genetic regions to their corresponding colors, or \code{NULL} if an invalid data type is provided.
  #'
  #' @export
  #'
  #' @examples
  #' \dontrun{
  #' path <- "data/"
  #' data <- "ibd_data.txt"
  #' info <- "info_data.csv"
  #' max <- 0.5
  #' prune <- TRUE
  #' network <- input("IBD", path, data, info, max, prune)
  #' }
  #'
  # Determine the type of network to create
  if (kind == "IBD") {
    network <-
      ibd_network(path, data, info, max, prune)
  } else{
    if (kind == "PCA") {
      network <-
        pca_network(path, data, info, max)
    } else{
      if (kind == "GRM") {
        network <-
          grm_network(path, data, info, max)
      } else{
        cat(
          "The permissible options for the data type are limited to:",
          "\n",
          "'IBD' (Identical by descent), 'PCA'(Principal Component Analysis) or 'GRM' (Genetic Relationship Matrix).\n"
        )
      }
    }
  }
  # Returns a list with the igraph object and the genetic regions and their corresponding colors.
  return(network)
}



# Function to perform clustering over a range of resolutions and save the results
get_lambda_results <- function(input, steps, path, name, lower_limit, upper_limit, use_leiden=FALSE, random_seed) {
  #' Perform Clustering Over a Range of Resolutions and Save Results
  #'
  #' This function performs clustering using the Louvain method across a range of resolutions and saves the resulting community memberships for each resolution to a specified file.
  #'
  #' @param input An igraph object representing the network on which clustering will be performed.
  #' @param steps An integer indicating the number of resolution steps to evaluate.
  #' @param path A character string specifying the path to the directory where the results will be saved.
  #' @param name A character string specifying the file name for saving the clustering results.
  #' @param lower_limit A numeric value indicating the lower bound of the resolution range for clustering.
  #' @param upper_limit A numeric value indicating the upper bound of the resolution range for clustering.
  #' @param use_leiden Logical. If TRUE, the Leiden algorithm is used; if FALSE, the Louvain algorithm is used.
  #' @param random_seed A numeric value to use as seed.
  #' @details
  #' The function iterates over a range of resolution values (generated using a logarithmic scale between \code{lower_limit} and \code{upper_limit}) and performs community detection clustering at each resolution. The community memberships for each resolution are saved to a file, and the number of communities for each resolution is recorded.
  #'
  #' @return A list containing the clustering results for each resolution, where each entry is a result of \code{cluster_louvain} or \code{cluster_leiden} applied with a different resolution.
  #'
  #' @export
  #'
  #' @examples
  #' \dontrun{
  #' net <- igraph::make_ring(10)
  #' path <- "results/"
  #' name <- "clustering_results.txt"
  #' lower_limit <- 0.1
  #' upper_limit <- 1
  #' steps <- 5
  #' results <- get_lambda_results(net, steps, path, name, lower_limit, upper_limit)
  #' }
  #'  
  # An empty list to store the clustering results for each resolution.
  list_cl <- list()
  # A numeric vector to store the number of communities for each resolution.
  ncom <- 0
  a <- 1
  # Iterates over a sequence of resolution parameters
  print(random_seed)
  for (i in logspace(lower_limit, upper_limit, n = steps)) {
    # Use seed
    if (!is.null(random_seed)) {
      set.seed(random_seed)
    }
    # The community detection  function is called with the input data and the current resolution to perform clustering.
    if (use_leiden) {
      cl <- cluster_leiden(input, resolution_parameter = i)
    } else {
      cl <- cluster_louvain(input, resolution = i)
    }
    # The clustering result cl is stored in list_cl with the resolution value as the key.
    list_cl[[as.character(i)]] <- cl
    # The membership of the clustering result is transposed and written to the specified file.
    write.table(
      t(cl$membership),
      sep = ",",
      file = file.path(path, name),
      append = TRUE,
      quote = F,
      row.names = F,
      col.names = F
    )
    # The number of communities (ncom[a]) is calculated and stored.
    ncom[a] <- max(cl$membership)
    a <- a + 1
  }
  if (length(list_cl) > 0 && !any(is.na(list_cl))) {
    print("The grouping into various resolutions has been successfully completed.")
  } else {
    stop("The list of groups is either empty or contains NA.")
  }
  # The function returns the list list_cl, which contains the clustering results for each resolution.
  return(list_cl)
}


# nodebased, relabel and pollock are a set of functions that were developed by Dr. Lewis for the research paper titled
# "The function of communities in protein interaction networks at multiple scales."
# The functions have been adapted from Matlab to R by our research group and some modifications were done 
# to avoid label the same communities with overlap 0 and to color comunities that are > min at least once. 
# For more detailed information, please refer to the publication by Lewis, Anna CF, et al., 2010.
# We are available to assist and address any inquiries related to these functions.

# This function `nodebased(Sloop)` is used to re-label communities in a way that preserves their identity.
# In such a way that we can know how individuals are accommodating and reaccommodating in the different communities throughout lambda.
# The function takes a matrix `Sloop`, which contains community assignments for each node at different parameter values
# (assumed to be the highest resolution parameter) in the rows.

nodebased <- function(Sloop) {
  #This R function, `nodebased(Sloop)`, is used to re-label communities while ensuring their continuity.
  #This re-labeling process allows us to ascertain how individuals are adapting and readapting within various communities across lambda.
  #The function requires a matrix, `Sloop`, which includes community assignments for each node at varying parameter values.
  #It is assumed that these parameter values represent the highest resolution.
  N <- nrow(Sloop)
  Ssort <- matrix(0, nrow = N, ncol = ncol(Sloop))
  S1 <- Sloop[N,]
  Ssort[N,] <- S1
  #The loop provided iterates from `i = 1` to `N - 1` in order to reassign labels to Ssort, starting from the N-1th row (N-i)
  # and moving towards the first row in reverse order.
  #The loop will continue until all parameter values (resolutions) have been processed.
  # As a result, `Ssort` will hold the re-labeled communities for each parameter value.
  cat("Relabeling communities. Step:\n")
  for (i in 1:(N - 1)) {
    cat(i, "\n")
    S2 <- Sloop[N - i,]
    currentmax <- max(Ssort)
    Snew <- relabel(i, S1, S2, currentmax)
    Ssort[N - i,] <- Snew
    S1 <- Snew
  }
  if (nrow(Ssort) > 0 && ncol(Ssort) > 0 && any(!is.na(Ssort))) {
    print("Communities have been relabeled.")
  } else {
    stop("The Ssort matrix is empty or contains NA.")
  }
  return(Ssort)
}
# This `relabel` function takes several arguments: `num`, `S1`, `S2`, and `currentmax`.
# The purpose of this function is to re-label communities in `S2` based on the mapping provided by `S1`,
# while ensuring that the re-labeled communities retain their identity.
relabel <- function(num, S1, S2, currentmax) {
  # Function for pairwise re-labeling
  # The purpose of the `relabel` function is to re-label communities in `S2` while preserving their identity.
  # This function takes several arguments, including `num`, `S1`, `S2`, and `currentmax`.
  # The communities in `S2` are relabeled based on the mapping provided by `S1`.
  
  # Create a matrix `Snew` with 1 row and the same number of columns as `S1`.
  #This will be used to store the re-labeled communities. It is filled with zeros.
  Snew <- matrix(0, nrow = 1, ncol = length(S1))
  # Calculate the number of unique communities in `S2` and store it in `nocoms2`.
  nocoms2 <- length(unique(S2))
  # Get the unique community labels in `S2`, sort them, and store them in `coms2`. Convert `coms2` to a vector.
  coms2 <- sort(unique(S2))
  # Calculate the number of unique communities in `S1` and store it in `nocoms1`.
  nocoms1 <- length(unique(S1))
  # Get the unique community labels in `S1`, sort them, and store them in `coms1`. Convert `coms1` to a vector.
  coms1 <- sort(unique(S1))
  # Initialize a variable `n` to 0. This will be used to keep track of the iteration count for filling matrix `l`.
  n <- 0
  # Create a matrix `C` filled with zeros, with dimensions `nocoms1` by `nocoms2`.
  # This matrix will store the overlap values between communities in `S1` and `S2`.
  C <- matrix(0, nrow = nocoms1, ncol = nocoms2)
  # Create an empty matrix `l` with dimensions `nocoms1 * nocoms2` by 3.
  # This matrix will be used to store the overlap values, corresponding community labels in `S1`,
  # and corresponding community labels in `S2`.
  l <- matrix(0, nrow = nocoms1 * nocoms2, ncol = 3)
  # Create two matrices `A` and `B` to store the community memberships of each node in `S1` and `S2`,
  # respectively, using a one-hot encoding scheme.
  A <- matrix(0, nrow = length(S1), ncol = nocoms1)
  for (i in 1:nocoms1)
    A[, i] <- (as.numeric(S1 == as.numeric(coms1[i])))
  B <- matrix(0, nrow = length(S2), ncol = nocoms2)
  for (i in 1:nocoms2)
    B[, i] <- (as.numeric(S2 == as.numeric(coms2[i])))
  # Calculate the overlap values between communities in `S1` and `S2`, and store them in matrix `C`.
  # Fill the matrix `l` with the overlap values and corresponding community labels.
  for (i in 1:nocoms1) {
    for (j in 1:nocoms2) {
      n <- n + 1
      # this bit particular to nodes based
      # Find the overlap of nodes in one of the communities in partition 1
      # with nodes in one of the communities in partition 2, for all pairs
      # of communitites
      C[i, j] <- sum(A[, i] & B[, j]) / sum(A[, i] | B[, j])
      l[n, 1] <- C[i, j]
      l[n, 2] <- as.numeric(coms1[i])
      l[n, 3] <- as.numeric(coms2[j])
    }
  }
  # Sort the rows of `l` based on the overlap values in descending order, and store the indices of the sorted rows in `ind`.
  val <- l[, 1]
  ind <- order(val, decreasing = TRUE)
  # Create an empty matrix `l2` to store the sorted rows of `l`.
  l2 <-
    l[ind,] # this has highest overlap in top row, etc. First column is value of overlap, second is label you want, third is what it was called in the partition you want to relabel
  # Initialize a variable `m` to 0. This will be used as a counter to iterate through the rows of `l2`.
  m <- 0
  assigned <-
    matrix(0, nrow = 1, ncol = 2) # just added to make coding easier, ignored later
  noverlap <- 0
  noverlap_coms <- c()
  # Enter a `while` loop to handle cases where there are more communities in `S2` than rows in `l2`.
  # Handle cases with more communities in S2
  # If there are more communities in `S2`, add new labels to the communities without labels in `assigned`,
  # where `assigned` is a matrix that stores assigned community labels.
  while (nrow(assigned) < (nocoms2 + 1)) {
    m = m + 1
    if (is.matrix(l2)) {
      filas_l2 = nrow(l2)
    } else{
      if (is.vector(l2)) {
        filas_l2 = 1
      } else{
        cat("An unexpected error occurred in relabel.\n")
      }
    }
    if (m > filas_l2) {
      currentmax = currentmax + noverlap
      unassigned = sort(setdiff(coms2, assigned[, 2])) # communities without labels
      for (k in 1:length(unassigned)) {
        assigned = rbind(assigned, c(currentmax + k, unassigned[k]))
      }
    } else {
      # If `m` is not greater than the number of rows in `l2`, proceed with the re-labeling process.
      if (is.matrix(l2)) {
        consider = l2[m, ]
      } else{
        if (is.vector(l2)) {
          consider = l2
        } else{
          cat("There is an error when trying to assing communities.\n")
        }
      }
      if ((length(intersect(consider[3], assigned[, 2])) == 0) &&
          (length(intersect(consider[2], assigned[, 1])) == 0)) {
        if (consider[1] == 0) {
          noverlap = noverlap + 1
          noverlap_coms <- c(noverlap_coms, consider[3])
          assigned = rbind(assigned, c(currentmax + noverlap, consider[3]))
        }
        else{
          assigned = rbind(assigned, consider[2:3])
        }
      }
    }
  }
  
  ## Now actually relabel
  for (j in 2:nrow(assigned)) {
    x = which(S2 == assigned[j, 2])
    Snew[x] = assigned[j, 1]
  }
  return(Snew)
}

# This `pollock` function takes three arguments: `Ssort`, `little`, and `r`.
# The purpose of this function is to reorder the node labels in `Ssort` to make sense of the data for plotting.
# It also assigns colors to communities based on their sizes, and provides a visualization of the communities with specified colors.
pollock <- function(Ssort, little, r, path, lower_limit, upper_limit ) {
  # to make pretty for plot. 
  # r is which layer (counting from bottom) will be perfectly sorted 
  # (choose on what looks prettiest. Try 0 as default.)
  
  # This code reorders the node labels. New labels are in 'order'. It only
  # colours communities which reach, for some value of the parameter, size
  # 'little'. These coms are 'bigcoms'.
  
  # First, reorder nodes to make sense of what is going on. NB there will be
  # no way to keep all communities together over all values of parameter in
  # 1D! Usually best to ensure that largest communitites are fully ordered,
  # having made sure that smaller communities are fully ordered first.
  
  # Calculate the number of rows in `Ssort` and store it in `N`.
  N <- nrow(Ssort)
  # Create an initial ordering of node labels in `order`, ranging from 1 to the number of columns in `Ssort`.
  order <- 1:ncol(Ssort)
  # Reorder the columns of `Ssort` and update the `order` based on the values of the `r`-th row of `Ssort`.
  for (l in 1:(N - r)) {
    # r will be the perfectly ordered partition
    ix <- order(order = Ssort[N + 1 - l,])
    Ssort <- Ssort[, ix]
    order <- order[ix]
  }
  # Now arrange to colour all communities that are never over size 'little' to
  # be painted a uniform colour, here, white
  # Create a vector `bigcoms` to store the community labels of communities that have a size greater than `little`.
  bigcoms <- integer(0)
  # Find the maximum value in `Ssort` and store it in `m`.
  m <- max(Ssort)
  # Create a matrix `sizes` filled with zeros, with dimensions `N` by `m`. This matrix will store the sizes of each community in `Ssort`.
  sizes <- matrix(0, nrow = N, ncol = m)
  # Calculate the sizes of each community in `Ssort` and store them in the `sizes` matrix.
  for (i in 1:N) {
    b <- sort(unique(Ssort[i,]))
    for (j in 1:length(b)) {
      sizes[i, b[j]] <- sum(Ssort[i,] == b[j])
    }
  }
  # Determine the community labels of communities with sizes greater than `little`, and store them in `bigcoms`.
  for (k in 1:m) {
    sizescom <- sizes[, k]
    if (sum(sizescom > little) >= 1) {
      bigcoms <- c(bigcoms, k)
    }
  }
  # Create a vector `smallcoms` that contains the community labels of communities with sizes less than or equal to `little`.
  smallcoms <- sort(setdiff(1:m, bigcoms))
  # Create a copy of `Ssort`, named `Ssort2`, and set the community labels in `smallcoms` to 0.
  Ssort2 <- Ssort
  for (i in 1:length(smallcoms)) {
    y <- which(Ssort == smallcoms[i])
    Ssort2[y] <- 0
  }
  
  # Rename coms, to take into account have lost all littlecoms
  # Create a vector `rep` that contains the unique community labels in `Ssort2`, sorted in ascending order.
  rep <-
    sort(unique(as.vector(Ssort2)))
  # Create a vector `a` that contains the integers from 1 to the length of `rep`.
  ### Ask if there are small comms
  if (rep[1] == 0) {
  a <- 1:length(rep)
  }else{
  a <- 2:(length(rep)+1)
  }
  im <- Ssort2
  # Replace the community labels in `Ssort2` with their corresponding indices in `a`.
  for (i in 1:length(rep)) {
    z <- which(Ssort2 == rep[i])
    im[z] <- a[i]
  }
  # Number of communities
  m <- length(unique(as.vector(im)))
  if (rep[1] == 0) {
    cat("The resolution matrix contains small communities that will appear in white.")
    # Select colors for communities.
    colores <- distinct_colors(
      (m - 1),
      minimal_saturation = 33,
      minimal_lightness = 20,
      maximal_lightness = 80
    )
    #### Select alternative colors1
    colores_alt1 <- distinct_colors(
      (m - 1),
      minimal_saturation = 33,
      minimal_lightness = 20,
      maximal_lightness = 80
    )
    colores_alt2 <- distinct_colors(
      (m - 1),
      minimal_saturation = 33,
      minimal_lightness = 20,
      maximal_lightness = 80
    )
    
    # Obtain colors in hex
    map <- colores$name
    colalt1 <- colores_alt1$name
    colalt2 <- colores_alt2$name
    # Mix the colors
    map <- map[sample(length(map))]
    colalt1 <- colalt1[sample(length(colalt1))]
    colalt2 <- colalt1[sample(length(colalt2))]
    
    # Store colors in vector ma
    ma <- vector("character", length = m)
    col_ma1 <- vector("character", length = m)
    col_ma2 <- vector("character", length = m)
    ma[2:m] <- map
    col_ma1[2:m] <- colalt1
    col_ma2[2:m] <- colalt2
    # Makes first one white. Little coms all have label 0, so will appear white.
    ma[1] <-
      c("#FFFFFF")
    col_ma1[1] <- c("#FFFFFF")
    col_ma2[1] <- c("#FFFFFF")
  }else{
    # Select colors for communities.
    colores <- distinct_colors(
      m,
      minimal_saturation = 33,
      minimal_lightness = 20,
      maximal_lightness = 80
    )
    #### Select alternative colors1
    colores_alt1 <- distinct_colors(
      m,
      minimal_saturation = 33,
      minimal_lightness = 20,
      maximal_lightness = 80
    )
    colores_alt2 <- distinct_colors(
      m,
      minimal_saturation = 33,
      minimal_lightness = 20,
      maximal_lightness = 80
    )
    
    # Obtain colors in hex
    map <- colores$name
    colalt1 <- colores_alt1$name
    colalt2 <- colores_alt2$name
    # Mix the colors
    map <- map[sample(length(map))]
    colalt1 <- colalt1[sample(length(colalt1))]
    colalt2 <- colalt1[sample(length(colalt2))]
    
    # Store colors in vector ma
    ma <- vector("character", length = m)
    col_ma1 <- vector("character", length = m)
    col_ma2 <- vector("character", length = m)
    ma[1:m] <- map
    col_ma1[1:m] <- colalt1
    col_ma2[1:m] <- colalt2
  }
  ### Save alternative colors into a file
  write(col_ma1, file=paste0(path, "/Browser_files/alternative_colors_1.txt"))
  write(col_ma2, file=paste0(path, "/Browser_files/alternative_colors_2.txt"))
  im2 <- apply(im, 2, rev)
  # Axis X values
  x_vals <- 1:ncol(Ssort)
  # Axis Y values
  y_vals <- 1:nrow(Ssort)
  # Store the plot as a .png image.
  png(file= paste(path,"/Community_detection/res_plot.png",sep="/"), res=300, width = 2000, height = 2000)
  image(
    x = x_vals,
    y = y_vals,
    z = t(im2),
    useRaster = TRUE,
    col = ma,
    xaxs = "i",
    yaxt = "n",
    xlab = "Individuals",
    ylab = " ",
    las = 1,
    cex.lab = 1.6,
    cex.axis= 1.6
    
  )
  mtext(
    side = 2,
    line = 2,
    text = "Resolution value",
    cex = 1.6,
  )
  # Y labels
  yticklabels <- rev(lower_limit:upper_limit)
  # Create the position vector for the y-axis ticks.
  yticks <- seq(1, nrow(Ssort), length.out = length(yticklabels))
  # Configure y-axis ticks and labels
  axis(
    side = 2,
    at = yticks,
    labels = yticklabels,
    las = 1,
    cex.axis = 1.6
  )
  dev.off()
  cat("Te resolution plot has been generated\n")
  # Returns a list with the individuals communities by level or resolution, the individal orders and the colors by community.
  return(list(im, order, ma))
}


# Function to plot community detection clustering results colored by superpopulation
plot_louvain_by_spop <- function(graph, cl_list, spop_color, R, lay, name) {
  #' Plot Community detection Results Colored by Superpopulation
  #'
  #' This function plots the Community detection results for a given resolution, with vertices colored by superpopulation, and saves the plot as an SVG file.
  #'
  #' @param graph An igraph object representing the network to be plotted.
  #' @param cl_list A list of Louvain/Leiden clustering results, where each entry corresponds to a different resolution.
  #' @param spop_color A data frame containing the superpopulation information and their corresponding colors. The first column should contain superpopulation names, and the second column should contain colors.
  #' @param R A numeric value representing the resolution for which to plot the Louvain clustering result.
  #' @param lay A matrix specifying the layout of the graph for the plot (typically generated with layout functions in igraph).
  #' @param name A character string specifying the file name (including path) for saving the SVG plot.
  #'
  #' @details
  #' The function extracts the clustering results for the specified resolution and plots the graph using the provided layout. Vertices are colored according to their superpopulation. The plot is saved as an SVG file, and a legend indicating the superpopulation colors is added to the plot.
  #'
  #' @return The clustering result (\code{cl}) for the specified resolution.
  #'
  #' @import igraph
  #' @export
  #'
  #' @examples
  #' \dontrun{
  #' net <- igraph::make_ring(10)
  #' cl_list <- list("1" = igraph::cluster_louvain(net))
  #' spop_color <- data.frame(genetic_region = c("pop1", "pop2"),
  #'                          color = c("red", "blue"))
  #' lay <- igraph::layout_with_fr(net)
  #' R <- 1
  #' name <- "louvain_plot.svg"
  #' plot_louvain_by_spop(net, cl_list, spop_color, R, lay, name)
  #' }
  #'
  # Extract the clustering results for the given resolution
  cl <- cl_list[[as.character(R)]]
  # Open an SVG device for plotting
  png(
    file = paste(path, name, sep = "/"),
    height = 2000,
    width = 2000,
    res=300
  )
  # Plot the graph with Louvain clustering results
  plot(
    cl,
    graph,
    layout = lay,
    vertex.size = 6,
    vertex.label = NA,
    edge.color = "gray",
    col = V(graph)$color
  )
  # Add a legend to the plot
  legend(
    1.1,
    1.1,
    legend = levels(as.factor(spop_color[, 1])) ,
    bty = "n",
    pch = 20,
    pt.cex = 3,
    cex = 3 ,
    horiz = FALSE,
    inset = c(0.1, 0.1),
    col = spop_color[, 2]
  )
  dev.off()
  # Return the clustering result
  return(cl)
}



#This function creates a heatmap that shows the overlap between communities and the defined populations and super populations.
overlap_hm <- function(path, fname2, cl_list, R, spop_color, name) {
  #' Heatmap of Overlap Between Communities and Populations
  #'
  #' This function creates a heatmap that visualizes the overlap between detected communities (from Louvain/Leiden clustering) and predefined populations or superpopulations. The heatmap is normalized by community size and can include row annotations based on superpopulation information.
  #'
  #' @param path A character string specifying the directory path where the input and output files are located.
  #' @param fname2 A character string specifying the file name of the additional information containing population and superpopulation data.
  #' @param cl_list A list of Louvain/Leiden clustering results for various resolutions.
  #' @param R A numeric value specifying the resolution at which to extract the clustering result.
  #' @param spop_color A data frame containing superpopulation names and their corresponding colors. The first column should contain the superpopulation names and the second column the associated colors.
  #' @param name A character string specifying the file name for saving the heatmap (including the path).
  #'
  #' @details
  #' The function reads in the population and superpopulation information, matches it with the Louvain clustering results, and calculates an overlap matrix that indicates the proportion of individuals in each population belonging to different communities. The overlap matrix is then normalized by the community size and visualized in a heatmap.
  #'
  #' The function uses \code{ComplexHeatmap} for plotting the heatmap and includes row annotations for superpopulations, which are colored according to the input color mapping.
  #'
  #' @return A heatmap is saved as a PNG file in the specified directory, and a confirmation message is printed to the console.
  #'
  #' @import ComplexHeatmap
  #' @export
  #'
  #' @examples
  #' \dontrun{
  #' cl_list <- list("1" = igraph::cluster_louvain(net))
  #' spop_color <- data.frame(genetic_region = c("pop1", "pop2"),
  #'                          color = c("red", "blue"))
  #' overlap_hm("output_path", "info_file.csv", cl_list, 1, spop_color, "overlap_heatmap.png")
  #' }
  #'
  # Check if files exist
  if (!file.exists(file.path(path, fname2))) {
    stop("The additional information file does not exist at the specified path.")
  }
  # Get the list of community names
  cl <- cl_list[[as.character(R)]]
  info <-
    read.csv(file.path(path, fname2), sep = "\t", head = TRUE)
  # Select relevant columns
  info = info[, c(1, 3, 5, 8)]
  # Filter the information to include only names in cl
  filtrado <- info[info[, 1] %in% cl$name, ]
  # Order the information according to the order in cl$name
  id_order <- match(filtrado[, 1], cl$name)
  id_ordenado <- filtrado[order(id_order),]
  # Combine the information with community assignments
  idv_comm_pop_spop <-
    cbind(id_ordenado, cl$membership)
  idv_comm_pop_spop <- idv_comm_pop_spop[, c(1, 2, 3, 5, 4)]
  # Calculate the community overlap table
  tabla <-
    xtabs( ~ idv_comm_pop_spop[, 4] + idv_comm_pop_spop[, 2], data = idv_comm_pop_spop)
  # Normalize the overlap table
  tabla_overlap <- tabla
  hil = 1:nrow(tabla_overlap)
  sum_by_comm = sapply(hil,  function(hil) {
    sum(tabla_overlap[hil, ])
  })
  for (hil in 1:nrow(tabla_overlap)) {
    for (col in 1:ncol(tabla_overlap)) {
      tabla_overlap[hil, col] = tabla[hil, col] / sum_by_comm[hil]
    }
  }
  # Assign row names
  rownames(tabla_overlap) <- paste0("C", 1:nrow(tabla_overlap))
  # Calculate the number of population levels
  num_2level <- length(unique(idv_comm_pop_spop[, 3]))
  # Calculate the population count by subpopulation
  list_pop_by_spop <-
    aggregate(
      idv_comm_pop_spop[, 1] ~ idv_comm_pop_spop[, 3] + idv_comm_pop_spop[, 2],
      data = idv_comm_pop_spop,
      FUN = length
    )
  # Create row annotation for subpopulations
  annotation_row_spop = data.frame(Spop = factor(list_pop_by_spop[, 1]))
  rownames(annotation_row_spop) = list_pop_by_spop[, 2]
  # Get subpopulation colors
  Spop_colores = spop_color[, 2]
  names(Spop_colores) = spop_color[, 1]
  # Create subpopulation color list
  ann_colors_spop = list(Spop = Spop_colores)
  # Convert overlap table to data matrix format
  mat_overlap <- as.data.frame.matrix(tabla_overlap)
  if (nrow(mat_overlap) > 0 &&
      ncol(mat_overlap) > 0 && any(!is.na(mat_overlap))) {
    print("The array for the heatmap has been created.")
  } else {
    stop("The array for the heatmap is empty or contains NA.")
  }
  # Generate the heatmap
  #https://jokergoo.github.io/2020/05/06/translate-from-pheatmap-to-complexheatmap/
  png(
    file = paste(path, name, sep = "/"),
    width = 3000,
    height = 3000,
    res = 300
  )
  HM <- ComplexHeatmap::pheatmap(
    t(mat_overlap),
    name = "Overlap",
    annotation_row = annotation_row_spop,
    annotation_colors = ann_colors_spop,
    display_numbers = TRUE,
    number_color = "black",
    fontsize = 7,
    fontsize_row = 8,
    fontsize_col = 8,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    border_color = "white",
    row_split = annotation_row_spop$Spop,
    legend = TRUE
  )
  draw(HM)
  dev.off()
  cat("The heatmap has been created\n")
}



# Function to calculate stability metrics for Louvain clustering resuls
stability_matrix_metrics <- function (network, R, path, n_runs = 100, use_leiden=FALSE)
{
  #' Stability Metrics for Community detection (Louvain or Leiden)
  #'
  #' This function calculates stability metrics for Louvain clustering results across multiple runs. It computes the Normalized Information Distance (NID) and Adjusted Rand Index (ARI) for pairwise comparisons of clusterings obtained from different runs, and also tracks the frequency with which pairs of nodes are assigned to the same community across runs.
  #'
  #' @param network An igraph object representing the network on which Louvain clustering is performed.
  #' @param R A numeric value specifying the resolution parameter for Louvain clustering.
  #' @param path A character string specifying the directory path where the output files will be saved.
  #' @param n_runs An integer value specifying the number of times to run the community detection algorithm. Default is 100.
  #' @param use_leiden A boolen flag to indicate whether to use the Leiden algorithm instead of Louvain. Default is FALSE (Louvain).
  #' @details
  #' The function runs the community detection algorithm `n_runs` times on the input `network`, each time with the same resolution `R`. It stores the clustering results (node memberships) for each run, and calculates pairwise comparisons between the runs using NID and ARI. The function also records how often pairs of nodes are assigned to the same cluster across the runs.
  #'
  #' The resulting pairwise NID and ARI matrices are returned, and a matrix indicating the number of times each pair of nodes is clustered together is written to a file.
  #'
  #' @return A list containing:
  #' \item{NID}{A matrix of pairwise Normalized Information Distance (NID) between runs.}
  #' \item{ARI}{A matrix of pairwise Adjusted Rand Index (ARI) between runs.}
  #'
  #' The common membership matrix is also saved to a file.
  #'
  #' @export
  #'
  #' @examples
  #' \dontrun{
  #' # Assuming 'net' is a pre-defined igraph object
  #' stability_metrics <- stability_matrix_metrics(net, R = 1.0, path = "output", n_runs = 100, use_leiden)
  #' }
  #'
  ## Select community detection algorithm
  # Initialize matrices to store NID and ARI pairwise comparisons
  mpair_NID <- matrix(nrow = n_runs, ncol = n_runs)
  mpair_ARI <- matrix(nrow = n_runs, ncol = n_runs)
  # Initialize a matrix to store community membership for each node in each run
  obj <- matrix(nrow = length(V(network)), ncol = n_runs)
  # Run the Louvain clustering algorithm for each iteration and store memberships
  for (i in 1:n_runs)
  { 
    if (use_leiden) {
      cl <- cluster_leiden(network, resolution_parameter = R)
    } else {
      cl <- cluster_louvain(network, resolution = R)
    }
    obj[, i] <- cl$membership # Store the cluster membership for each node
  }
  # Calculate NID and ARI between every pair of cluster results
  for (i in 1:n_runs)
  {
    for (j in i:n_runs)
    {
      if (i == j)
      {
        # Diagonal elements are set to NA because they are comparisons of the same clustering
        mpair_NID[i, j] <- NA
        mpair_ARI[i, j] <- NA
      }
      else
      {
        # Compute NID and ARI for each pair of runs
        mpair_NID[i, j] <- NID(obj[, i], obj[, j])
        mpair_ARI[i, j] <- ARI(obj[, i], obj[, j])
        
      }
    }
  }
  spop_membership <- as.numeric(factor(V(network)$spop))
  spop_ari <- apply(obj, 2, function(col) ARI(col, spop_membership))
  spop_nid <- apply(obj, 2, function(col) NID(col, spop_membership))
  n_com <- apply(obj, 2, function(col) max(col))
  range_nc <- range(n_com)
  mean_nc <- mean(n_com)
  std_nc <- sd(n_com)
  vector_nc <- c(mean_nc, std_nc, range_nc)
  names(vector_nc) <- c("mean", "std", "lower_range", "upper_range")
  # Initialize a matrix to store the number of times each pair of nodes shares the same cluster
  ind_mem <-
    matrix(nrow = length(V(network)),
           ncol = length(V(network)),
           data = 0)
  # Count the number of times each pair of nodes is in the same cluster over all runs
  for (i in 1:n_runs)
  {
    for (j in 1:(length(V(network)) - 1))
    {
      for (k in j:length(V(network)))
      {
        if (obj[j, i] == obj[k, i]) {
          ind_mem[j, k] = ind_mem[j, k] + 1 # Increment the count for pairs in the same cluster
          ind_mem[k, j] = ind_mem[k, j] + 1 # Symmetric update
        }
      }
    }
  }
  # Set row and column names for the membership matrix based on node names
  rownames(ind_mem) <- cl$names
  colnames(ind_mem) <- cl$names
  # Adjust the diagonal values (self-comparisons)
  diag(ind_mem) <- diag(ind_mem) / 2
  # Write the common membership matrix to a file
  write.table(
    ind_mem,
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE,
    file = paste(path, paste("/Stability_metrics/Common_membership_file_R", R, sep = ""), sep =
                   "/")
  )
  # Return the results as a list containing NID and ARI matrices
  return(list(NID = mpair_NID, ARI = mpair_ARI, spop_nid = spop_nid, spop_ari= spop_ari, ncom= vector_nc))
}


df_metrics <- function(graph, steps, path, lower_limit, upper_limit, use_leiden=FALSE) {
  #' Compute and Aggregate Stability Metrics Across Multiple Resolutions
  #'
  #' This function computes stability metrics (NID and ARI) across multiple resolution (\eqn{\lambda}) values for a given network. The stability metrics are computed using parallel processing and aggregated for each resolution. The results are split into data frames for each metric (NID and ARI).
  #'
  #' @param graph An igraph object representing the network for which stability metrics will be computed.
  #' @param steps An integer specifying the number of \eqn{\lambda} (resolution) values to evaluate.
  #' @param path A character string specifying the directory path where output files will be saved.
  #' @param lower_limit A numeric value representing the lower bound for the \eqn{\lambda} range.
  #' @param upper_limit A numeric value representing the upper bound for the \eqn{\lambda} range.
  #'
  #' @details
  #' The function calculates stability metrics for a range of resolution (\eqn{\lambda}) values. Stability metrics, such as Normalized Information Distance (NID) and Adjusted Rand Index (ARI), are computed for each resolution. The function leverages parallel processing to compute these metrics efficiently over a large number of resolutions.
  #'
  #' The stability metrics are stored as data frames and are returned for both NID and ARI. Each data frame includes the stability values for all pairs of clusterings at different \eqn{\lambda} values.
  #'
  #' @return A list containing two data frames:
  #' \item{ARI}{A data frame of pairwise Adjusted Rand Index (ARI) values across different resolutions.}
  #' \item{NID}{A data frame of pairwise Normalized Information Distance (NID) values across different resolutions.}
  #'
  #' @export
  #'
  #' @examples
  #' \dontrun{
  #' # Assuming 'graph' is a pre-defined igraph object
  #' stability_metrics_df <- df_metrics(graph, steps = 50, path = "output_directory", lower_limit = -2, upper_limit = 2, use_leiden)
  #' }
  #'  
  # Define the number of cores to use for parallel processing
  num_cores <- detectCores()
  cluster <- makeCluster(num_cores)
  registerDoParallel(cluster)
  
  # Define the range of lambda values
  lambda_values <- logspace(lower_limit, upper_limit, n = steps)
  
  # Define the function to be executed in parallel
  
  calculate_stability <- function(lambda) {
    # Compute the stability matrix for the given lambda and metric
    smn <- stability_matrix_metrics(graph, R = lambda, path = path, use_leiden = use_leiden)
    # Convert the upper triangular part of the stability matrix to a data frame
    df_NID <- as.data.frame(smn$NID[upper.tri(smn$NID, diag = FALSE)])
    df_ARI <- as.data.frame(smn$ARI[upper.tri(smn$ARI, diag = FALSE)])
    df_spop_ari <- as.data.frame(smn$spop_ari)
    df_spop_nid <- as.data.frame(smn$spop_nid)
    ncoms <-as.data.frame(t(as.data.frame(smn$ncom)))
    rownames(ncoms) <- NULL
    ncoms$lambda <- lambda
    
    colnames(df_NID) <- "value"
    colnames(df_ARI) <- "value"
    colnames(df_spop_ari) <- "value"
    colnames(df_spop_nid) <- "value"
    
    if (nrow(df_NID) > 0) 
    {
      df_NID$lambda <- lambda
      df_NID$metric <- "NID"
      df_NID$stability <- "Internal"
    }
    
    if (nrow(df_ARI) > 0) 
    {
      df_ARI$lambda <- lambda
      df_ARI$metric <- "ARI" 
      df_ARI$stability <- "Internal"
    }
    df_spop_ari$lambda <- lambda
    df_spop_nid$lambda <- lambda
    df_spop_ari$metric <- "ARI"
    df_spop_nid$metric <- "NID"
    df_spop_ari$stability <- "Spop"
    df_spop_nid$stability <- "Spop"
    single_df <- rbind(df_NID, df_ARI,df_spop_ari, df_spop_nid )
    single_df <- merge(single_df, ncoms, by = "lambda", all.y = T)
    return(single_df)
  }

  
  # Perform parallel computation using foreach
  results <- foreach(lambda = lambda_values, .combine = 'rbind', .export = c("stability_matrix_metrics"), .packages = c("igraph", "aricode")) %dopar% {
    calculate_stability(lambda)
  }
  
  stopCluster(cluster)
  
  ### get NID and ARI
  
  split_results <- split(x = results, f = results$metric)
  
  
  return(list(ARI=split_results$ARI, NID=split_results$NID))
  
}

# Function to plot and save metric analysis results
plots_metric <- function(DF_lambda, metric_name, path, name, lower_limit, upper_limit, steps){
  #' Generate Plots for Stability Metric comparing replicates and 
  #'
  #' This function generates and saves a boxplot
  #' of a specified stability metric across resolution values.
  #'
  #' @param DF_lambda A data frame containing the data to be plotted. It should include columns for \code{lambda} (resolution values) and \code{value} (metric values).
  #' @param metric_name A character string specifying the name of the metric to plot. Accepted values are \code{"nid"} or \code{"ari"} (case-insensitive).
  #' @param path A character string specifying the directory path to save the plot images.
  #' @param name A character string for naming the output files. The name will be included in the file name.
  #' @param lower_limit A numeric value specifying the lower limit of the x-axis (resolution values).
  #' @param upper_limit A numeric value specifying the upper limit of the x-axis (resolution values).
  #' @param steps An integer specifying the number of steps for the x-axis breaks between \code{lower_limit} and \code{upper_limit}.
  #'
  #' @return This function does not return any objects. It saves the plot as .png files in the specified directory and prints a confirmation message.
  #'
  #' @details The function generates four types of plots:
  #'   \itemize{
  #'     \item Boxplot of the metric values across \code{lambda} values.
  #'   }
  #' The y-axis label will correspond to the metric name specified by \code{metric_name}.
  #'
  #' @import ggplot2
  #' @import dplyr
  #'
  #' @examples
  #' \dontrun{
  #' DF_lambda <- data.frame(lambda = c(1, 2, 3), value = runif(30))
  #' plots_metric(DF_lambda, "ARI", "./plots", "stability_metric", 1, 3, 10)
  #' }
  #' 
  #' @export
  
  if (metric_name=="nid" || metric_name=="NID")
  {
    axis_name_y="Normalized Information Distance (NID)"
  }
  if (metric_name=="ari" || metric_name=="ARI"){
    
    axis_name_y="Adjusted Rand Index (ARI)"
  }
  write.csv(DF_lambda, file = paste(path,paste0(name,"_summary.csv"), sep="/"), row.names = FALSE)
  ## Number of communities
  unique_ranges <- DF_lambda %>%
    select(lambda, lower_range, upper_range, stability) %>%
    distinct() %>%
    filter(stability == "Internal")
  ### Text
  ann_text <- unique_ranges %>%
    mutate(
      log_lambda = log10(lambda),
      y_position = 1.3,  # Position above the plot
      label = paste0(lower_range,"\n",upper_range)
      
    )
  # Create the boxplot
  boxplot_2 <- ggplot(DF_lambda, aes(x = log10(lambda), y = value, color = stability, group = interaction(lambda, stability))) +
    geom_boxplot() +
    labs(x = "Resolution value", y =axis_name_y ) +
    theme_minimal() +
    theme(
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      legend.position = "right",
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 20)
    ) +
    scale_x_continuous(breaks = seq(lower_limit, upper_limit), labels = seq(lower_limit, upper_limit)) +
    scale_color_manual(
    name = "Comparison",  # New title for the legend
    values = c("#6A5ACD", "#556B2F"),  # Colors for each legend value
    labels = c("Between Communities", "Communities vs Superpopulations")  # New names for the values
  )
  ggsave(filename = paste(path,paste0(name,"_boxplot_2.png"),sep="/"), plot = boxplot_2, width = 13, height = 10, units = "in",
         bg = "white",
         dpi = 300)
  
  cat(
    "The boxplot for the stability metric",
    metric_name,
    "has been created.\n"
  )
  Internal_data <- DF_lambda %>% filter(stability == "Internal")
  Spop_data <- DF_lambda %>% filter(stability == "Spop")
  
  stats <- Internal_data %>%
    group_by(lambda) %>%
    summarise(
      n = n(),
      median_val = median(value, na.rm = TRUE),
      variance = var(value, na.rm = TRUE),
    )
  stats2 <- Spop_data %>%
    group_by(lambda) %>%
    summarise(
      n = n(),
      median_val = median(value, na.rm = TRUE),
      variance = var(value, na.rm = TRUE),
    )
  stats <- stats %>%
    mutate(Resolution = log10(lambda)) %>%
    select(-lambda)  # Remove the original `lambda` column if no longer needed
  stats2 <- stats2 %>%
    mutate(Resolution = log10(lambda)) %>%
    select(-lambda)  # Remove the original `lambda` column if no longer needed
  # Save the modified data to a CSV file
  write.csv(stats2, paste(path,"Stability_comparing_against_Spop.csv", sep="/"), row.names = FALSE)
  write.csv(stats, paste(path,"Stability_between_replicates.csv", sep="/"), row.names = FALSE)
}




# Shiny App: Interactive Map input
info_map <-
  function(path,
           im,
           cl,
           ord,
           info,
           steps,
           outputname,
           lower_limit,
           upper_limit) {
    #' Generate Population Information Table for Community Analysis
    #'
    #' This function loads metadata on individuals, processes it according to specified resolution values, and creates a formatted output table that can be used in Shiny applications for visualizing population and community data.
    #'
    #' @param path A character string specifying the directory path where the metadata and output files are stored.
    #' @param im A matrix where each element represents the community assignment for an individual at a specific resolution level.
    #' @param cl An object containing community detection results, which includes individual names.
    #' @param ord An ordering vector for individuals as used in visualization.
    #' @param info A character string specifying the filename of the metadata file (CSV format) containing individual information.
    #' @param steps An integer specifying the number of steps for the resolution values between \code{lower_limit} and \code{upper_limit}.
    #' @param outputname A character string specifying the output file name for the generated table.
    #' @param lower_limit A numeric value specifying the lower limit for resolution values.
    #' @param upper_limit A numeric value specifying the upper limit for resolution values.
    #'
    #' @return This function does not return a formated table . It saves a formatted table to the specified directory and outputs a message upon successful completion.
    #'
    #' @details The function performs several steps:
    #'   \itemize{
    #'     \item Loads metadata on individuals and filters for relevant individuals in the \code{cl$names}.
    #'     \item Orders individuals according to the order vector and creates a table for communities per individual.
    #'     \item Computes community ratios by lambda and population, adjusting column and row orders as needed.
    #'     \item Merges the formatted table with metadata to include population and region information.
    #'   }
    #' After processing, the function outputs the data in tab-separated format, suitable for visualizing population information alongside community assignments in applications like Shiny.
    #'
    #' @import dplyr
    #' @importFrom utils read.csv write.table
    #'
    #' @examples
    #' \dontrun{
    #' info_map(path = "./data",
    #'          im = community_matrix,
    #'          cl = community_data,
    #'          ord = order_vector,
    #'          info = "individual_metadata.csv",
    #'          steps = 10,
    #'          outputname = "community_info_table.txt",
    #'          lower_limit = 1,
    #'          upper_limit = 3)
    #' }
    #' 
    #' @export
    
    # Load the metadata of the individuals
    meta_file = file.path(path, info)
    # Check if files exist
    if (!file.exists(meta_file)) {
      stop("The additional information file does not exist at the specified path.")
    }
    metadata = read.csv(meta_file, sep = "\t", header = T)
    # Leave only individuals and their pop that are in cl$names
    # Remember that they enter the nodebases function in order of cl$names
    # and that in Pollock they are assigned a different place to improve the
    # visualization and that order comes in the order vector that comes out of Pollock
    met_indv_pop <- metadata[, c(1, 4, 3)]
    filtrado <-
      met_indv_pop[met_indv_pop$indv %in% as.character(cl$names),]
    # Sort like Pollock using cl$names and ord
    cl_names_ord <- cl$names[ord[, 1]]
    id_order <- match(filtrado[, 1], as.character(cl_names_ord))
    id_ordenado <- filtrado[order(id_order),]
    # Create the table for communities per individual with their info
    Lambda <-
      rep(logspace(lower_limit, upper_limit, n = steps),
          times = nrow(id_ordenado))
    # Join the populations of each individual already ordered with lambda
    pp <- rep(id_ordenado$population, each = steps)
    pp_lmbd <- cbind(pp, Lambda)
    # Since pp are ordered by indv in im, we join im with pp_lmbd
    # Thus each individual will have his community according to a lambda
    # First obtain the communities in a vector
    imcols = ncol(im)
    imhils = nrow(im)
    comms <- matrix(data = 0,
                    ncol = 1,
                    nrow = nrow(pp_lmbd))
    count = 1
    for (i in 1:imcols) {
      for (j in 1:imhils) {
        comms[count] = im[j, i]
        count = count + 1
      }
    }
    colnames(comms) <- "com"
    # Join the population of each individual that already comes with its lambda with the
    # community that you have obtained in that lambda
    pp_lmbd_comms <- cbind(pp_lmbd[, c(1, 2)], comms)
    # Get community ratios by lambda and pop
    tabla <- ftable(as.data.frame(pp_lmbd_comms))
    prop_by_row <- prop.table(as.matrix(tabla), margin = 1)
    # When the proportion table is created, the populations and lambdas are disordered
    # And the communities are not in ascending order
    # Sort communities in ascending order passing cols to rows
    tprop_by_row <- t(prop_by_row)
    numCom <- as.numeric(rownames(tprop_by_row))
    comms_order <- order(numCom)
    tprop_by_row_ord <- tprop_by_row[comms_order, ]
    # Return rows to cols and proceed to remove the names of the pops and lambdas
    prop_by_row_ord <- t(tprop_by_row_ord)
    prop_pp_lmbd <- rownames(prop_by_row_ord)
    split_pp_lmd <- strsplit(prop_pp_lmbd, "_")
    prop_pp_lmbd_df <- as.data.frame(do.call(rbind, split_pp_lmd))
    colnames(prop_pp_lmbd_df) <- c("population", "Lambda")
    # We put together the rest of the information from the metadata populations
    # First, remove duplicates from metadata
    metadata_nodup <-
      metadata[!duplicated(metadata$population, fromLast = FALSE),]
    
    pp_merge <-
      merge(prop_pp_lmbd_df,
            metadata_nodup,
            by = "population",
            all.x = TRUE)
    # Sort columns
    pp_merge2 <-
      cbind(
        pp_merge$population,
        pp_merge$pop3code,
        pp_merge$genetic_region,
        pp_merge$poject,
        pp_merge$latitude,
        pp_merge$longitude,
        pp_merge$Lambda
      )
    # Rename
    colnames(pp_merge2) <-
      c("Pop",
        "Pop3code",
        "Genetic_region",
        "Project",
        "Latitud",
        "Longitud",
        "Lambda")
    # Finally, unite the information of the populations with their communities
    prop <- prop_by_row_ord
    colnames(prop) <- paste0("C", min(im):max(im))
    # Convert zeros to NA
    prop[prop == 0] <- "NA"
    final_table <- cbind(pp_merge2, prop)
    # Sort by lambda and pop
    indice_orden_2 <-
      order(final_table[, 1], as.numeric(final_table[, 6]))
    final_table_order <- final_table[indice_orden_2,]
    if (nrow(final_table_order) > 0 &&
        ncol(final_table_order) > 0 &&
        any(!is.na(final_table_order))) {
      print("The array for the shiny app has been created.")
    } else {
      stop("The array for the shiny app is empty or contains NA.")
    }
    # Save
    final_table_file = file.path(path, outputname)
    write.table(
      final_table_order,
      final_table_file,
      sep = "\t",
      quote = F,
      row.names = F,
      col.names = T
    )
   return(final_table_order)
  }

# Heatmap: by communities
overlap_by_comm <-
  function(path, info, graph, im, ord, step, name, use_leiden=FALSE) {
    #' Generate Heatmap of Community Overlap by Population
    #'
    #' This function processes population metadata, community assignments, and graph structure to create a heatmap illustrating the overlap of communities by population. The output is a PNG heatmap visualizing the relative overlap of individuals within communities across different populations.
    #'
    #' @param path A character string specifying the directory path where the metadata and output files are stored.
    #' @param info A character string specifying the filename of the metadata file (CSV format) containing individual information.
    #' @param graph An igraph object representing the network structure of individuals.
    #' @param im A matrix of community assignments for individuals across multiple resolutions.
    #' @param ord A vector indicating the order of individuals for visualization purposes.
    #' @param step An integer representing the specific step of community assignments to analyze.
    #' @param name A character string specifying the name of the output PNG file for the heatmap.
    #' @param use_leiden Logical. If TRUE, uses the Leiden algorithm; if FALSE, uses the Louvain algorithm.
    #'
    #' @return This function does not return an object. It saves a heatmap PNG file to the specified directory and outputs a message upon successful completion.
    #'
    #' @details The function performs the following steps:
    #'   \itemize{
    #'     \item Loads metadata for individuals and filters for relevant entries in the order defined by \code{ord}.
    #'     \item Assigns communities to individuals based on their metadata and community matrix \code{im}.
    #'     \item Computes the overlap between communities and populations, normalizes values, and prepares a matrix for visualization.
    #'     \item Creates a heatmap with ComplexHeatmap, where each row represents a subpopulation, columns represent populations, and color intensity indicates the overlap.
    #'   }
    #' ComplexHeatmap is used to generate the heatmap image, which is saved to the specified directory in PNG format.
    #'
    #' @import ComplexHeatmap
    #' @importFrom utils read.csv
    #' 
    #' @examples
    #' \dontrun{
    #' overlap_by_comm(
    #'   path = "./data",
    #'   info = "individual_metadata.csv",
    #'   graph = network_graph,
    #'   im = community_matrix,
    #'   ord = order_vector,
    #'   step = 3,
    #'   name = "community_overlap_heatmap.png",
    #'   use_leiden = TRUE
    #' )
    #' }
    #'
    #' @export
    
    # Check if files exist
    if (!file.exists(file.path(path, info))) {
      stop("The additional information file does not exist at the specified path.")
    }
    # Take the individuals from cl$name
    if (use_leiden) {
      cl <- cluster_leiden(graph, resolution_parameter = 1)
    } else {
      cl <- cluster_louvain(graph, resolution = 1)
    }
    #the value is not important
    cl_name <- cl$name
    # Order the individuals as Pollock ordered them
    cl_name_ord <- cl_name[ord]
    # Metadata to have the populations
    info <-
      read.csv(file.path(path, info), sep = "\t", head = TRUE)
    info = info[, c(1, 3, 5, 8)]
    # Leave only individuals in cl_name_ord and in that order
    filtrado <- info[info[, 1] %in% cl_name_ord, ]
    id_order <- match(filtrado[, 1], cl_name_ord)
    id_ordenado <- filtrado[order(id_order),]
    # Get the communities of im
    imcols = ncol(im)
    comm <- matrix(data = 0,
                   ncol = 1,
                   nrow = ncol(im))
    count = 1
    for (j in 1:imcols) {
      comm[count, 1] = im[step, j]
      count = count + 1
    }
    # Unite id with communities
    comm <- paste0("C", comm)
    idv_comm_pop_spop <-
      cbind(id_ordenado, comm)
    idv_comm_pop_spop <- idv_comm_pop_spop[, c(1, 2, 3, 5, 4)]
    # Create the overlay table
    tabla <-
      xtabs( ~ idv_comm_pop_spop[, 4] + idv_comm_pop_spop[, 2], data = idv_comm_pop_spop)
    tabla_overlap <- tabla
    hil = 1:nrow(tabla_overlap)
    sum_by_comm = sapply(hil,  function(hil) {
      sum(tabla_overlap[hil, ])
    })
    for (hil in 1:nrow(tabla_overlap)) {
      for (col in 1:ncol(tabla_overlap)) {
        tabla_overlap[hil, col] = tabla[hil, col] / sum_by_comm[hil]
      }
    }
    # Give names depending on the number of communities
    num_2level <- length(unique(idv_comm_pop_spop[, 3]))
    list_pop_by_spop <-
      aggregate(
        idv_comm_pop_spop[, 1] ~ idv_comm_pop_spop[, 3] + idv_comm_pop_spop[, 2],
        data = idv_comm_pop_spop,
        FUN = length
      )
    # Graph with ComplexHeatmap. Info in:
    # https://jokergoo.github.io/2020/05/06/translate-from-pheatmap-to-complexheatmap/
    annotation_row_spop = data.frame(Spop = factor(list_pop_by_spop[, 1]))
    rownames(annotation_row_spop) = list_pop_by_spop[, 2]
    spop_color <-
      idv_comm_pop_spop[!duplicated(idv_comm_pop_spop[, 5], fromLast = FALSE),]
    spop_color <- spop_color[, c(3, 5)]
    Spop_colores = spop_color[, 2]
    names(Spop_colores) = spop_color[, 1]
    ann_colors_spop = list(Spop = Spop_colores)
    mat_overlap <- as.data.frame.matrix(tabla_overlap)
    if (nrow(mat_overlap) > 0 &&
        ncol(mat_overlap) > 0 && any(!is.na(mat_overlap))) {
      print("The array for the heatmap has been created.")
    } else {
      stop("The array for the heatmap is empty or contains NA.")
    }
    png(
      file = paste(path, name, sep = "/"),
      width = 3000,
      height = 3000,
      res = 300
    )
    HM <- ComplexHeatmap::pheatmap(
      t(mat_overlap),
      name = "Overlap",
      annotation_row = annotation_row_spop,
      annotation_colors = ann_colors_spop,
      display_numbers = TRUE,
      number_color = "black",
      fontsize_number = 6,
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      border_color = "white",
      row_split = annotation_row_spop$Spop,
      legend = TRUE
    )
    draw(HM)
    dev.off()
  }

# Functions to obtain the networks
# Fix that there are two edges between two communities
# 1)plot_louvain_by_comm()
# Function to graph coloring communities
# Function needed to create the table for the Shiny map
plot_louvain_by_comm <-
  function(ord, im, step, graph, mapbig, lay, name, path, use_leiden=FALSE) {
    #' Plot Louvain or Leiden Communities on Graph
    #'
    #' This function generates a plot of communities detected by the Louvain or Leiden algorithm, assigning each community a unique color and saving the plot as an PNG file.
    #'
    #' @param ord A vector indicating the order of individuals for visualization purposes.
    #' @param im A matrix of community assignments across multiple resolutions for each individual.
    #' @param step An integer specifying the row of \code{im} to use for community assignments.
    #' @param graph An igraph object representing the network of individuals.
    #' @param mapbig A data frame or matrix with community indices and their assigned colors.
    #' @param lay A layout matrix or function for plotting the igraph network.
    #' @param name A character string specifying the name of the output SVG file.
    #' @param path A character string specifying the directory path where the SVG file will be saved.
    #' @param use_leiden Logical. If TRUE, the Leiden algorithm is used; if FALSE, the Louvain algorithm is used.
    #'
    #' @return An igraph object with the community and color attributes added to each vertex. Saves the community plot as an SVG file in the specified directory.
    #'
    #' @details The function performs the following steps:
    #'   \itemize{
    #'     \item Selects the clustering algorithm (Louvain or Leiden) based on the \code{use_leiden} parameter and applies it to the graph.
    #'     \item Orders community assignments according to \code{ord}, associates colors from \code{mapbig} to communities, and applies them to the graph.
    #'     \item Verifies that the igraph object is valid and has vertices and edges, then generates a community plot using the specified layout.
    #'     \item Saves the plot as an SVG file with the filename specified by \code{name}.
    #'   }
    #'
    #' @importFrom igraph is.igraph vcount ecount V
    #' @importFrom utils write.table
    #'
    #' @examples
    #' \dontrun{
    #' plot_louvain_by_comm(
    #'   ord = order_vector,
    #'   im = community_matrix,
    #'   step = 1,
    #'   graph = network_graph,
    #'   mapbig = color_map,
    #'   lay = layout_matrix,
    #'   name = "community_plot.svg",
    #'   path = "./plots",
    #'   use_leiden = TRUE
    #' )
    #' }
    #'
    #' @export
    
    # Here step refers to the row of im from which we will get the communities
    # Sort like Pollock using cl$names and ord
    if (use_leiden) {
      cl <- cluster_leiden(graph, resolution_parameter = 1)
    } else {
      cl <- cluster_louvain(graph, resolution = 1)
    }
    # Any r is ok, we just need names
    cl_names <- cl$names
    cl_names_ord <- cl_names[ord[, 1]]
    # Put together id and its comm
    imcols = ncol(im)
    # get min from im
    min_im <- min(im)
    id_comm <- matrix(data = 0,
                      ncol = 2,
                      nrow = ncol(im))
    count = 1
    for (j in 1:imcols) {
      id_comm[count, 2] = im[step, j]
      count = count + 1
    }
    id_comm[, 1] <- cl_names_ord
    colnames(id_comm) <- c("ID", "Comm")
    comms <- sort(unique(as.numeric(id_comm[, 2])))
    # Get colors of those communities
    if(min_im==2){
    colors <- mapbig[comms-1, 1]
    }else{
     colors <- mapbig[comms, 1]
    }
    comms_colors <- cbind(comms, colors)
    # Order the communities with their color from smallest to largest
    colnames(comms_colors) <- c("Comm", "Colors")
    # Create a matrix of individuals with their community and color
    id_colors <- merge(id_comm, comms_colors, by = "Comm")
    id_colors2 <-
      cbind(id_colors[, 2], id_colors[, 1], id_colors[, 3])
    colnames(id_colors2) <- c("ID", "Comm", "Colors")
    # Organize the information with the individuals in a graph
    vert <- V(graph)$name
    id_order <- match(id_colors2[, 1], vert)
    id_ordenado <- id_colors2[order(id_order), ]
    # Assign the community and color attributes to the graph
    V(graph)$carac <- id_ordenado[, 2]
    V(graph)$color <- id_ordenado[, 3]
    # Check if it is an igraph object
    if (!is.igraph(graph)) {
      stop("The igraph object in plot_louvain_by_comm has NOT been created.")
    }
    # Check if the graph has vertices
    if (vcount(graph) == 0) {
      stop("The graph in plot_louvain_by_comm has no vertices.")
    }
    # Check if the graph has edges
    if (ecount(graph) == 0) {
      stop("The graph in plot_louvain_by_comm has no edges.")
    }
    #Plot
    png(
      file = paste(path, name, sep = "/"),
      width = 2000,
      height = 2000,
      res=300
    )
    plot(
      graph,
      layout = lay,
      vertex.size = 6,
      vertex.label = NA,
      edge.color = "gray",
      col = V(graph)$color
    )
    dev.off()
    
    cat("the ", name, " plot has been created.\n")
    return(graph)
  }

#It is almost the same as plot_louvain_by_comm(), but it does not plot.
graph_louvain_by_comm <-
  function(ord, im, step, graph, mapbig, use_leiden=FALSE) {
    #' Assign Community Colors to Graph Vertices
    #'
    #' This function assigns communities and their associated colors to vertices in an igraph object, based on the Louvain or Leiden community detection algorithm.
    #'
    #' @param ord A vector indicating the order of individuals for visualization or assignment.
    #' @param im A matrix of community assignments across multiple resolutions for each individual.
    #' @param step An integer specifying the row of \code{im} to use for community assignments.
    #' @param graph An igraph object representing the network of individuals.
    #' @param mapbig A data frame or matrix with community indices and their assigned colors.
    #' @param use_leiden Logical. If TRUE, the Leiden algorithm is used; if FALSE, the Louvain algorithm is used.
    #'
    #' @return An igraph object with added vertex attributes \code{carac} (community) and \code{color} for each individual.
    #'
    #' @details The function performs the following steps:
    #'   \itemize{
    #'     \item Selects the clustering algorithm (Louvain or Leiden) based on the \code{use_leiden} parameter and applies it to the graph.
    #'     \item Orders community assignments according to \code{ord} and applies community colors from \code{mapbig} to the vertices of the graph.
    #'     \item Verifies that the igraph object has been successfully created and contains vertices and edges.
    #'   }
    #'
    #' @importFrom igraph is.igraph vcount ecount V
    #'
    #' @examples
    #' \dontrun{
    #' graph_louvain_by_comm(
    #'   ord = order_vector,
    #'   im = community_matrix,
    #'   step = 1,
    #'   graph = network_graph,
    #'   mapbig = color_map,
    #'   use_leiden = TRUE
    #' )
    #' }
    #'
    #' @export
    
    # Here step refers to the row of im from which we will get the communities
    # Sort like Pollock using cl$names and ord
    if (use_leiden) {
      cl <- cluster_leiden(graph, resolution_parameter = 1)
    } else {
      cl <- cluster_louvain(graph, resolution = 1)
    }
    # Any r is ok, we just need names
    cl_names <- cl$names
    cl_names_ord <- cl_names[ord[, 1]]
    # Put together id and its comm
    imcols = ncol(im)
    min_im <- min(im)
    id_comm <- matrix(data = 0,
                      ncol = 2,
                      nrow = ncol(im))
    count = 1
    for (j in 1:imcols) {
      id_comm[count, 2] = im[step, j]
      count = count + 1
    }
    id_comm[, 1] <- cl_names_ord
    colnames(id_comm) <- c("ID", "Comm")
    comms <- sort(unique(as.numeric(id_comm[, 2])))
    # Get the colors of those communities
    if (min_im==2){
    colors <- mapbig[comms-1, 1]
    }else{
    colors <- mapbig[comms, 1]
    }
    comms_colors <- cbind(comms, colors)
    # Order the communities with their color from smallest to largest
    colnames(comms_colors) <- c("Comm", "Colors")
    # Create a matrix of individuals with their community and color
    id_colors <- merge(id_comm, comms_colors, by = "Comm")
    id_colors2 <-
      cbind(id_colors[, 2], id_colors[, 1], id_colors[, 3])
    colnames(id_colors2) <- c("ID", "Comm", "Colors")
    # Organize the information with the individuals in a graph
    vert <- V(graph)$name
    id_order <- match(id_colors2[, 1], vert)
    id_ordenado <- id_colors2[order(id_order), ]
    # Assign the community and color attributes to the graph
    V(graph)$carac <- id_ordenado[, 2]
    V(graph)$color <- id_ordenado[, 3]
    # Check if it is an igraph object
    if (!is.igraph(graph)) {
      stop("The igraph object in plot_louvain_by_comm has NOT been created.")
    }
    # Check if the graph has vertices
    if (vcount(graph) == 0) {
      stop("The graph in plot_louvain_by_comm has no vertices.")
    }
    # Check if the graph has edges
    if (ecount(graph) == 0) {
      stop("The graph in plot_louvain_by_comm has no edges.")
    }
    return(graph)
  }

# Here we need the graph that is generated in
# plot_louvain_by_comm()
# Get average position of x and y from lay
# vertex_attr_names(graph) helps see what attributes are there
media_graph <- function(graph_by_comm, lay) {
  #' Calculate Average Position of Community Centers
  #'
  #' This function computes the average x and y positions for each community in a graph layout, based on the provided community assignments.
  #'
  #' @param graph_by_comm An igraph object with a \code{carac} vertex attribute representing community assignments.
  #' @param lay A matrix or data frame representing the layout of the graph, where each row corresponds to the x and y coordinates of each vertex.
  #'
  #' @return A data frame with the average x and y coordinates for each community, providing a summary of each community's central position in the graph layout.
  #'
  #' @details This function:
  #'   \itemize{
  #'     \item Extracts the \code{carac} community assignment for each vertex.
  #'     \item Aggregates the x and y positions across vertices within each community.
  #'     \item Returns a data frame with the average x and y coordinates for each community.
  #'   }
  #' 
  #' @examples
  #' \dontrun{
  #' media_graph(graph_by_comm = community_graph, lay = layout_matrix)
  #' }
  #'
  #' @importFrom igraph V
  #'
  #' @export
  
  pos <- cbind(V(graph_by_comm)$carac, lay)
  pos_x <- aggregate(as.numeric(pos[, 2]) ~ pos[, 1],
                     data = pos,
                     FUN = mean)
  pos_y <- aggregate(as.numeric(pos[, 3]) ~ pos[, 1],
                     data = pos,
                     FUN = mean)
  means_table <- cbind(pos_x, pos_y[, 2])
  if (nrow(means_table) == 0 ||
      ncol(means_table) == 0 || all(is.na(means_table))) {
    stop("The means_table array is empty or contains only NA.")
  }
  return(means_table)
}

# Number of nodes per community
num_vert <- function(graph_by_comm) {
  #' Calculate Number of Nodes per Community
  #'
  #' This function computes the number of nodes (vertices) in each community within a graph object.
  #'
  #' @param graph_by_comm An igraph object where the \code{carac} vertex attribute indicates the community assignment for each node.
  #'
  #' @return A matrix with two columns:
  #'   \itemize{
  #'     \item \code{Comm}: The community identifier.
  #'     \item \code{Frecuencia}: The number of nodes in each community.
  #'   }
  #' The rows are ordered by community identifiers.
  #'
  #' @details This function:
  #'   \itemize{
  #'     \item Extracts community information from the \code{carac} attribute in the graph vertices.
  #'     \item Computes the frequency of nodes in each community and returns an ordered matrix.
  #'   }
  #'
  #' @examples
  #' \dontrun{
  #' num_vert(graph_by_comm = community_graph)
  #' }
  #'
  #' @importFrom igraph V
  #'
  #' @export
  
  # Get the frequency table
  tabla_frecuencias <- table(V(graph_by_comm)$carac)
  # Convert the frequency table into a data frame
  df_frecuencias <- as.data.frame(tabla_frecuencias)
  df_frecuencias <- as.matrix(df_frecuencias)
  # Rename columns
  colnames(df_frecuencias) <- c("Comm", "Frecuencia")
  ord <- order(as.numeric(df_frecuencias[, 1]))
  df_frecuencias_ord <- df_frecuencias[ord,]
  # Make freq_ordenado a matrix
  freq_ordenado <- matrix(0, nrow = length(ord), ncol = 2)
  for (i in 1:length(ord)) {
    freq_ordenado[i, 1] = df_frecuencias[i, 1]
    freq_ordenado[i, 2] = df_frecuencias[i, 2]
  }
  if (nrow(freq_ordenado) == 0 ||
      ncol(freq_ordenado) == 0 || all(is.na(freq_ordenado))) {
    stop("The freq_ordenado array is empty or contains only NA.")
  }
  return(freq_ordenado)
}

# Connection density
# Attributes: edge_attr_names(graph)
dens_con <- function(graph_by_comm) {
  #' Calculate Connection Density Between Communities
  #'
  #' This function computes the density of connections between different communities in a graph.
  #'
  #' @param graph_by_comm An igraph object where each vertex has a community assignment stored in the \code{carac} attribute.
  #'
  #' @return A matrix with three columns:
  #'   \itemize{
  #'     \item \code{from}: The originating community of an inter-community edge.
  #'     \item \code{to}: The target community of an inter-community edge.
  #'     \item \code{weight}: The count of edges between the \code{from} and \code{to} communities.
  #'   }
  #' The matrix excludes intra-community connections.
  #'
  #' @details This function:
  #'   \itemize{
  #'     \item Merges community information with edge lists to identify inter-community connections.
  #'     \item Removes self-loop connections within the same community.
  #'     \item Returns a non-zero weighted matrix showing connection density between communities.
  #'   }
  #'
  #' @examples
  #' \dontrun{
  #' dens_con(graph_by_comm = community_graph)
  #' }
  #'
  #' @importFrom igraph V get.edgelist
  #'
  #' @export
  
  # Obtain pairs of individuals from the graph and name
  conections <- get.edgelist(graph_by_comm)
  colnames(conections) <- c("ID1", "ID2")
  # Obtain a list of individuals and the community to which they belong
  indv_comm <- cbind(V(graph_by_comm)$name, V(graph_by_comm)$carac)
  colnames(indv_comm) <- c("ID1", "Comm")
  # Get the communities from ID1
  merge_1 <- merge(conections, indv_comm, by = "ID1")
  colnames(merge_1) <- c("ID1", "ID2", "Comm1")
  # We change name to obtain communities for ID2
  colnames(indv_comm) <- c("ID2", "Comm")
  merge_2 <- merge(merge_1, indv_comm, by = "ID2")
  # Sort and name the matrix
  merged <- merge_2[, c(2, 1, 3, 4)]
  colnames(merged) <- c("ID1", "ID2", "Comm1", "Comm2")
  # Now we are only left with Com1 and Com2
  comms <- cbind(merged$Comm1, merged$Comm2)
  # Alphabetically sort the entity names in each row
  relaciones_ordenadas <- t(apply(comms, 1, function(row) {
    if (row[1] < row[2]) {
      return(c(row[1], row[2]))
    } else {
      return(c(row[2], row[1]))
    }
  }))
  # Convert the resulting matrix into a data frame
  relaciones_ordenadas <- as.data.frame(relaciones_ordenadas)
  # Obtain frequency table
  vector1 <- relaciones_ordenadas[, 1]
  vector2 <- relaciones_ordenadas[, 2]
  combined_levels <- unique(c(vector1, vector2))
  freq_table <-
    table(
      factor(vector1, levels = combined_levels),
      factor(vector2, levels = combined_levels)
    )
  # Shows the frequency table
  # There are connections within communities
  # We don't need those and we make them 0
  for (i in 1:ncol(freq_table)) {
    freq_table[i, i] = 0
  }
  tabla_freq2 <- as.data.frame(freq_table)
  # In the frequency table all against all are integrated
  # So there are many zeros left. We remove them.
  if (nrow(tabla_freq2) > 1) {
    nonzeros2 <- which(tabla_freq2[, 3] != 0)
    tabla_freq_nonzeros2 <- tabla_freq2[nonzeros2, ]
    tabla_freq_nonzeros2 <- as.matrix(tabla_freq_nonzeros2)
    colnames(tabla_freq_nonzeros2) <- c("from", "to", "weight")
  } else{
    tabla_freq_nonzeros2 <- tabla_freq2
    tabla_freq_nonzeros2 <- as.matrix(tabla_freq_nonzeros2)
    colnames(tabla_freq_nonzeros2) <- c("from", "to", "weight")
  }
  if (nrow(tabla_freq_nonzeros2) == 0 ||
      ncol(tabla_freq_nonzeros2) == 0 ||
      all(is.na(tabla_freq_nonzeros2))) {
    stop("The tabla_freq_nonzeros2 array is empty or contains only NA.")
  }
  return(tabla_freq_nonzeros2)
}

# Make the graph
graph_comms <-
  function(graph_by_comm,
           tabla_freq_nonzeros,
           freq_ordenado) {
    #' Create a Community-Level Graph
    #'
    #' This function constructs a new graph where each vertex represents a community from the original graph,
    #' with edges weighted by the frequency of connections between communities.
    #'
    #' @param graph_by_comm An igraph object with community assignments stored as vertex attributes.
    #' @param tabla_freq_nonzeros A matrix with three columns:
    #'   \itemize{
    #'     \item \code{from}: Originating community of an inter-community edge.
    #'     \item \code{to}: Target community of an inter-community edge.
    #'     \item \code{weight}: Frequency of edges between the communities.
    #'   }
    #' @param freq_ordenado A matrix with two columns:
    #'   \itemize{
    #'     \item \code{Comm}: Community identifier.
    #'     \item \code{Frecuencia}: Number of nodes in each community.
    #'   }
    #'
    #' @return An igraph object representing the community-level graph, with the following attributes:
    #'   \itemize{
    #'     \item \code{carac}: Community identifier for each vertex.
    #'     \item \code{size}: Size of each community based on \code{freq_ordenado}.
    #'     \item \code{color}: Color assigned to each community.
    #'     \item \code{weight}: Frequency of connections between communities, used to define edge width.
    #'   }
    #'
    #' @details This function:
    #'   \itemize{
    #'     \item Creates a community-level graph based on edge frequency between communities.
    #'     \item Assigns attributes to vertices and edges based on community properties and connection weights.
    #'   }
    #'
    #' @examples
    #' \dontrun{
    #' graph_comms(graph_by_comm = community_graph, tabla_freq_nonzeros = edge_matrix, freq_ordenado = community_sizes)
    #' }
    #'
    #' @importFrom igraph V E graph_from_data_frame
    #'
    #' @export
    
    
    # Obtain the names of the communities from smallest to largest
    vert <- sort(as.numeric(unique(V(graph_by_comm)$carac)))
    # Get the colors of the communities
    color <- unique(V(graph_by_comm)$color)
    # Obtain the order of the colors to match the position of the communities
    orden <- order(as.numeric(unique(V(graph_by_comm)$carac)))
    color <- color[orden]
    # Create a graph with this information
    graph_only_comms <-
      graph_from_data_frame(tabla_freq_nonzeros,
                            directed = F,
                            vertices = vert)
    # Assign attributes to graph
      if(nrow(freq_ordenado) > 1) {               
      order_freq <-
        match(freq_ordenado[, 1], vert)
        freq_ordenado2 <-
        as.matrix(freq_ordenado[order(order_freq), ])   
    } else{freq_ordenado2 <- freq_ordenado}
    V(graph_only_comms)$carac <- vert
    V(graph_only_comms)$size <- as.numeric(freq_ordenado2[, 2])
    V(graph_only_comms)$color <- color
    E(graph_only_comms)$weight <-
      as.numeric(E(graph_only_comms)$weight)
    E(graph_only_comms)$width <- E(graph_only_comms)$weight
    # Check if it is an igraph object
    if (!is.igraph(graph_only_comms)) {
      stop("The igraph object in graph_comms has NOT been created.")
    }
    # Check if the graph has vertices
    if (vcount(graph_only_comms) == 0) {
      stop("The graph in graph_comms has no vertices.")
    }
    # Check if the graph has edges
    if (ecount(graph_only_comms) == 0) {
      stop("The graph in graph_comms has no edges.")
    }
    return(graph_only_comms)
  }

# Function to count unique numbers per row
contar_unicos <- function(fila) {
  longitud <- length(unique(fila))
  return(longitud)
}
# 3D layout tights
media_graph_3D <- function(graph_by_comm, lay_3D) {
  #' Calculate Average 3D Positions for Each Community
  #'
  #' This function computes the mean 3D coordinates (x, y, z) for each community in a 3D layout,
  #' providing a summary table with the average positions for visualization or further analysis.
  #'
  #' @param graph_by_comm An igraph object with community assignments stored as vertex attributes.
  #' @param lay_3D A 3D layout matrix or dataframe with three columns representing the x, y, and z coordinates for each vertex in \code{graph_by_comm}.
  #'
  #' @return A matrix containing the mean 3D positions of each community, with the following columns:
  #'   \itemize{
  #'     \item \code{pos[,1]}: Community identifier.
  #'     \item \code{pos_x[,2]}: Mean x-coordinate.
  #'     \item \code{pos_y[,2]}: Mean y-coordinate.
  #'     \item \code{pos_z[,2]}: Mean z-coordinate.
  #'   }
  #'
  #' @details This function:
  #'   \itemize{
  #'     \item Extracts community identifiers and 3D coordinates for each vertex.
  #'     \item Aggregates the coordinates by community and computes mean x, y, and z positions.
  #'     \item Returns a summary matrix of these means.
  #'
  #' @examples
  #' \dontrun{
  #' lay_3D <- layout_with_fr(graph, dim = 3) # example 3D layout
  #' media_graph_3D(graph_by_comm = community_graph, lay_3D = lay_3D)
  #' }
  #'
  #' @importFrom igraph V
  #' @importFrom stats aggregate
  #'
  #' @export
  
  pos <- cbind(V(graph_by_comm)$carac, lay_3D)
  pos_x <- aggregate(as.numeric(pos[, 2]) ~ pos[, 1],
                     data = pos,
                     FUN = mean)
  pos_y <- aggregate(as.numeric(pos[, 3]) ~ pos[, 1],
                     data = pos,
                     FUN = mean)
  pos_z <- aggregate(as.numeric(pos[, 4]) ~ pos[, 1],
                     data = pos,
                     FUN = mean)
  means_table_3D <- cbind(pos_x, pos_y[, 2], pos_z[, 2])
  if (nrow(means_table_3D) == 0 ||
      ncol(means_table_3D) == 0 || all(is.na(means_table_3D))) {
    stop("The means_table_3D array is empty or contains only NA.")
  }
  return(means_table_3D)
}

# The plot_by_density function takes as input a graph with communities (without community 1),
# a table of means, a file path, and a name for the output file.
# The function sorts the nodes in the graph based on the mean table,
# adjusts the size of the nodes and the width of the edges based on their relative ranking,
# and then graphs and saves the resulting display.
plot_by_density <-
  function(graph_only_comms_no1,
           means_table,
           path,
           plotname) {
    #' Plot Graph by Community Density
    #'
    #' The `plot_by_density` function takes as input a graph with communities (excluding community 1),
    #' a table of average positions (means), a file path, and a file name for the output plot.
    #' The function arranges the nodes in the graph based on the provided mean positions, 
    #' adjusts the size of the nodes and the width of the edges based on their relative ranking, 
    #' and then plots and saves the graph as an image file.
    #'
    #' @param graph_only_comms_no1 An igraph object representing the community graph excluding community 1.
    #' @param means_table A data frame containing the mean positions for each community.
    #' @param path A character string specifying the file path where the plot image will be saved.
    #' @param plotname A character string specifying the name of the output image file.
    #'
    #' @details
    #' This function sorts the nodes in `graph_only_comms_no1` according to the community ordering in `means_table`.
    #' It then adjusts the node sizes and edge widths based on their rank within the graph. The resulting 
    #' plot is saved as a PNG file at the specified location.
    #'
    #' @return
    #' This function does not return a value; it saves the plot as a PNG image.
    #'
    #' @examples
    #' # Assuming `graph_only_comms_no1`, `means_table`, `path`, and `plotname` are defined:
    #' plot_by_density(graph_only_comms_no1, means_table, "path/to/save", "output_plot.png")
    #'
    #' @export
    
    # Sort means_table
    comm_order <-
      match(as.numeric(V(graph_only_comms_no1)$carac), as.numeric(means_table[, 1]))
    means_table_ord_filt <- means_table[comm_order,]
    lay_mean <- as.matrix(means_table_ord_filt[, c(2, 3)])
    # Adjusts the size of the nodes and the width of the edges based on their relative ranking
    vsize = (rank(V(graph_only_comms_no1)$size) / length(V(graph_only_comms_no1)$size)) *
      50
    ewidth = (rank(E(graph_only_comms_no1)$width) / length(E(graph_only_comms_no1)$width)) *
      20
    if (length(rank(V(graph_only_comms)$size)) == 1) {
      ecolor = NA
    } else{
      ecolor = "gray"
    }
    # Graph and save
    png(file.path(path, plotname), res = 300, width = 2000, height = 2000)
    plot(
      graph_only_comms_no1,
      layout = lay_mean,
      vertex.size = vsize,
      edge.width = ewidth,
      edge.color = ecolor,
      col = V(graph_only_comms_no1)$color,
    )
    dev.off()
    cat("The ", plotname, " plot has been created.\n")
  }

make_res_plot <- function(lay_3D_colors, im, path, plot_name, lower_limit, upper_limit){
  #' Generate Resolution Plot for Community Analysis
  #'
  #' The `make_res_plot` function generates a resolution plot based on a 3D layout color mapping, 
  #' a matrix of community data, and specified resolution limits. This plot is saved as a PNG image.
  #'
  #' @param lay_3D_colors A matrix with colors corresponding to the communities, where each row 
  #' represents a unique color for each community in 3D layout.
  #' @param im A matrix representing community information, with rows corresponding to 
  #' resolution levels and columns to individuals.
  #' @param path A character string specifying the directory path where the plot image will be saved.
  #' @param plot_name A character string specifying the name of the output PNG file.
  #' @param lower_limit An integer specifying the minimum resolution value for the y-axis.
  #' @param upper_limit An integer specifying the maximum resolution value for the y-axis.
  #'
  #' @details
  #' This function reverses the rows of `im` for proper orientation and uses the values to generate a plot
  #' with individuals on the x-axis and resolution values on the y-axis. Each individual's color is determined
  #' by their respective community.
  #'
  #' @return
  #' This function does not return a value; it saves the plot as a PNG image.
  #'
  #' @examples
  #' # Assuming `lay_3D_colors`, `im`, `path`, `plot_name`, `lower_limit`, and `upper_limit` are defined:
  #' make_res_plot(lay_3D_colors, im, "path/to/save", "res_plot.png", 1, 10)
  #'
  #' @export
  
  im2 <- apply(im, 2, rev)
  # fix axis values
  x_vals <- 1:ncol(im)
  y_vals <- 1:nrow(im)
  png(file = paste(path, plot_name, sep = "/"), res=300, width = 2000, height = 2000)
  image(
    x = x_vals,
    y = y_vals,
    z = t(im2),
    useRaster = TRUE,
    col = lay_3D_colors[, 2],
    xaxs = "i",
    yaxt = "n",
    xlab = "Individuals",
    ylab = " ",
    las = 1,
    cex.lab = 1.6,
    cex.axis = 1.6
  )
  mtext(
    side = 2,
    line = 2,
    text = "Resolution value",
    cex = 1.6,
  )
  yticklabels <- rev(lower_limit:upper_limit)
  # Crear el vector de posiciones para los ticks del eje y
  yticks <- seq(1, nrow(im), length.out = length(yticklabels))
  # Configurar los ticks y etiquetas del eje y
  axis(
    side = 2,
    at = yticks,
    labels = yticklabels,
    las = 1,
    cex.axis = 1.6
  )
  dev.off()
}

plot_indv_no_comms1 <- function(path, name_no_comm1, graph_by_comm_no1, lay_no_comm1) {
  #' Plot Graph of Individuals Excluding Community 1
  #'
  #' The `plot_indv_no_comms1` function generates a 2D plot of a graph that excludes individuals from 
  #' community 1 and saves it as an SVG file.
  #'
  #' @param path A character string specifying the directory path where the plot will be saved.
  #' @param name_no_comm1 A character string specifying the file name for the SVG plot.
  #' @param graph_by_comm_no1 An igraph object representing the graph of individuals, excluding those 
  #' in community 1.
  #' @param lay_no_comm1 A matrix representing the layout coordinates for each vertex in `graph_by_comm_no1`.
  #'
  #' @details
  #' This function creates a 2D SVG plot of the graph, positioning vertices according to the provided layout. 
  #' Vertex size is set to 6, edges are colored gray, and each vertex is colored according to its community 
  #' attribute in the graph.
  #'
  #' @return
  #' This function does not return a value; it saves the plot as an SVG file.
  #'
  #' @examples
  #' # Assuming `path`, `name_no_comm1`, `graph_by_comm_no1`, and `lay_no_comm1` are defined:
  #' plot_indv_no_comms1("path/to/save", "graph_no_comm1.svg", graph_by_comm_no1, lay_no_comm1)
  #'
  #' @export
  
  png(
    file = paste(path, name_no_comm1, sep = "/"),
    width = 2000,
    height = 2000,
    res = 300
  )
  plot(
    graph_by_comm_no1,
    layout = lay_no_comm1,
    vertex.size = 6,
    vertex.label = NA,
    edge.color = "gray",
    col = V(graph_by_comm_no1)$color
  )
  dev.off()
  cat("the ", name_no_comm1, " plot has been created.\n")
}

make_legend <-
  function(colors,
           legendcols,
           val_width,
           val_height,
           val_pt.cex,
           val_cex) {
    
    #' Create Legend for Community Colors
    #'
    #' The `make_legend` function generates a legend with specified colors and labels for each community
    #' and saves it as a PNG image.
    #'
    #' @param colors A matrix or data frame containing colors to use in the legend. Each row represents 
    #' a color associated with a community.
    #' @param legendcols An integer specifying the number of columns in the legend.
    #' @param val_width An integer indicating the width of the PNG image for the legend.
    #' @param val_height An integer indicating the height of the PNG image for the legend.
    #' @param val_pt.cex A numeric value controlling the size of the points in the legend.
    #' @param val_cex A numeric value controlling the text size in the legend.
    #'
    #' @details
    #' This function saves a PNG file displaying a legend where each entry represents a community with its
    #' corresponding color. The colors are filled in circular symbols with black borders, and the labels are 
    #' numbered according to the communities.
    #'
    #' @return
    #' This function does not return a value; it saves the legend as a PNG file.
    #'
    #' @examples
    #' # Assuming `colors`, `legendcols`, `val_width`, `val_height`, `val_pt.cex`, and `val_cex` are defined:
    #' make_legend(colors, 3, 300, 200, 1.5, 1)
    #'
    #' @export
    
    
    png(
      file = paste(path, name_legend, sep = "/"),
      width = val_width,
      height = val_height,
      res = 300
    )
    if (colors[1,]=="#FFFFFF"){
    labels <- c(1:nrow(colors))
    }else{
    labels <- c(2:(nrow(colors)+1))
    }

    plot.new()
    legend(
      "center",
      legend = labels,
      pt.bg = colors[, 1],
      bty = "n",
      pch = 21,
      pt.cex = val_pt.cex,
      cex = val_cex,
      horiz = FALSE,
      inset = c(0.1, 0.1),
      col = "black",
      ncol = legendcols,
    )
    dev.off()
  }

###############################################################################
################################SCRIPT#########################################
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
use_leiden <- "-L" %in% args
args <- args[args != "-L"]

kind = args[1]
path = args[2]
data = args[3]
info = args[4]
max = as.numeric(args[5])
steps = as.numeric(args[6])
lambda = args[7]
prune = as.numeric(args[8])
min_comms = as.numeric(args[9])
lower_limit = as.numeric(args[10])
upper_limit = as.numeric(args[11])

# Optional random_seed
if (length(args) >= 12 && !is.na(as.numeric(args[12]))) {
  random_seed = as.numeric(args[12])
  if (!is.na(random_seed)) {
    print(paste("Random seed set to:", random_seed))
  } else {
    print("Invalid random seed provided.")
  }
} else {
  print("No random seed provided. Using default behavior.")
  random_seed <- NULL
}


r = 0 # r is which layer (counting from bottom) will be perfectly sorted


# Creating the graph
network <- input(kind, path, data, info, max, prune)
graph <- network[[1]]
spop_color <- network[[2]]
if (!is.null(random_seed)) {
  set.seed(random_seed)
}
lay <- layout_with_fr(graph)


# Plot individual network with predifined groups
NetworksIndividuals_dir = file.path(path, "NetworksIndividuals")
if (!file.exists(NetworksIndividuals_dir)) {
  # If it does not exist, create the directory
  dir.create(NetworksIndividuals_dir)
  cat("'NetworksIndividuals' folder created successfully.\n")
}
Requested_RIndex_dir = file.path(path, "NetworksIndividuals", "Requested_RIndex")
if (!file.exists(Requested_RIndex_dir)) {
  # If it does not exist, create the directory
  dir.create(Requested_RIndex_dir)
  cat("'Requested_RIndex' folder created successfully.\n")
}
png(
  paste(
    NetworksIndividuals_dir,
    "Network_individuals_PredefinedGroups.png",
    sep = "/"
  )
)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(
  graph,
  vertex.size = 6,
  layout = lay,
  vertex.label = NA
)
legend(
  1.3, 1,
  legend = unique(spop_color[, 1]) ,
  bty = "n",
  pch = 20,
  pt.cex = 1,
  cex = 1,
  horiz = FALSE,
  inset = c(0.1, 0.1),
  col = spop_color[, 2]
)
dev.off()

Community_detection_dir = file.path(path, "Community_detection")
if (!file.exists(Community_detection_dir)) {
  # If it does not exist, create the directory
  dir.create(Community_detection_dir)
  cat("'Community_detection' folder created successfully.\n")
}

# Running the Louvain algorithm at different resolution values
cl_list <-
  get_lambda_results(graph,
                     steps,
                     path,
                     "/Community_detection/Result_matrix_raw.txt",
                     lower_limit = lower_limit,
                     upper_limit = upper_limit,
                     use_leiden = use_leiden,
                     random_seed = random_seed)

# Plot: Resolution plot
Mfile <- file.path(path, "/Community_detection/Result_matrix_raw.txt")
M <- as.matrix(read.csv(Mfile, sep = ",", header = FALSE))
res = nodebased(M)

resfile = file.path(path, "/Community_detection/Result_matrix_relabeled_allcomm.txt")
write.table(
  as.data.frame(res),
  resfile,
  sep = ",",
  quote = F,
  row.names = FALSE,
  col.names = FALSE
)

### Browser file
Browser_files_dir = file.path(path, "Browser_files")
if (!file.exists(Browser_files_dir)) {
  # If it does not exist, create the directory
  dir.create(Browser_files_dir)
  cat("'Browser_files' folder created successfully.\n")
}
### Pollock
resfile = file.path(path, "/Community_detection/Result_matrix_relabeled_allcomm.txt")
Ssort <- as.matrix(read.csv(resfile, sep = ",", header = FALSE))

res_list <- pollock(Ssort, min_comms, r, path, upper_limit=upper_limit, lower_limit=lower_limit)
res_mat <- res_list[[1]]
res_indv <- res_list[[2]]
res_colors <- res_list[[3]]



#### Give ids instead of index
if (use_leiden) {
  cl <- cluster_leiden(graph, resolution_parameter = 1)
} else {
  cl <- cluster_louvain(graph, resolution = 1)
}
cl_names <- cl$names
cl_names_ord <- cl_names[res_indv]


# Save the principal objects
im_file = file.path(path, "Browser_files/Resolution_matrix.txt")
or_file = file.path(path, "Community_detection/individual_index_order.txt")
or_file_ids = file.path(path, "Community_detection/individual_ID_order.txt")
ma_file = file.path(path, "Browser_files/Distinctive_colors.txt")
write.table(
  res_mat,
  im_file,
  sep = ",",
  quote = F,
  row.names = F,
  col.names = F
)
write.table(
  res_indv,
  or_file,
  sep = ",",
  quote = F,
  row.names = F,
  col.names = F
)
write.table(
  cl_names_ord,
  or_file_ids,
  sep = ",",
  quote = F,
  row.names = F,
  col.names = F 
)
write.table(
  res_colors,
  ma_file,
  sep = ",",
  quote = F,
  row.names = F,
  col.names = F
)

### Stability plots
Stability_metrics_dir = file.path(path, "Stability_metrics")
if (!file.exists(Stability_metrics_dir)) {
  # If it does not exist, create the directory
  dir.create(Stability_metrics_dir)
  cat("'Stability_metrics' folder created successfully.\n")
}

results_stability <- df_metrics(
  graph,
  steps = steps,
  path = path,
  lower_limit = lower_limit,
  upper_limit = upper_limit,
  use_leiden = use_leiden
)

plots_metric(
  DF_lambda = results_stability$NID,
  metric_name = "NID",
  path = paste0(path, "/Stability_metrics/"),
  name = "Stability_NID",
  lower_limit = lower_limit,
  upper_limit = upper_limit,
  steps = steps
)

plots_metric(  
  DF_lambda = results_stability$ARI,
  metric_name = "ARI",
  path = paste0(path, "/Stability_metrics/"),
  name = "Stability_ARI",
  lower_limit = lower_limit,
  upper_limit = upper_limit,
  steps = steps
)

# Make for individuals network, and Heatmaps 
# Read Lambda values
lambda_path = file.path(path, lambda)
lambda_values <-
  as.vector(unlist(read.table(
    lambda_path, sep = ",", header = F
  )))
for (lambda_index in lambda_values) {
  # Get lambda value
  Lambda <- logspace(lower_limit, upper_limit, steps)[lambda_index]
  NetworksIndividuals_dir = file.path(path, "NetworksIndividuals")
  if (!file.exists(NetworksIndividuals_dir)) {
    # If it does not exist, create the directory
    dir.create(NetworksIndividuals_dir)
    cat("'NetworksIndividuals' folder created successfully.\n")
  }
  Requested_RIndex_dir = file.path(NetworksIndividuals_dir, "Requested_RIndex")
  if (!file.exists(Requested_RIndex_dir)) {
    # If it does not exist, create the directory
    dir.create(Requested_RIndex_dir)
    cat("'Requested_RIndex' folder created successfully.\n")
  }
  # Plot: Network
  plot_louvain_by_spop(
    graph = graph,
    cl_list = cl_list,
    R = Lambda,
    lay = lay,
    spop_color = spop_color,
    name = paste0(
      "NetworksIndividuals/Requested_RIndex/NetworkIndividuals_PredefinedGroups_CommShaded_RIndex",
      lambda_index,
      ".png"
    )
  )
  # Heatmap
  overlap_hm(
    path,
    info,
    cl_list = cl_list,
    R = Lambda,
    spop_color,
    name = paste0("heatmap_", lambda_index, ".png")
  )
  
  # Heatmap with comm 1
  overlap_by_comm(
    path,
    info,
    graph,
    res_mat,
    res_indv ,
    step = lambda_index,
    paste0("heatmap_white_comm_", lambda_index, ".png"),
    use_leiden = use_leiden
  )
  
  # Heatmap membership
  hmp <-
    read.table(file = paste(
      path,
      paste("/Stability_metrics/Common_membership_file_R", Lambda, sep = ""),
      sep = "/"
    ), header = TRUE, check.names= FALSE)
  
  heatmap_info <-
    read.csv(file.path(path, info), sep = "\t", head = TRUE)
  
  heatmap_info <-
    heatmap_info[heatmap_info[, 1] %in% colnames(hmp), ]
  
  d <- data.frame(heatmap_info[, 5], row.names = heatmap_info[, 1])
  colnames(d) <- "Groups"
  Spop_colores <- spop_color[, 2]
  names(Spop_colores) <- spop_color[, 1]
  colors_spop = list(Groups = Spop_colores)
  
  png(
    file = paste(
      path,
      paste("/Stability_metrics/Common_membership_R", lambda_index, ".png", sep = ""),
      sep = "/"
    ),
    res = 300,
    width = 3000,
    height = 3000
  )
  ph_com <-
    pheatmap(
      as.matrix(hmp),
      use_raster = TRUE,
      annotation_row = d,
      show_rownames = FALSE,
      show_colnames = FALSE,
      name = "Common membership",
      show_row_den = FALSE,
      annotation_colors = colors_spop,
      fontsize = 10,
      annotation_col = d
    )
  draw(ph_com)
  dev.off()
}

### Input: Interactive Map:
# Read files if necessary:
imfile = file.path(path, "Browser_files/Resolution_matrix.txt")
im <- read.csv2(imfile, sep = ",", header = F)
order_file = file.path(path, "/Community_detection/individual_index_order.txt")
ord <- read.csv2(order_file, sep = ",", header = F)
cl <- cl_list[[as.character(Lambda)]]
ma_file = file.path(path, "Browser_files/Distinctive_colors.txt")
mapbig <- read.csv2(ma_file, sep = ",", header = F)

# Calling the function: shiny info
Browser_files_dir = file.path(path, "Browser_files")
if (!file.exists(Browser_files_dir)) {
  # If it does not exist, create the directory
  dir.create(Browser_files_dir)
  cat("'Browser_files' folder created successfully.\n")
}
info_obj <- info_map(
  path = path,
  cl = cl,
  im = im,
  ord = ord,
  steps = steps,
  outputname = "Browser_files/shinny_info.txt",
  info = info,
  lower_limit = lower_limit,
  upper_limit = upper_limit
)

###################################################################################################################
############################################# GET SIMILAR COLORS ###########################################

# With the graph we obtain the layouts that we are going to use
if (!is.null(random_seed)) {
  set.seed(random_seed)
}
lay_3D <- layout_with_fr(graph, dim = 3)

# We first look at which step are the largest number of communities
# Find number of communities by step
num_comm <- apply(im, 1, contar_unicos)

# Find the step with the largest number of communities
posiciones <- which(num_comm == max(num_comm))

# Select the smallest resolution to have more individuals per community
posicion_final <- min(posiciones)
step = posicion_final

# We get the graph of individuals by comm
# Adds colors by community
# For the chosen step:
graph_by_comm <-
  graph_louvain_by_comm(ord,
                        im,
                        step,
                        graph,
                        mapbig,
                        use_leiden = use_leiden)
# We obtain the average position of individuals per community in 3D
means_table_3D <- media_graph_3D(graph_by_comm, lay_3D)

# We order the communities from smallest to largest in means_table_3D
orden <- order(as.numeric(means_table_3D[, 1]))
means_table_3D_ord <- means_table_3D[orden,]

# We check if communities are missing

comm_faltantes <- length(unique(as.vector(as.matrix(im)))) - max(num_comm)
if (comm_faltantes != 0) {
  min_im <- min(im)
  B <- logical(length(unique(as.vector(as.matrix(im)))))
  if (min_im == 2) {
  B[as.numeric(means_table_3D_ord[, 1]) - 1] <- TRUE
  } else {
  B[as.numeric(means_table_3D_ord[, 1])] <- TRUE
  }
  for (i in sort(unique(as.vector(as.matrix(im))))) {
  # Adjust index based on minimum value if needed
  if (min_im == 2){
  index <- i-1
  }else{
  index <- i
  }
  if (B[i] != TRUE){
      # For the community that is missing
      posiciones_i <- which(im == i, arr.ind = TRUE)
      # Select the smallest position
      posicion_final_i <- min(posiciones_i[, 1])
      # graph_by_comm for the i community
      graph_by_comm_i <-
        graph_louvain_by_comm(
          ord,
          im,
          posicion_final_i,
          graph,
          mapbig,
          use_leiden = use_leiden
        )
      #Medias para la comunidad i
      means_table_3D_i <- media_graph_3D(graph_by_comm_i, lay_3D)
      medias_i <-
        which(means_table_3D_i[, 1] == i, arr.ind = TRUE)
      means_table_3D_ord <-
        rbind(means_table_3D_ord, means_table_3D_i[medias_i, ])
    }
  }
}

# We organize communities from smallest to largest
orden_2 <- order(as.numeric(means_table_3D_ord[, 1]))
means_table_3D_ord <- means_table_3D_ord[orden_2,]

lay_3D_colors <-
  data2cielab(as.data.frame(means_table_3D_ord[, 2:4]))

# Comunnity 1 in "white"
if (min(im)==1){
lay_3D_colors[1, 2] = "#FFFFFF"
}


Community_detection_dir = file.path(path, "Community_detection")
if (!file.exists(Community_detection_dir)) {
  # If it does not exist, create the directory
  dir.create(Community_detection_dir)
  cat("'Community_detection' folder created successfully.\n")
}
plot_name = "/Community_detection/resolution_plot_similarC.png"

# Create the resolution plot
make_res_plot(lay_3D_colors, im, path, plot_name, lower_limit, upper_limit)

# Colors for shiny app
Browser_files_dir = file.path(path, "Browser_files")
if (!file.exists(Browser_files_dir)) {
  # If it does not exist, create the directory
  dir.create(Browser_files_dir)
  cat("'Browser_files' folder created successfully.\n")
}
write.table(
  lay_3D_colors$V2,
  file.path(path, "Browser_files/Similar_colors.txt"),
  sep = ",",
  quote = F,
  row.names = F,
  col.names = F
)


##############################################################################################
################## Plot: Network by comm: normal colors / similar colors #####################

############# COLORS ##############
## Normal colors: mapbig

## Similar colors
colors_cielab <-
  matrix(0, nrow = length(lay_3D_colors[, 2]), ncol = 1)
colors_cielab[, 1] <- lay_3D_colors[, 2]

## Listas
graph_distintive_list <- list()
graph_similar_list <- list()

NetworksIndividuals_dir = file.path(path, "NetworksIndividuals")
if (!file.exists(NetworksIndividuals_dir)) {
  # If it does not exist, create the directory
  dir.create(NetworksIndividuals_dir)
  cat("'NetworksIndividuals' folder created successfully.\n")
}

# One legend to rule them all
legendcols = ceiling(length(unique(as.vector(im))) / 30)
# Distinctive colors
name_legend = paste0(
  "/NetworksIndividuals/Network_Individuals_Communities_DistinctiveColors_legend.png"
)
make_legend(mapbig, legendcols, val_width = 3500, val_height = 3500, val_pt.cex = 1.4, val_cex = 1.4)
# Similar colors
name_legend = paste0(
  "/NetworksIndividuals/Network_Individuals_Communities_SimilarColors_legend.png"
)
make_legend(colors_cielab, legendcols, val_width = 3500, val_height = 3500, val_pt.cex = 1.4, val_cex = 1.4)

# x 1:50 matching total of lambda values
for (x in 1:steps) {
  ######################################################################## Distinctive colors
  graph_by_comm <-
    plot_louvain_by_comm(
      ord,
      im,
      step = x,
      graph,
      mapbig,
      lay,
      name = paste0(
        "/NetworksIndividuals/Network_Individuals_Communities_DistinctiveColors_whiteComm_RIndex",
        x,
        ".png"
      ),
      path = path,
      use_leiden = use_leiden
    )
  
  # We remove community 1
  if (any(V(graph_by_comm)$carac == 1)) {
    vertices_to_remove <- V(graph_by_comm)[V(graph_by_comm)$carac == 1]
    lay_no_comm1 <- lay[-vertices_to_remove,]
    graph_by_comm_no1 <- delete_vertices(graph_by_comm, vertices_to_remove)
  } else{
    graph_by_comm_no1 <- graph_by_comm
    lay_no_comm1 <- lay
  }
  
  name_no_comm1 <- paste0(
    "/NetworksIndividuals/Network_Individuals_Communities_DistinctiveColors_RIndex",
    x,
    ".png"
  )
  plot_indv_no_comms1(path, name_no_comm1, graph_by_comm_no1, lay_no_comm1)
  
  # Make community network
  means_table <- media_graph(graph_by_comm, lay)
  freq_ordenado <- num_vert(graph_by_comm)
  tabla_freq_nonzeros <- dens_con(graph_by_comm)
  graph_only_comms <-
    graph_comms(graph_by_comm, tabla_freq_nonzeros, freq_ordenado)
  
  ## Plot
  Community_networks_dir = file.path(path, "Community_networks")
  if (!file.exists(Community_networks_dir)) {
    # If it does not exist, create the directory
    dir.create(Community_networks_dir)
    cat("'Community_networks' folder created successfully.\n")
  }
  plot_by_density(
    graph_only_comms,
    means_table,
    path,
    plotname = paste0(
      "/Community_networks/CommunityNetworks_2D_DistintiveC_whiteComm_RIndex",
      x ,
      ".png"
    )
  )
  
  # We remove community 1
  if (any(V(graph_only_comms)$carac == 1)) {
    graph_only_comms_no1 <- delete_vertices(graph_only_comms, 1)
  } else{
    graph_only_comms_no1 <- graph_only_comms
  }
  
  plot_by_density(
    graph_only_comms_no1,
    means_table,
    path,
    plotname = paste0(
      "/Community_networks/CommunityNetworks_2D_DistintiveC_RIndex",
      x ,
      ".png"
    )
  )
  
  # 3D Network
  graph_3d = graph_only_comms_no1
  # We normalize vertices and edges
  V(graph_3d)$size <-
    (rank(V(graph_3d)$size) / length(V(graph_3d)$size)) * 50
  E(graph_3d)$width <-
    (rank(E(graph_3d)$width) / length(E(graph_3d)$width)) * 20
  # 3D layout selected for communities in the graph
  coords <-
    means_table_3D_ord[as.numeric(V(graph_only_comms_no1)$name),]
  
  # We add the coordinates to a single object
  graph_3d$coords <- coords[,-1]
  
  graph_distintive_list[[x]] <- graph_3d
  
  # Delete objects
  rm(
    graph_by_comm,
    means_table,
    freq_ordenado,
    tabla_freq_nonzeros,
    graph_only_comms,
    graph_only_comms_no1,
    graph_3d,
    coords
  )
  
  ############################################################################ Plots with similar colors
  # Network of nodes colored by community:
  graph_by_comm <-
    plot_louvain_by_comm(
      ord,
      im,
      step = x,
      graph,
      colors_cielab,
      lay,
      name = paste0(
        "/NetworksIndividuals/Network_Individuals_Communities_whiteComm_SimilarColors_RIndex",
        x,
        ".png"
      ),
      path = path,
      use_leiden = use_leiden
    )
  
  # We remove community 1
  if (any(V(graph_by_comm)$carac == 1)) {
    vertices_to_remove <- V(graph_by_comm)[V(graph_by_comm)$carac == 1]
    lay_no_comm1 <- lay[-vertices_to_remove,]
    graph_by_comm_no1 <- delete_vertices(graph_by_comm, vertices_to_remove)
  } else{
    graph_by_comm_no1 <- graph_by_comm
    lay_no_comm1 <- lay
  }
  
  name_no_comm1 <- paste0(
    "/NetworksIndividuals/Network_Individuals_Communities_SimilarColors_RIndex",
    x,
    ".png"
  )
  plot_indv_no_comms1(path, name_no_comm1, graph_by_comm_no1, lay_no_comm1)
  
  # Make community Network
  means_table <- media_graph(graph_by_comm, lay)
  freq_ordenado <- num_vert(graph_by_comm)
  tabla_freq_nonzeros <- dens_con(graph_by_comm)
  graph_only_comms <-
    graph_comms(graph_by_comm, tabla_freq_nonzeros, freq_ordenado)
  
  ## Plot
  Community_networks_dir = file.path(path, "Community_networks")
  if (!file.exists(Community_networks_dir)) {
    # If it does not exist, create the directory
    dir.create(Community_networks_dir)
    cat("'Community_networks' folder created successfully.\n")
  }
  plot_by_density(
    graph_only_comms,
    means_table,
    path,
    plotname = paste0(
      "/Community_networks/CommunityNetworks_2D_SimilarC_whiteComm_RIndex",
      x ,
      ".png"
    )
  )
  
  # Remove Community 1
  if (any(V(graph_only_comms)$carac == 1)) {
    graph_only_comms_no1 <- delete_vertices(graph_only_comms, 1)
  } else{
    graph_only_comms_no1 <- graph_only_comms
  }
  
  plot_by_density(
    graph_only_comms_no1,
    means_table,
    path,
    plotname = paste0("/Community_networks/Community_Networks_2D_SimilarC_RIndex", x , ".png")
  )
  
  # 3D Network
  graph_3d = graph_only_comms_no1
  # We normalize vertices and edges
  V(graph_3d)$size <-
    (rank(V(graph_3d)$size) / length(V(graph_3d)$size)) * 50
  E(graph_3d)$width <-
    (rank(E(graph_3d)$width) / length(E(graph_3d)$width)) * 20
  # 3D layout selected for communities in the graph
  coords <-
    means_table_3D_ord[as.numeric(V(graph_only_comms_no1)$name),]
  
  # We add the coordinates to a single object
  graph_3d$coords <- coords[,-1]
  
  graph_similar_list[[x]] <- graph_3d
  
  # Delete objects
  rm(
    graph_by_comm,
    means_table,
    freq_ordenado,
    tabla_freq_nonzeros,
    graph_only_comms,
    graph_only_comms_no1,
    graph_3d,
    coords
  )
}

Community_networks_dir = file.path(path, "Community_networks")
if (!file.exists(Community_networks_dir)) {
  # If it does not exist, create the directory
  dir.create(Community_networks_dir)
  cat("'Community_networks' folder created successfully.\n")
}
saveRDS(
  graph_distintive_list,
  file.path(Community_networks_dir, "3Dplots_distintive.rds")
)
saveRDS(graph_similar_list,
        file.path(Community_networks_dir, "3Dplots_similar.rds"))

save(graph, res_mat, res_colors , lay_3D_colors, info_obj, file = paste(path, "Summary_objects.RData", sep = "/"))
writeLines(capture.output(sessionInfo()), "session_info.txt")

