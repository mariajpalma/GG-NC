#### Explore the range
#### This auxiliar script makes a initial exploration
#### of the communities using louvain or leide
#### to give you a logarithmic space while keeping a defined % of your data
#### in communities bigger than a set mininum. 

## load libraries

library(igraph)
library(ggplot2)
library(pracma)


### load parameters
kind <- readline(prompt = "Enter the kind of data you are using (PCA/GRM/IBD): ")
path <- readline(prompt = "Enter the path where your files are found:")
file_name <- readline(prompt = "Enter your file name or the base name in the case of GRM: ")
percentage <- readline(prompt = "Enter the percentage of data you want to keep in communities bigger thant the min: ")
min_size <- readline(prompt = "Enter the minimum size of the communities: ")
algorithm <- readline(prompt = "Enter 1 if you want to use Leiden algorithm if not use 0 for Louvain: ")
seed <- readline(prompt = "Enter a number to use as seed: ")
steps <- readline(prompt = "Enter the number of steps used to explore: ")
metadata <- readline(prompt = "Enter the your file name for metadata, not used, but required: ")
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

ibd_network <- function (path, fname1, fname2, max, prune) {
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
  relatives_1 <- ibd[which(ibd$IBD_CM_SUM >= max),]$ID1
  relatives_2 <- ibd[which(ibd$IBD_CM_SUM >= max),]$ID2
  relatives <- c(relatives_1, relatives_2)
  # Calculate frequency of each individual
  freq <- table(relatives)
  # Identify hubs (individuals with IBD connections above threshold)
  hubs <- names(freq[freq > 1]) 
  # Remove hubs from IBD data
  ibd <- ibd[!(ibd$ID1 %in% hubs), ]
  # Select unrelated individuals based on IBD threshold
  Nohub <- ibd[which(ibd$IBD_CM_SUM >= max),]$ID1
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
    set.seed(3)
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

if(kind=="GRM" || kind=="PCA" )
{
  m_parameter <-0
}else
{
  m_parameter <-1000000
}

graph <- input(kind = kind, path = path, data = file_name, info = metadata, max = m_parameter, prune = 1)[[1]]

### Explore values for Louvain or Leiden algorithm

min_resolution <- 1e-6  # Lowest resolution (10^-6)
max_resolution <- 10     # Highest resolution (10^1)
num_steps <- 8           # Number of steps (powers of 10)
log_min <- logspace(log10(min_resolution), log10(max_resolution), num_steps)

n_comm <- rep(0, num_steps)
stability_found <- FALSE
a <- 1  # Start from the highest resolution and go down


# Loop to explore from higher to lower resolutions (descending)
for (i in log_min) {
  cl <- if (as.integer(algorithm)) {
    set.seed(as.integer(seed))
    cluster_leiden(graph, resolution_parameter = i)
    
  } else {
    set.seed(as.integer(seed))
    cluster_louvain(graph, resolution = i)
  }
  
  current_comm <- max(cl$membership)  # Current number of communities
  n_comm[a] <- current_comm
  
  a <- a +1
}

rle_values <- rle(n_comm)

if (rle_values$lengths[1]==1) {
  if (rle_values$values[1]==1){
    selected_min <- log_min[1]
    message(paste("Stability was found at "), log10(log_min[1]), "con ", rle_values$values[1], " comunidad" )
  }
  else{
    stop("No stability found in the number of communities across the tested resolutions.")
  }
  
}else
{
  selected_min <- log_min[rle_values$lengths[1]]
  message(paste("Stability was found at "), log10(log_min[rle_values$lengths[1]]), " con ", rle_values$values[1], " comunidades" )
}



# Adaptive upper limit based on the selected minimum resolution
max_res <- log10(selected_min) + 10  # Adjusted to allow a wider range for exploration
exploration_space <- logspace(log10(selected_min), max_res, as.integer(steps))

# Explore higher resolutions to find upper limit based on community sizes
flag <- 0
a <- 1

while (flag != 1) {
  if (a > length(exploration_space)) {
    stop(paste("The upper limit exceeds the range for resolution starting from", round(log10(selected_min))))
  }
  
  i <- exploration_space[a]
  cl <- if (as.integer(algorithm)) {
    set.seed(as.integer(seed))
    cluster_leiden(graph, resolution_parameter = i)
  } else {
    set.seed(as.integer(seed))
    cluster_louvain(graph, resolution = i)
  }
 
  sizes <- table(cl$membership)
  p <- sum(sizes[names(which(sizes > as.integer(min_size)))]) / sum(sizes) ## individuals in communities > min
  
  if (p * 100 < as.integer(percentage)) {
    flag <- 1
    message(paste("The upper limit that meets your requirements is near ", log10(exploration_space[a-1]), " we suggest ", round(log10(exploration_space[a-1]))))
  }
  
  a <- a + 1
}



