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
ibd_network <- function (path, fname1, fname2, max, prune){
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
  # Write hubs to a text file for information
  write.table(
    hubs,
    file = file.path(path, "output/hubs.txt"),
    sep = ",",
    row.names = F,
    col.names = F
  )
  # Remove hubs from IBD data
  ibd <- ibd[!(ibd$ID1 %in% hubs), ]
  # Select unrelated individuals based on IBD threshold
  Nohub <- ibd[which(ibd$IBD_CM_SUM >= max),]$ID1
  # Write unrelated individuals to a text file
  write.table(
    Nohub,
    file = file.path(path, "output/Nohub.txt"),
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
  # Extract unique genetic regions and corresponding colors
  info_nodup <-
    info[!duplicated(info$genetic_region, fromLast = FALSE), ]
  spop_color <- cbind(info_nodup$genetic_region, info_nodup$color)
  colnames(spop_color) <- c("genetic_region", "color")
  # Create an igraph from the unrelated data
  ibd_graph <-
    graph_from_data_frame(unrelated_df, directed = F, vertices = id_ordenado)
  print(length(V(ibd_graph)))
  # Set colors of vertices in the graph
  V(ibd_graph)$color <- V(ibd_graph)$color
   # Set widths of edges based on weight
  E(ibd_graph)$width <-
    (E(ibd_graph)$weight / max(E(ibd_graph)$weight))*50
  # Prune the network if requested
  if (prune)
  {
    ibd_graph <- prune_network(ibd_graph)
  }
  # Return the graph and the colors by superpopulation per genetic region
  return(list(ibd_graph, spop_color))
}

prune_network <- function (net){
  while (sum(degree(net) == 1)> 0){
    net <- delete_vertices(net, names(which(degree(net) == 1)))
    net <- delete_vertices(net, names(which(degree(net) == 0)))
  }
  
  clusters <- cluster_louvain(net, resolution = 0.01)
  print(table(clusters$membership))
  to_remove <- names(which(table(clusters$membership)< 25))
  print(to_remove)
  ind_list <- c()
  for (i in to_remove){
    comm <- which(clusters$membership==i)
    ind_list <- c(ind_list,clusters$names[comm])
    net <- delete_vertices(net, clusters$names[comm])
  }
  print(ind_list)
  return(net)
}

# This function reads PCA data and additional information about individuals (such as population labels) from files.
# It then calculates a correlation matrix, creates an adjacency matrix, and constructs a network graph based on
# a specified threshold. Population labels and colors are assigned to vertices in the graph, and the function returns
# the igraph object along with information about genetic regions and their colors.
pca_network <- function (path, fname1, fname2, max) {
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
  vert[,1] <- V(net_pca)$name
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
  # Read the GRM binary files
  datafile = file.path(path, fname1)
  GRM = ReadGRMBin(datafile)
  # Extract the diagonal and upper triangle of the matrix
  diagonal <- GRM[[1]]  # Diagonal elements of the GRM
  triangulo_superior <- GRM[[2]]  # Upper triangle elements of the GRM
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
  # Returns the igraph object and the the genetic regions and their corresponding colors.
  return(list(net_grm, spop_color))
}
# This function makes an igraph object depending on the kind and values of the input data.
input <- function(kind, path, data, info, max, prune) {
  # Verify if the output directory already exists
  directorio = file.path(path,"output")
  if (!file.exists(directorio)) {
    # If it does not exist, create the directory
    dir.create(directorio)
    cat("Carpeta 'output' creada correctamente.\n")
  } else {
    cat("La carpeta 'output' ya existe.\n")
  }
  # Determine the type of network to create
  if(kind == "IBD") {
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
          "'IBD' (Identical by descent), 'PCA'(Principal Component Analysis) or 'GRM' (Genetic Relationship Matrix)."
        )
      }
    }
  }
  # Returns a list with the igraph object and the genetic regions and their corresponding colors.
  return(network)
}
# Function to perform clustering over a range of resolutions and save the results
get_lambda_results <- function(input, steps, path, name, lower_limit, upper_limit) {
  # An empty list to store the clustering results for each resolution.
  list_cl <- list()
  # A numeric vector to store the number of communities for each resolution.
  ncom <- 0
  a <- 1
  # Iterates over a sequence of resolution parameters
  for (i in logspace(lower_limit, upper_limit, n = steps)) {
    cat(a, "\n")
    # The cluster_louvain function is called with the input data and the current resolution to perform clustering.
    cl <- cluster_louvain(input, resolution = i)
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
  # The function returns the list list_cl, which contains the clustering results for each resolution.
  return(list_cl)
}
# nodebased, relabel and pollock ar a set of functions that were developed by Dr. Lewis for the research paper titled
# "The function of communities in protein interaction networks at multiple scales."
# The functions have been adapted from Matlab to R by our research group.
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
  for (i in 1:(N - 1)) {
    cat(i, "\n")
    S2 <- Sloop[N - i,]
    currentmax <- max(Ssort)
    Snew <- relabel(i, S1, S2, currentmax)
    Ssort[N - i,] <- Snew
    S1 <- Snew
  }
  return(Ssort)
}
# This `relabel` function takes several arguments: `num`, `S1`, `S2`, and `currentmax`.
# The purpose of this function is to re-label communities in `S2` based on the mapping provided by `S1`,
# while ensuring that the re-labeled communities retain their identity.
relabel <- function(num, S1, S2, currentmax){
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
         cat("SUPER ERROR")
       }
     }
      if (m > filas_l2) {
        currentmax= currentmax + noverlap
        unassigned = sort(setdiff(coms2, assigned[, 2])) # communities without labels
        for (k in 1:length(unassigned)) {
          assigned = rbind(assigned, c(currentmax + k, unassigned[k]))
        }
      } else { # If `m` is not greater than the number of rows in `l2`, proceed with the re-labeling process.
        if (is.matrix(l2)) {
          consider = l2[m, ]
        } else{
          if (is.vector(l2)) {
            consider = l2
          } else{
            cat("SUPER ERROR")
          }
        }
        if ((length(intersect(consider[3], assigned[, 2])) == 0) &&
            (length(intersect(consider[2], assigned[, 1])) == 0)) {
          if (consider[1]==0){
            noverlap=noverlap +1
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
pollock <- function(Ssort, little, r, path) {
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
    if (sum(sizescom > little) > 1) {
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
  a <- 1:length(rep)
  im <- Ssort2
  # Replace the community labels in `Ssort2` with their corresponding indices in `a`.
  for (i in 1:length(rep)) {
    z <- which(Ssort2 == rep[i])
    im[z] <- a[i]
  }
  # Number of communities
  m <- max(im)
  # Select colors for communities.
  colores <- distinct_colors(
    (max(im) - 1),
    minimal_saturation = 33,
    minimal_lightness = 20,
    maximal_lightness = 80
  )
  # Obtain colors in hex
  map <- colores$name
  # Mix the colors
  map <- map[sample(length(map))]
  # Store colors in vector ma
  ma <- vector("character", length = m)
  ma[2:m] <- map
   # Makes first one white. Little coms all have label 0, so will appear white.
  ma[1] <-
    c("#FFFFFF")
  im2 <- apply(im, 2, rev)
  # Axis X values
  x_vals <- 1:ncol(Ssort)
  # Axis Y values
  y_vals <- 1:nrow(Ssort)
  # Store the plot as a .png image.
  png(file= paste(path,"res_plot.png",sep="/"))
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
    las = 1
  )
  mtext(
    side = 2,
    line = 2,
    text = "Resolution value",
    cex = 1,
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
    las = 1
  )
  dev.off()
  # Returns a list with the individuals communities by level or resolution, the individal orders and the colors by community.
  return(list(im, order, ma))
}

# Function to plot Louvain clustering results colored by superpopulation
plot_louvain_by_spop <- function(graph, cl_list, spop_color, R, lay, name) {
  # Extract the clustering results for the given resolution
  cl <- cl_list[[as.character(R)]]
  # Open an SVG device for plotting
  svg(file= paste(path,name,sep="/"), height = 20, width = 30)
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
overlap_hm <- function(path, fname2, cl_list, R, spop_color,name){
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
    cbind(id_ordenado,cl$membership)
  idv_comm_pop_spop <- idv_comm_pop_spop[,c(1,2,3,5,4)]
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
  # Generate the heatmap
  #https://jokergoo.github.io/2020/05/06/translate-from-pheatmap-to-complexheatmap/
  png(file= paste(path,name,sep="/"),width=1000, height=1000)
  HM <- ComplexHeatmap::pheatmap(
    t(mat_overlap),
    name = "Overlap",
    annotation_row = annotation_row_spop,
    annotation_colors = ann_colors_spop,
    display_numbers = TRUE,
    number_color = "black",
    fontsize_number = 5,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    border_color = "white",
    row_split = annotation_row_spop$Spop,
    legend = TRUE
  )
  draw(HM)
  dev.off()
}

# Function to calculate stability metrics for Louvain clustering results

stability_matrix_metrics <- function (network, R, metric, path)
{
  mpair <- matrix(nrow = 100, ncol = 100)
  obj <- matrix(nrow=length(V(network)),ncol = 100)
  cl <-  cluster_louvain(network, resolution = R)
  obj[,1] <- cl$membership
  for (i in 2:100)
  {
    cl <- cluster_louvain(network, resolution = R)
    obj[,i] <- cl$membership
  }
  for (i in 1:100)
  {
    for (j in i:100)
    {
      if (i==j)
      {
       	mpair[i,j] <- NA
      }
      else
      {
       	if (metric=="NID")
        {
          mpair[i,j] <-NID(obj[,i],obj[,j])
        }
	if (metric=="ARI")
        {
          mpair[i,j] <- compare(obj[,i],obj[,j], "adjusted.rand")
        }
      }
    }
  }
    ind_mem <- matrix(nrow = length(V(network)), ncol = length(V(network)), data = 0)
    for (i in 1:100)
    {
      for (j in 1:(length(V(network))-1))
      {
       	for (k in j:length(V(network)))
        {
          if (obj[j,i]==obj[k,i]){
            ind_mem[j,k]=ind_mem[j,k]+1
            ind_mem[k,j]=ind_mem[k,j]+1
          }
	}
      }
    }
    rownames(ind_mem) <-cl$names
    colnames(ind_mem) <-cl$names
    diag(ind_mem) <- diag(ind_mem)/2
    write.table(ind_mem, quote = FALSE, row.names = TRUE, col.names = TRUE, file=paste(path,paste("Common_membership_file_R",R, sep=""),sep="/"))
    return(mpair)
}

# Function to compute and aggregate stability metrics for various resolution (lambda) values
df_metrics <- function(graph, metric_name, steps, path, lower_limit, upper_limit) {
  # Initialize a data frame to store metric values and corresponding lambda values
  DF_lambda <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(DF_lambda) <- c(toupper(metric_name), "lambda")
  
  # Define the number of cores to use for parallel processing
  num_cores <- detectCores()
  cluster <- makeCluster(num_cores)
  registerDoParallel(cluster)
  
  # Define the range of lambda values
  lambda_values <- logspace(lower_limit, upper_limit, n = steps)
  
  # Define the function to be executed in parallel
  calculate_stability <- function(lambda) {
    # Compute the stability matrix for the given lambda and metric
    smn <- stability_matrix_metrics(graph, R = lambda,toupper(metric_name), path = path)
    # Convert the upper triangular part of the stability matrix to a data frame
    tempo <- as.data.frame(smn[upper.tri(smn, diag = FALSE)])
    tempo$lambda <- lambda # Add the lambda value as a column
    colnames(tempo) <- c(toupper(metric_name), "lambda")
    return(tempo)
  }
  
  # Use parallel processing to calculate stability metrics for each lambda value
  DF_lambda <- foreach(lambda = lambda_values, .combine = rbind, .export = "stability_matrix_metrics", .packages = c("igraph", "aricode", "pracma" )) %dopar% {
    calculate_stability(lambda)
  }
  # Stop the parallel cluster
  stopCluster(cluster)
  # A data frame containing the stability metric values and corresponding lambda values for the specified range of lambda values.
  # Each row represents a pair of clusterings compared by the specified metric, along with the lambda value used for clustering.
  return(DF_lambda)
}
# Function to plot and save metric analysis results
plots_metric <- function(DF_lambda, metric_name, path, name){
  ### Boxplot
  colnames(DF_lambda) <- c("metric", "lambda")
  box <- ggplot(DF_lambda, aes( x= lambda, y=metric, group=lambda)) + 
    geom_boxplot() + scale_x_continuous(trans='log10') + ggtitle(metric_name)
  ggsave(filename = paste(path,paste0(name,"_boxplot.png"),sep="/"), plot = box, width = 13, height = 7, units = "in",
         bg = "white",
         dpi = 300)
  ####  Mean plot
  DF_lambda_m <- DF_lambda %>% group_by(lambda) %>% 
    summarise(mean=mean(metric))
  # Convert tibble to df
  DF_lambda_m  <- DF_lambda_m  %>% as.data.frame()
  mean_plot <- ggplot(DF_lambda_m, aes( x= lambda, y=mean)) + 
    geom_point() + scale_x_continuous(trans='log10') + ggtitle(paste0(metric_name,"_mean"))
  ggsave(filename = paste(path,paste0(name,"_mean.png"),sep="/"), plot = mean_plot, width = 13, height = 7, units = "in",
         bg = "white",
         dpi = 300)
  #### Median plot
  DF_lambda_m <- DF_lambda %>% group_by(lambda) %>% 
    summarise(median=median(metric))
  # Convert tibble to df
  DF_lambda_m  <- DF_lambda_m  %>% as.data.frame()
  median_plot <- ggplot(DF_lambda_m, aes( x= lambda, y=median)) + 
    geom_point() + scale_x_continuous(trans='log10') + ggtitle(paste0(metric_name,"_median"))
  ggsave(filename = paste(path,paste0(name,"_median.png"),sep="/"), plot = median_plot, width = 13, height = 7, units = "in",
         bg = "white",
         dpi = 300)
  #### Variance plot
  DF_lambda_m <- DF_lambda %>% group_by(lambda) %>% 
    summarise(variance=var(metric))
  DF_lambda_m <- DF_lambda_m  %>% as.data.frame()
  variance_plot <- ggplot(DF_lambda_m, aes( x= lambda, y=variance)) + 
    geom_point() + scale_x_continuous(trans='log10') + ggtitle(paste0(metric_name,"_variance"))
  ggsave(filename = paste(path,paste0(name,"_variance.png"),sep="/"), plot=variance_plot, width = 13, height = 7, units = "in",
         bg = "white",
         dpi = 300)
}

# Shiny App: Interactive Map input 
# outputname es el nombre del archivo final que entrará a la carpeta output/
info_map <- function(path, im, cl, ord, info, steps, outputname, lower_limit, upper_limit) {
  # Load the metadata of the individuals
  meta_file = file.path(path, info)
  metadata = read.csv(meta_file, sep = "\t", header = T)
  # Leave only individuals and their pop that are in cl$names
  # Remember that they enter the nodebases function in order of cl$names
  # and that in Pollock they are assigned a different place to improve the
  # visualization and that order comes in the order vector that comes out of Pollock
  met_indv_pop <- metadata[,c(1, 4,3)]
  filtrado <-
    met_indv_pop[met_indv_pop$indv %in% as.character(cl$names), ]
  # Sort like Pollock using cl$names and ord
  cl_names_ord <- cl$names[ord[, 1]]
  id_order <- match(filtrado[, 1], as.character(cl_names_ord))
  id_ordenado <- filtrado[order(id_order),]
  # Create the table for communities per individual with their info
  Lambda <-
    rep(logspace(lower_limit, upper_limit, n = steps), times = nrow(id_ordenado))
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
  prop_by_row_file = file.path(path, "output/comms_prop.txt")
  write.table(
    prop_by_row,
    prop_by_row_file,
    sep = "\t",
    quote = F,
    row.names = F,
    col.names = T
  )
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
    cbind(pp_merge$population,
          pp_merge$pop3code,
          pp_merge$genetic_region,
          pp_merge$poject,
          pp_merge$latitude,
          pp_merge$longitude,
          pp_merge$Lambda)
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
  colnames(prop) <- paste0("C", 1:max(im))
  # Convert zeros to NA
  prop[prop == 0] <- "NA"
  final_table <- cbind(pp_merge2, prop)
  # Sort by lambda and pop
  indice_orden_2 <-
    order(final_table[, 1], as.numeric(final_table[, 6]))
  final_table_order <- final_table[indice_orden_2,]
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
}

# Heatmap: by communities
overlap_by_comm <- function(path, info, graph, im, ord, step, name){
  # Take the individuals from cl$name
  cl <- cluster_louvain(graph, r = 1) #No importa el valor de r
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
    comm[count,1] = im[step, j]
    count = count + 1
  }
  # Unite id with communities
  comm <- paste0("C",comm)
  idv_comm_pop_spop <-
    cbind(id_ordenado,comm)
  idv_comm_pop_spop <- idv_comm_pop_spop[,c(1,2,3,5,4)]
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
  spop_color <- idv_comm_pop_spop[!duplicated(idv_comm_pop_spop[,5], fromLast = FALSE),]
  spop_color <- spop_color[,c(3,5)]
  Spop_colores = spop_color[, 2]
  names(Spop_colores) = spop_color[, 1]
  ann_colors_spop = list(Spop = Spop_colores)
  mat_overlap <- as.data.frame.matrix(tabla_overlap)
  png(file= paste(path,name,sep="/"),width=1000, height=1000)
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
plot_louvain_by_comm <- function(ord, im, step, graph, mapbig, lay, name, path) {
  # Here step refers to the row of im from which we will get the communities
  # Sort like Pollock using cl$names and ord
  cl <- cluster_louvain(graph, resolution = 1) # Any r is ok, we just need names
  cl_names <- cl$names
  cl_names_ord <- cl_names[ord[,1]]
  # Put together id and its comm
  imcols = ncol(im)
  id_comm <- matrix(data = 0,
                    ncol = 2,
                    nrow = ncol(im))
  count = 1
  for (j in 1:imcols) {
    id_comm[count,2] = im[step, j]
    count = count + 1
  }
  id_comm[,1] <- cl_names_ord
  colnames(id_comm) <- c("ID","Comm")
  comms <- sort(unique(as.numeric(id_comm[,2])))
  # Get the colors of those communities
  colors <- mapbig[comms,1]
  comms_colors <- cbind(comms,colors)
  # Order the communities with their color from smallest to largest
  colnames(comms_colors) <- c("Comm","Colors")
  # Create a matrix of individuals with their community and color
  id_colors <- merge(id_comm,comms_colors, by = "Comm")
  id_colors2 <- cbind(id_colors[,2],id_colors[,1],id_colors[,3])
  colnames(id_colors2) <- c("ID","Comm","Colors")
  # Organize the information with the individuals in a graph
  vert <- V(graph)$name
  id_order <- match(id_colors2[,1], vert)
  id_ordenado <- id_colors2[order(id_order),]
  # Assign the community and color attributes to the graph
  V(graph)$carac <- id_ordenado[,2]
  V(graph)$color <- id_ordenado[,3]
  # Graph
  svg(file= paste(path,name,sep="/"), width = 30, height = 20)
  plot(
    graph,
    layout = lay,
    vertex.size = 6,
    vertex.label = NA,
    edge.color = "gray",
    col = V(graph)$color
  )
  legend(
    1.1,
    1.2,
    legend = comms_colors[,1] ,
    bty = "n",
    pch = 20,
    pt.cex = 3,
    cex = 3 ,
    horiz = FALSE,
    inset = c(0.1, 0.1),
    col = comms_colors[, 2],
  )
  dev.off()
  return(graph)
}

# Here we need the graph that is generated in
# plot_louvain_by_comm()
# Get average position of x and y from lay
# vertex_attr_names(graph) helps see what attributes are there
media_graph <- function(graph_by_comm, lay) {
  pos <- cbind(V(graph_by_comm)$carac, lay)
  pos_x <- aggregate(as.numeric(pos[, 2]) ~ pos[, 1],
                     data = pos,
                     FUN = mean)
  pos_y <- aggregate(as.numeric(pos[, 3]) ~ pos[, 1],
                     data = pos,
                     FUN = mean)
  means_table <- cbind(pos_x, pos_y[, 2])
  
  return(means_table)
}
# Number of nodes per community
num_vert <- function(graph_by_comm) {
  # Get the frequency table
  tabla_frecuencias <- table(V(graph_by_comm)$carac)  
  # Convert the frequency table into a data frame
  df_frecuencias <- as.data.frame(tabla_frecuencias)
  df_frecuencias <- as.matrix(df_frecuencias)
  # Rename columns
  colnames(df_frecuencias) <- c("Comm", "Frecuencia")
  ord <- order(as.numeric(df_frecuencias[,1]))
  df_frecuencias_ord <- df_frecuencias[ord,]
  # Make freq_ordenado a matrix
  freq_ordenado <- matrix(0, nrow = length(ord), ncol = 2)
  for (i in 1:length(ord)) {
    freq_ordenado[i,1] = df_frecuencias[i,1]
    freq_ordenado[i,2] = df_frecuencias[i,2]
  }
  return(freq_ordenado)
}

# Connection density
# Attributes: edge_attr_names(graph)
dens_con <- function(graph_by_comm) {
  # Obtain pairs of individuals from the graph and name
  conections <- get.edgelist(graph_by_comm)
  colnames(conections) <- c("ID1", "ID2")
  # Obtain a list of individuals and the community to which they belong
  indv_comm <- cbind(V(graph_by_comm)$name, V(graph_by_comm)$carac)
  colnames(indv_comm) <- c("ID1", "Comm")
  # Get the communities from ID1
  merge_1 <- merge(conections, indv_comm, by = "ID1")
  colnames(merge_1) <- c("ID1", "ID2","Comm1")
  # We change name to obtain communities for ID2
  colnames(indv_comm) <- c("ID2", "Comm")
  merge_2 <- merge(merge_1, indv_comm, by = "ID2")
  # Sort and name the matrix
  merged <- merge_2[,c(2,1,3,4)]
  colnames(merged) <- c("ID1", "ID2","Comm1","Comm2")
  # Now we are only left with Com1 and Com2
  comms <- cbind(merged$Comm1,merged$Comm2)
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
  vector1 <- relaciones_ordenadas[,1]
  vector2 <- relaciones_ordenadas[,2]
  combined_levels <- unique(c(vector1, vector2))
  freq_table <- table(factor(vector1, levels = combined_levels), factor(vector2, levels = combined_levels))
  # Shows the frequency table
  # There are connections within communities
  # We don't need those and we make them 0
  for(i in 1:ncol(freq_table)){
    freq_table[i,i] = 0
  }
  tabla_freq2 <- as.data.frame(freq_table)
  # In the frequency table all against all are integrated
  # So there are many zeros left. We remove them.
  if(nrow(tabla_freq2) > 1) {
    nonzeros2 <- which(tabla_freq2[, 3] != 0)
    tabla_freq_nonzeros2 <- tabla_freq2[nonzeros2, ]
    tabla_freq_nonzeros2 <- as.matrix(tabla_freq_nonzeros2)
    colnames(tabla_freq_nonzeros2) <- c("from", "to", "weight")
  }else{
    tabla_freq_nonzeros2 <- tabla_freq2
    tabla_freq_nonzeros2 <- as.matrix(tabla_freq_nonzeros2)
    colnames(tabla_freq_nonzeros2) <- c("from", "to", "weight")
  }
  return(tabla_freq_nonzeros2)
}

# Make the graph and plot
graph_comms <- function(graph_by_comm, tabla_freq_nonzeros, freq_ordenado) {
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
  V(graph_only_comms)$carac <- vert
  V(graph_only_comms)$size <- as.numeric(freq_ordenado[, 2])
  V(graph_only_comms)$color <- color
  E(graph_only_comms)$weight <- as.numeric(E(graph_only_comms)$weight)
  E(graph_only_comms)$width <- E(graph_only_comms)$weight
  
  return(graph_only_comms)
}

# Function to count unique numbers per row
contar_unicos <- function(fila) {
  longitud <- length(unique(fila))
  return(longitud)
}
# 3D layout tights
media_graph_3D <- function(graph_by_comm, lay_3D) {
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
  
  return(means_table_3D)
}

# The plot_by_density function takes as input a graph with communities (without community 1),
# a table of means, a file path, and a name for the output file.
# The function sorts the nodes in the graph based on the mean table,
# adjusts the size of the nodes and the width of the edges based on their relative ranking,
# and then graphs and saves the resulting display.
plot_by_density <- function(graph_only_comms_no1, means_table, path, plotname) {
  # Sort means_table
  comm_order <- match(as.numeric(V(graph_only_comms_no1)$carac),as.numeric(means_table[, 1]))
  means_table_ord_filt <- means_table[comm_order,]
  lay_mean <- as.matrix(means_table_ord_filt[,c(2,3)])
  # Adjusts the size of the nodes and the width of the edges based on their relative ranking
  vsize = (rank(V(graph_only_comms_no1)$size) / length(V(graph_only_comms_no1)$size)) *
    50
  ewidth = (rank(E(graph_only_comms_no1)$width) / length(E(graph_only_comms_no1)$width)) *
    20
  if(length(rank(E(graph_only_comms_no1)$width)) == 1){
    ecolor = NA
  }else{ecolor = "gray"}
  # Graph and save
  png(file.path(path, plotname))
  plot(
    graph_only_comms_no1,
    layout = lay_mean,
    vertex.size = vsize,
    edge.width = ewidth, 
    edge.color = ecolor,
    col = V(graph_only_comms_no1)$color,
  )
  dev.off()
}
# The function takes a network as input and iteratively removes vertices
# with specific characteristics to clean up the network
prune_network <- function (net){
  # Removes vertices with a degree of 1 or 0.
  while (sum(degree(net) == 1)> 0){
    net <- delete_vertices(net, names(which(degree(net) == 1)))
    net <- delete_vertices(net, names(which(degree(net) == 0)))
  }
  # Removes small communities based on a specified size threshold.
  clusters <- cluster_louvain(net, resolution = 0.01)
  to_remove <- names(which(table(clusters$membership)< 25))
  for (i in to_remove){
    comm <- which(clusters$membership==i)
    net <- delete_vertices(net, clusters$names[comm])
  }
  return(net)
}

###############################################################################
################################SCRIPT#########################################
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
kind = args[1]
path = args[2]
data = args[3]
info = args[4]
max = as.numeric(args[5])
steps = as.numeric(args[6])
#Lambda = as.numeric(args[7]) #Resolución para cl
Lambda = args[7]
Prune = as.numeric(args[8])
min_comms = as.numeric(args[9])
lower_limit = as.numeric(args[10])
upper_limit = as.numeric(args[11])
r = 0 #r is which layer (counting from bottom) will be perfectly sorted



# Creating graph 
network <- input(kind, path, data, info, max, Prune)
graph <- network[[1]]
spop_color <- network[[2]]
lay <- layout_with_fr(graph)
tic()
cl_list <- get_lambda_results(graph, steps, path, "output/M_IBD_50.txt", lower_limit=lower_limit, upper_limit=upper_limit)
toc()

set.seed(2)
lay=layout_with_fr(graph)
svg(paste(path,"Graph_Simple.svg",sep="/"))
plot(graph, vertex.size = 3, layout=lay, vertex.label=NA)
dev.off()

# Read Lambda values 
lambda_path = file.path(path, Lambda)
lambda_values <- as.vector(unlist(read.table(lambda_path, sep = ",", header = F)))

# Plot: Stability plot
Mfile <- file.path(path, "output/M_IBD_50.txt")
M <- as.matrix(read.csv(Mfile, sep = ",", header = FALSE))
tic()
res = nodebased(M)
toc()

resfile = file.path(path, "output/res_50.txt")
write.table(
  as.data.frame(res),
  resfile,
  sep = ",",
  quote = F,
  row.names = FALSE,
  col.names = FALSE
)

### Pollock 
resfile = file.path(path, "output/res_50.txt")
Ssort <- as.matrix(read.csv(resfile, sep = ",", header = FALSE))

res_list <- pollock(Ssort, min_comms, r, path)
res_mat <- res_list[[1]]
res_indv <- res_list[[2]]
res_colors <- res_list[[3]]

plots_metric(df_metrics(graph, "nid", steps = steps, path = path, lower_limit=lower_limit, upper_limit=upper_limit), metric_name = "NID",
             path = path, name = "TEST2_NID")

for (lambda_index in lambda_values){
  # Get lambda value 
  Lambda <- logspace(lower_limit,upper_limit,steps)[lambda_index]
  
  # Plot: Network 
  plot_louvain_by_spop(graph = graph, cl_list= cl_list, R = Lambda, lay = lay, 
                       spop_color = spop_color, name = paste0("louvain_spop_", lambda_index, ".svg"))
  # Heatmap
  overlap_hm(path, info, cl_list = cl_list, 
             R= Lambda, spop_color, name = paste0("heatmap_",lambda_index, ".png"))
  
  # Heatmap by com
  overlap_by_comm(path, info, graph, res_mat, res_indv , step= lambda_index, paste0("heatmap_white_comm_",lambda_index, ".png"))
  
  # Heatmap membership
  hmp <- read.table(file=paste(path,paste("Common_membership_file_R",Lambda, sep=""),sep="/"), header = TRUE )
  heatmap_info <-read.csv(file.path(path, info), sep = "\t", head = TRUE)
  heatmap_info <- heatmap_info[heatmap_info[, 1] %in% colnames(hmp), ]
  
  d <- data.frame(heatmap_info[,5], row.names=heatmap_info[,1])
  colnames(d) <- "Groups"
  Spop_colores <- spop_color[, 2]
  names(Spop_colores) <- spop_color[, 1]
  colors_spop = list(Groups = Spop_colores)
  
  png(file= paste(path,paste("Common_membership_R",lambda_index, ".png", sep=""),sep="/"),width=1000, height=1000)
  ph_com <- pheatmap(as.matrix(hmp),use_raster=TRUE, annotation_row = d, show_rownames = FALSE, 
           show_colnames = FALSE, name = "Common membership", show_row_den=FALSE,
           annotation_colors = colors_spop, fontsize = 6, annotation_col = d)
  draw(ph_com)
  dev.off()
}

im_file = file.path(path, "output/im_50.txt")
or_file = file.path(path, "output/order_50.txt")
ma_file = file.path(path, "output/mapbig_50.txt")
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
  res_colors,
  ma_file,
  sep = ",",
  quote = F,
  row.names = F,
  col.names = F
)

### Input: Interactive Map:
# Leer los archivos de ser necesario: 
imfile = file.path(path, "output/im_50.txt")
im <- read.csv2(imfile, sep = ",", header = F)
order_file = file.path(path, "output/order_50.txt")
ord <- read.csv2(order_file, sep = ",", header = F)
cl <- cl_list[[as.character(Lambda)]]  ########## Afecta el cambio?
ma_file = file.path(path, "output/mapbig_50.txt")
mapbig <- read.csv2(ma_file, sep = ",", header = F)

# Llamando a la funcion: shiny info
info_map(path= path, cl = cl, im = im, ord = ord, steps = steps, 
         outputname = "shinny_info.txt", info=info, lower_limit=lower_limit, upper_limit=upper_limit)

###################################################################################################################
############################################# OBTENER COLORES SIMILARES ###########################################

#Ya con el grafo obtenemos lo layouts que vamos a usar
set.seed(1)
lay_2D <- layout_with_fr(graph, dim = 2)
set.seed(1)
lay_3D <- layout_with_fr(graph, dim = 3)

# Buscamos primero en qué resolución se encuentran la mayor cantidad de comunidades
# Encontrar número de comunidades por resolución
num_comm <- apply(im, 1, contar_unicos)

# Encontrar las resoluciones donde hay más comunidades
posiciones <- which(num_comm == max(num_comm))

# Seleccionar la resolución más pequeña para tener más individuos por comunidad
posicion_final <- min(posiciones)
step = posicion_final

# Obtenemos el grafo de individuos por comm
# Aquí no importa si el lay es 2D porque sólo queremos el grafo
# Para la lambda elegida:
graph_by_comm <- plot_louvain_by_comm(ord, im, step, graph, mapbig, lay_2D,
                                      name = "plot_louvain_final_step.svg", path = path)

# Obtenemos la posición media de los individuos por comunidad en 3D
means_table_3D <- media_graph_3D(graph_by_comm, lay_3D)

# Ordenamos las comunidades de menor a mayor en means_table_3D
orden <- order(as.numeric(means_table_3D[,1]))
means_table_3D_ord <- means_table_3D[orden,]

#Verificamos si faltan comunidades
comm_faltantes <- max(im) - max(num_comm)
if(comm_faltantes != 0) {
  B <- logical(max(im))
  B[as.numeric(means_table_3D_ord[, 1])] <- TRUE
  #all_comm <- means_table_3D_ord
  for (i in 1:max(im)) {
    if (B[i] != TRUE) {
      # Para la comunidad que nos falta
      posiciones_i <- which(im == i, arr.ind = TRUE)
      # Seleccionar la posición más pequeña
      posicion_final_i <- min(posiciones_i[, 1])
      #graph_by_comm para la comunidad i
      graph_by_comm_i <-
        plot_louvain_by_comm(ord, im, posicion_final_i, graph, mapbig, lay_2D,
                             name = "plot_louvain_by_com_comunidades_faltantes.svg", path = path)
      #Medias para la comunidad i
      means_table_3D_i <- media_graph_3D(graph_by_comm_i, lay_3D)
      medias_i <-
        which(means_table_3D_i[, 1] == i, arr.ind = TRUE)
      means_table_3D_ord <-
        rbind(means_table_3D_ord, means_table_3D_i[medias_i, ])
    }
  }
} else{
  cat("Aleluya! Están todas las comunidades")
}

#Organizamos a las comunidades de menor a mayor
orden_2 <- order(as.numeric(means_table_3D_ord[,1]))
means_table_3D_ord <- means_table_3D_ord[orden_2,]

lay_3D_colors <- data2cielab(as.data.frame(means_table_3D_ord[,2:4]))

### Gráfica de resolución
lay_3D_colors[1,2] = "#FFFFFF"
im2 <- apply(im, 2, rev)
# Fijar los valores de los ejes
x_vals <- 1:ncol(im)
y_vals <- 1:nrow(im)

# Crear la imagen con valores de ejes fijos
png(file= paste(path,"resolution_plot_similarC.png",sep="/"))
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
  las = 1
)
mtext(
  side = 2,
  line = 2,
  text = "Resolution value",
  cex = 1,
)
yticklabels <- rev(lower_limit:upper_limit)
# Crear el vector de posiciones para los ticks del eje y
yticks <- seq(1, nrow(im), length.out = length(yticklabels))
# Configurar los ticks y etiquetas del eje y
axis(
  side = 2,
  at = yticks,
  labels = yticklabels,
  las = 1
)
dev.off()

# Colors for shiny app
write.table(
  lay_3D_colors$V2,
  file.path(path, "output/similar_shinyCol.txt"),
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
colors_cielab <- matrix(0, nrow = length(lay_3D_colors[,2]), ncol = 1)
colors_cielab[,1] <- lay_3D_colors[, 2]

## Listas
graph_distintive_list <- list()
graph_similar_list <- list()

# x 1:50 matching total of lambda values 
for (x in 1:steps){
  ######################################################################## Colores distintivos 
  graph_by_comm <- plot_louvain_by_comm(ord, im, step= x, graph, mapbig, lay, 
                                        name = paste0("plot_louvain_comm_L", x, ".svg"), path = path)

  # Red de nodos coloreados por comunidad:
  name_net <- paste0("Network_ColorNodesDistintive", x, ".svg")
  svg(file= paste(path,name_net,sep="/"), width = 30, height = 20)
  plot(graph_by_comm, vertex.size = 3, layout=lay, vertex.label=NA, mark.col= V(graph_by_comm)$color)
  dev.off()
  
  means_table <- media_graph(graph_by_comm, lay)
  freq_ordenado <- num_vert(graph_by_comm)
  tabla_freq_nonzeros <- dens_con(graph_by_comm)
  graph_only_comms <- graph_comms(graph_by_comm, tabla_freq_nonzeros, freq_ordenado)
  
  ## Plot 
  plot_by_density(graph_only_comms, means_table, 
                  path, plotname= paste0("plot_comm_L", x , ".png"))
  
  # Quitamos la comunidad 1
  if(any(V(graph_only_comms)$carac == 1)) {
    graph_only_comms_no1 <- delete_vertices(graph_only_comms, 1)
  } else{
    graph_only_comms_no1 <- graph_only_comms
  }
  
  plot_by_density(graph_only_comms_no1, means_table, path, 
                  plotname= paste0("plot_comm_DistintiveC_L", x , ".png"))
  
  # RED 3D
  graph_3d = graph_only_comms_no1
  # Normalizamos vértices y edges
  V(graph_3d)$size <- (rank(V(graph_3d)$size) / length(V(graph_3d)$size)) * 50
  E(graph_3d)$width <- (rank(E(graph_3d)$width) / length(E(graph_3d)$width)) * 20
  # Layout en 3D seleccionado para las comunidades en el grafo
  coords <- means_table_3D_ord[as.numeric(V(graph_only_comms_no1)$name),]
  
  # Añadimos las coordenadas a un solo objeto
  graph_3d$coords <- coords[,-1]
  
  graph_distintive_list[[x]] <- graph_3d
  
  # Delete objects 
  rm(graph_by_comm, means_table, freq_ordenado, tabla_freq_nonzeros, graph_only_comms, graph_only_comms_no1,
     graph_3d, coords, name_net)
  
  ############################################################################ Plots with similar colors 
  graph_by_comm <- plot_louvain_by_comm(ord, im, step= x, graph, colors_cielab, lay_2D,
                                        name = paste0("plot_louvain_similarC_L", x, ".svg"), path = path)
  
  # Red de nodos coloreados por comunidad:
  name_net <- paste0("Network_ColorNodesSimilar", x, ".svg")
  svg(file= paste(path,name_net,sep="/"), width = 30, height = 20)
  plot(graph_by_comm, vertex.size = 3, layout=lay, vertex.label=NA, mark.col= V(graph_by_comm)$color)
  dev.off()
 
   means_table <- media_graph(graph_by_comm, lay_2D)
  freq_ordenado <- num_vert(graph_by_comm)
  tabla_freq_nonzeros <- dens_con(graph_by_comm)
  graph_only_comms <- graph_comms(graph_by_comm, tabla_freq_nonzeros, freq_ordenado)
  
  # Quitamos la comunidad 1
  if(any(V(graph_only_comms)$carac == 1)) {
    graph_only_comms_no1 <- delete_vertices(graph_only_comms, 1)
  } else{
    graph_only_comms_no1 <- graph_only_comms
  }
  
  plot_by_density(graph_only_comms_no1, means_table, path, 
                  plotname= paste0("plot_comm_similarC_L", x , ".png"))
  
  # RED 3D
  graph_3d = graph_only_comms_no1
  # Normalizamos vértices y edges
  V(graph_3d)$size <- (rank(V(graph_3d)$size) / length(V(graph_3d)$size)) * 50
  E(graph_3d)$width <- (rank(E(graph_3d)$width) / length(E(graph_3d)$width)) * 20
  # Layout en 3D seleccionado para las comunidades en el grafo
  coords <- means_table_3D_ord[as.numeric(V(graph_only_comms_no1)$name),]
  
  # Añadimos las coordenadas a un solo objeto
  graph_3d$coords <- coords[,-1]
  
  graph_similar_list[[x]] <- graph_3d
  
  # Delete objects 
  rm(graph_by_comm, means_table, freq_ordenado, tabla_freq_nonzeros, graph_only_comms, graph_only_comms_no1, 
     graph_3d, coords, name_net)
}

saveRDS(graph_distintive_list, file.path(path,"3Dplots_distintive.rds"))
saveRDS(graph_similar_list, file.path(path,"3Dplots_similar.rds"))

