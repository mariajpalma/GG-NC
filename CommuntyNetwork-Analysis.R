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

ReadGRMBin = function(prefix, AllN = F, size = 4) {
  sum_i = function(i) {
    return(sum(1:i))
  }
  BinFileName = paste(prefix, ".grm.bin", sep = "")
  NFileName = paste(prefix, ".grm.N.bin", sep = "")
  IDFileName = paste(prefix, ".grm.id", sep = "")
  id = read.table(IDFileName)
  n = dim(id)[1]
  BinFile = file(BinFileName, "rb")
  
  grm = readBin(BinFile,
                n = n * (n + 1) / 2,
                what = numeric(0),
                size = size)
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
  i = sapply(1:n, sum_i)
  return(list(
    diag = grm[i],
    off = grm[-i],
    id = id,
    N = N
  ))
}

grm_network <- function (path, fname1, fname2, max) {
  datafile = file.path(path, fname1)
  GRM = ReadGRMBin(datafile)
  # Definir la diagonal y el triángulo superior de la matriz
  diagonal <- GRM[[1]]  # Reemplaza con tus valores
  triangulo_superior <- GRM[[2]]  # Reemplaza con tus valores
  id <- GRM[[3]]
  id <- id[, 2]
  
  # Calcular el número de filas y columnas de la matriz
  n <- length(diagonal)
  matriz_simetrica <- matrix(0, n, n)
  
  # Llenar la matriz con los valores de la diagonal
  matriz_simetrica[upper.tri(matriz_simetrica)] <-
    triangulo_superior
  matriz_simetrica <- matriz_simetrica + t(matriz_simetrica)
  rownames(matriz_simetrica) <- GRM$id[, 2]
  colnames(matriz_simetrica) <- GRM$id[, 2]
  net_grm <-
    graph.adjacency(
      adjmatrix = matriz_simetrica,
      weighted = T,
      diag = F,
      mode = "lower"
    )
  net_grm <- delete.edges(net_grm, which(E(net_grm)$weight < max))
  
  #Asignar colores y spop
  infofile <- file.path(path, fname2)
  info <- read.csv(file.path(path, fname2), sep = "\t", head = TRUE)
  indv_spop_color <- info[, c(1, 5, 8)]
  #vert <- as.data.frame(net_grm, what = "vertices")
  vert <- names(V(net_grm))
  
  filtrado <- indv_spop_color[indv_spop_color[, 1] %in% vert,]
  id_order <- match(filtrado[, 1], vert)
  id_ordenado <- filtrado[order(id_order), ]
  
  info_nodup <-
    info[!duplicated(info$genetic_region, fromLast = FALSE), ]
  spop_color <- cbind(info_nodup$genetic_region, info_nodup$color)
  colnames(spop_color) <- c("genetic_region", "color")
  
  #Agregar en vertices los nombres de las superpoblaciones y sus colores como atributos.
  #Ponerle los colores en info.
  V(net_grm)$spop <- id_ordenado[, 2]
  V(net_grm)$color <- id_ordenado[, 3]
  
  return(list(net_grm, spop_color))
}

input <- function(kind, path, data, info, max, prune) {
  # Verificar si la carpeta ya existe
  directorio = file.path(path,"output")
  if (!file.exists(directorio)) {
    # Si no existe, crear la carpeta
    dir.create(directorio)
    cat("Carpeta 'output' creada correctamente.\n")
  } else {
    cat("La carpeta 'output' ya existe.\n")
  }
  
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
  return(network)
}

get_lambda_results <- function(input, steps, path, name) {
  list_cl <- list()
  ncom <- 0
  a <- 1
  for (i in logspace(-2, 2, n = steps)) {
    cat(a, "\n")
    cl <- cluster_louvain(input, resolution = i)
    list_cl[[as.character(i)]] <- cl
    
    write.table(
      t(cl$membership),
      sep = ",",
      file = file.path(path, name),
      append = TRUE,
      quote = F,
      row.names = F,
      col.names = F
    )
    ncom[a] <- max(cl$membership)
    a <- a + 1
  }
  return(list_cl)
}

nodebased <- function(Sloop) {
  #This R function, `nodebased(Sloop)`, is utilized to re-label communities while ensuring their continuity.
  #This re-labeling process allows us to ascertain how individuals are adapting and readapting within various communities across lambda.
  #The function requires a matrix, `Sloop`, which includes community assignments for each node at varying parameter values.
  #It is assumed that these parameter values represent the highest resolution.
  N <- nrow(Sloop)
  Ssort <- matrix(0, nrow = N, ncol = ncol(Sloop))
  S1 <- Sloop[N,]
  Ssort[N,] <- S1
  #The loop provided iterates from `i = 1` to `N - 1` in order to reassign labels to Ssort, starting from the N-1th row (N-i) and moving towards the first row in reverse order.
  #The loop will continue until all parameter values (resolutions) have been processed. As a result, `Ssort` will hold the re-labeled communities for each parameter value.
  for (i in 1:(N - 1)) {
    cat(i, "\n")
    S2 <- Sloop[N - i,]
    currentmax <- max(Ssort)
    Snew <- relabel(i, S1, S2, currentmax)
    Ssort[N - i,] <- Snew
    S1 <- Snew
  }
  return(Ssort) #8940.427 6892.974 18566.67 1.035 0.2 GRM
}

relabel <- function(num, S1, S2, currentmax){
  #Function for pairwise re-labeling
  #The purpose of the `relabel` function is to re-label communities in `S2` while preserving their identity.
  #This function takes several arguments, including `num`, `S1`, `S2`, and `currentmax`.
  #The communities in `S2` are relabeled based on the mapping provided by `S1`.
  Snew <- matrix(0, nrow = 1, ncol = length(S1))
  nocoms2 <- length(unique(S2))
  coms2 <- sort(unique(S2))
  nocoms1 <- length(unique(S1))
  coms1 <- sort(unique(S1))
  n <- 0
  C <- matrix(0, nrow = nocoms1, ncol = nocoms2)
  l <- matrix(0, nrow = nocoms1 * nocoms2, ncol = 3)
  A <- matrix(0, nrow = length(S1), ncol = nocoms1)
  for (i in 1:nocoms1)
    A[, i] <- (as.numeric(S1 == as.numeric(coms1[i])))
  B <- matrix(0, nrow = length(S2), ncol = nocoms2)
  for (i in 1:nocoms2)
    B[, i] <- (as.numeric(S2 == as.numeric(coms2[i])))
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
  val <- l[, 1]
  ind <- order(val, decreasing = TRUE)
  l2 <-
    l[ind,] # this has highest overlap in top row, etc. First column is value of overlap, second is label you want, third is what it was called in the partition you want to relabel
  m <- 0
  assigned <-
    matrix(0, nrow = 1, ncol = 2) # just added to make coding easier, ignored later
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

  noverlap <- 0
  noverlap_coms <- c()
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
        # more coms. Unusual if looking at decreasing res parameter, but does happen.
        # In this case give unlabeled coms brand new labels
        currentmax= currentmax + noverlap
        unassigned = sort(setdiff(coms2, assigned[, 2])) # communities without labels
        for (k in 1:length(unassigned)) {
          assigned = rbind(assigned, c(currentmax + k, unassigned[k]))
        }
      } else {
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
          # first column assigned gives community label from old partition, second column gives label of new one.
          ## Esta documentación está al revés
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
}

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
  N <- nrow(Ssort)
  order <- 1:ncol(Ssort)
  for (l in 1:(N - r)) {
    # r will be the perfectly ordered partition
    ix <- order(order = Ssort[N + 1 - l,])
    Ssort <- Ssort[, ix]
    order <- order[ix]
  }
  
  # Now arrange to colour all communities that are never over size 'little' to
  # be painted a uniform colour, here, white
  bigcoms <- integer(0)
  m <- max(Ssort)
  sizes <- matrix(0, nrow = N, ncol = m)
  for (i in 1:N) {
    b <- sort(unique(Ssort[i,]))
    for (j in 1:length(b)) {
      sizes[i, b[j]] <- sum(Ssort[i,] == b[j])
    }
  }
  for (k in 1:m) {
    sizescom <- sizes[, k]
    if (sum(sizescom > little) > 1) {
      bigcoms <- c(bigcoms, k)
    }
  }
  smallcoms <- sort(setdiff(1:m, bigcoms))
  Ssort2 <- Ssort
  for (i in 1:length(smallcoms)) {
    y <- which(Ssort == smallcoms[i])
    Ssort2[y] <- 0
  }
  
  # rename coms, to take into account have lost all littlecoms
  rep <-
    sort(unique(as.vector(Ssort2))) #Súper impoartante el as.vector
  a <- 1:length(rep)
  im <- Ssort2
  for (i in 1:length(rep)) {
    z <- which(Ssort2 == rep[i])
    im[z] <- a[i]
  }
  
  m <- max(im)
  cat("num comm: ", max(im), "\n")
  colores <- distinct_colors(
    (max(im) - 1),
    minimal_saturation = 33,
    minimal_lightness = 20,
    maximal_lightness = 80
  )
  map <- colores$name #Colores en hex
  map <- map[sample(length(map))]
  ma <- vector("character", length = m)
  ma[2:m] <- map
  ma[1] <-
    c("#FFFFFF") # makes first one white. Little coms all have label 0, so will appear white
  im2 <- apply(im, 2, rev)
  
  # Fijar los valores de los ejes
  x_vals <- 1:ncol(Ssort)
  #y_vals <- seq(-2, 2, length.out = 50)
  y_vals <- 1:nrow(Ssort)
  # Crear la imagen con valores de ejes fijos
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
  yticklabels <- rev(-2:2)
  # Crear el vector de posiciones para los ticks del eje y
  yticks <- seq(1, nrow(Ssort), length.out = 5)
  # Configurar los ticks y etiquetas del eje y
  axis(
    side = 2,
    at = yticks,
    labels = yticklabels,
    las = 1
  )
  dev.off()
  
  #order como lista de los individuos
  #cl_name <- read.csv(file.path(path, cl_name_file), sep = ",", head = FALSE)
  #cl_name_ord <- cl_name[order]
  
  return(list(im, order, ma))
}

plot_louvain_by_spop <- function(graph, cl_list, spop_color, R, lay, name) {
  cl <- cl_list[[as.character(R)]]
  svg(file= paste(path,name,sep="/"), height = 20, width = 30)
  plot(
    cl,
    graph,
    layout = lay,
    vertex.size = 6,
    vertex.label = NA,
    edge.color = "gray",
    col = V(graph)$color
  )
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

# Stability plots 
stability_matrix_metrics <- function (network, R, metric)
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
  for (i in 1:100){
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
  return(mpair)
}

df_metrics <- function(graph, metric_name, steps) {
  DF_lambda <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(DF_lambda) <- c(toupper(metric_name), "lambda")
  
  # Define the number of cores to use for parallel processing
  num_cores <- detectCores()
  cluster <- makeCluster(num_cores)
  registerDoParallel(cluster)
  
  # Define the range of lambda values
  lambda_values <- logspace(-2, 2, n = steps)
  
  # Define the function to be executed in parallel
  calculate_stability <- function(lambda) {
    smn <- stability_matrix_metrics(graph, R = lambda,toupper(metric_name))
    tempo <- as.data.frame(smn[upper.tri(smn, diag = FALSE)])
    tempo$lambda <- lambda
    colnames(tempo) <- c(toupper(metric_name), "lambda")
    return(tempo)
  }
  
  # Perform parallel computation using foreach
  DF_lambda <- foreach(lambda = lambda_values, .combine = rbind, .export = "stability_matrix_metrics", .packages = c("igraph", "aricode")) %dopar% {
    calculate_stability(lambda)
  }
  
  stopCluster(cluster)
  
  return(DF_lambda)
}

plots_metric <- function(DF_lambda, metric_name, path, name){
  ### boxplot
  colnames(DF_lambda) <- c("metric", "lambda")
  box <- ggplot(DF_lambda, aes( x= lambda, y=metric, group=lambda)) + 
    geom_boxplot() + scale_x_continuous(trans='log10') + ggtitle(metric_name)
  ggsave(filename = paste(path,paste0(name,"_boxplot.png"),sep="/"), plot = box, width = 13, height = 7, units = "in",
         bg = "white",
         dpi = 300)
  ####  Mean
  DF_lambda_m <- DF_lambda %>% group_by(lambda) %>% 
    summarise(mean=mean(metric))
  # Convert tibble to df
  DF_lambda_m  <- DF_lambda_m  %>% as.data.frame()
  mean_plot <- ggplot(DF_lambda_m, aes( x= lambda, y=mean)) + 
    geom_point() + scale_x_continuous(trans='log10') + ggtitle(paste0(metric_name,"_mean"))
  ggsave(filename = paste(path,paste0(name,"_mean.png"),sep="/"), plot = mean_plot, width = 13, height = 7, units = "in",
         bg = "white",
         dpi = 300)
  #### Median
  DF_lambda_m <- DF_lambda %>% group_by(lambda) %>% 
    summarise(median=median(metric))
  # Convert tibble to df
  DF_lambda_m  <- DF_lambda_m  %>% as.data.frame()
  median_plot <- ggplot(DF_lambda_m, aes( x= lambda, y=median)) + 
    geom_point() + scale_x_continuous(trans='log10') + ggtitle(paste0(metric_name,"_median"))
  ggsave(filename = paste(path,paste0(name,"_median.png"),sep="/"), plot = median_plot, width = 13, height = 7, units = "in",
         bg = "white",
         dpi = 300)
  #### Variance
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
#outputname es el nombre del archivo final que entrará a la carpeta output/

info_map <- function(path, im, cl, ord, info, steps, outputname) {
  #Cargar el metadata de los individuos
  meta_file = file.path(path, info)
  metadata = read.csv(meta_file, sep = "\t", header = T)
  #Dejar sólo a los individuos y su pop que estén en cl$names
  #Recuerda que a la función nodebases entran en orden del cl$names
  #y que en Pollock se les asigna un lugar diferente para mejorar la
  #visualización y ese orden viene en el vector order que sale de Pollock
  met_indv_pop <- metadata[,c(1, 4,3)]
  filtrado <-
    met_indv_pop[met_indv_pop$indv %in% as.character(cl$names), ]
  #Ordenar como en pollock usando cl$names y ord
  cl_names_ord <- cl$names[ord[, 1]]
  id_order <- match(filtrado[, 1], as.character(cl_names_ord))
  id_ordenado <- filtrado[order(id_order),]
  #Crear la tabla para las comunidades por individuo con su info
  Lambda <-
    rep(logspace(-2, 2, n = steps), times = nrow(id_ordenado))
  #Unimos las poblaciones de cada individuo ya ordenado con lambda
  pp <- rep(id_ordenado$population, each = steps)
  pp_lmbd <- cbind(pp, Lambda)
  #Ya que las pp están ordenadas por indv en im, unimos im con pp_lmbd
  #Así cada individuo tendrá su comunidad según una lambda
  #Primero obtenemos las comunidades en un vector
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
  #Unimos la población de cada individuo que ya viene con su lambda con la
  #comunidad que haya obtenido en esa lambda
  pp_lmbd_comms <- cbind(pp_lmbd[, c(1, 2)], comms)
  #Obtenemos proporciones de comundiades por lambda y pop
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
  #Al crearse la tabla de proporciones, se desordenan las poblaciones y las lambdas
  #Y las comunidades no vienen ordenadas en forma ascendente
  #Ordenamos comunidades en orden ascendete pasando cols a rows
  tprop_by_row <- t(prop_by_row)
  numCom <- as.numeric(rownames(tprop_by_row))
  comms_order <- order(numCom)
  tprop_by_row_ord <- tprop_by_row[comms_order, ]
  #Regresamos rows a cols y procedemos a sacar los nombres de las pops y los lambdas
  prop_by_row_ord <- t(tprop_by_row_ord)
  prop_pp_lmbd <- rownames(prop_by_row_ord)
  split_pp_lmd <- strsplit(prop_pp_lmbd, "_")
  prop_pp_lmbd_df <- as.data.frame(do.call(rbind, split_pp_lmd))
  colnames(prop_pp_lmbd_df) <- c("population", "Lambda")
  #Juntamos el resto de la información de las poblaciones del metadata
  #Primero, de metadata quitamos duplicados
  metadata_nodup <-
    metadata[!duplicated(metadata$population, fromLast = FALSE),]
  
  #metadata_nodup <- metadata_nodup[, c(-1,-7)]
  pp_merge <-
    merge(prop_pp_lmbd_df,
          metadata_nodup,
          by = "population",
          all.x = TRUE)
  #ordenamos las columnas
  pp_merge2 <-
    cbind(pp_merge$population,
          pp_merge$pop3code,
          pp_merge$genetic_region,
          pp_merge$poject,
          pp_merge$latitude,
          pp_merge$longitude,
          pp_merge$Lambda)
  #Renombramos
  
  colnames(pp_merge2) <-
    c("Pop",
      "Pop3code",
      "Genetic_region",
      "Project",
      "Latitud",
      "Longitud",
      "Lambda")
  #Por fin, unimos la info de las poblaciones con sus comunidades
  prop <- prop_by_row_ord
  colnames(prop) <- paste0("C", 1:max(im))
  #Convertimos los valores 0 a NA
  prop[prop == 0] <- "NA"
  final_table <- cbind(pp_merge2, prop)
  #Ordenamos por lambda y luego por pop
  indice_orden_2 <-
    order(final_table[, 1], as.numeric(final_table[, 6]))
  final_table_order <- final_table[indice_orden_2,]
  #Guardamos
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

# Heatmap: por comunidades
overlap_by_comm <- function(path, info, graph, im, ord, step, name){
  #Tomamos los individuos de cl$name
  cl <- cluster_louvain(graph, r = 1) #No importa el valor de r
  cl_name <- cl$name
  #Ordenamos los individuos como los ordenó Pollock
  cl_name_ord <- cl_name[ord]
  #Metadata para tener las poblaciones
  info <-
    read.csv(file.path(path, info), sep = "\t", head = TRUE)
  info = info[, c(1, 3, 5, 8)]
  #Dejamos sólo individuos en cl_name_ord y en ese orden
  filtrado <- info[info[, 1] %in% cl_name_ord, ]
  id_order <- match(filtrado[, 1], cl_name_ord)
  id_ordenado <- filtrado[order(id_order),]
  #Obtenemos las comunidades de im
  imcols = ncol(im)
  comm <- matrix(data = 0,
                 ncol = 1,
                 nrow = ncol(im))
  count = 1
  for (j in 1:imcols) {
    comm[count,1] = im[step, j]
    count = count + 1
  }
  #Unimos id con comunidades
  comm <- paste0("C",comm)
  idv_comm_pop_spop <-
    cbind(id_ordenado,comm)
  idv_comm_pop_spop <- idv_comm_pop_spop[,c(1,2,3,5,4)]
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
  #Damos nombres dependiendo de la cantidad de comunidades
  num_2level <- length(unique(idv_comm_pop_spop[, 3]))
  list_pop_by_spop <-
    aggregate(
      idv_comm_pop_spop[, 1] ~ idv_comm_pop_spop[, 3] + idv_comm_pop_spop[, 2],
      data = idv_comm_pop_spop,
      FUN = length
    )
  annotation_row_spop = data.frame(Spop = factor(list_pop_by_spop[, 1]))
  rownames(annotation_row_spop) = list_pop_by_spop[, 2]
  spop_color <- idv_comm_pop_spop[!duplicated(idv_comm_pop_spop[,5], fromLast = FALSE),]
  spop_color <- spop_color[,c(3,5)]
  Spop_colores = spop_color[, 2]
  names(Spop_colores) = spop_color[, 1]
  ann_colors_spop = list(Spop = Spop_colores)
  mat_overlap <- as.data.frame.matrix(tabla_overlap)
  #https://jokergoo.github.io/2020/05/06/translate-from-pheatmap-to-complexheatmap/
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

#Funciones ara obtener las redes
#Corrige que haya dos edges entre dos comunidades
#1)plot_louvain_by_comm()
#Función para graficar coloreando comunidades
#Función necesaria para crear la tabla para el mapa de Shiny
plot_louvain_by_comm <- function(ord, im, step, graph, mapbig, lay, name, path) {
  #Aquí step se refiere a la hilera de im de donde sacaremos las comunidades
  #Ordenar como en pollock usando cl$names y ord
  cl <- cluster_louvain(graph, resolution = 1) #Cualquier r sirve, sólo necesitamos names
  cl_names <- cl$names
  cl_names_ord <- cl_names[ord[,1]]
  #Juntamos id y su comm
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
  #Obtenemos los colores de esas comundidades
  colors <- mapbig[comms,1]
  comms_colors <- cbind(comms,colors)
  #Ordenamos las comunidades con su color de menor a mayor
  #comm_color_ord <- order(as.numeric(comms_colors[, 1]))
  #comms_colors <- comms_colors[comm_color_ord,]
  colnames(comms_colors) <- c("Comm","Colors")
  #Creamos matriz de individuos con su comundidad y color
  id_colors <- merge(id_comm,comms_colors, by = "Comm")
  id_colors2 <- cbind(id_colors[,2],id_colors[,1],id_colors[,3])
  colnames(id_colors2) <- c("ID","Comm","Colors")
  #Ordenamos la información con los individuos en en grafo
  vert <- V(graph)$name
  id_order <- match(id_colors2[,1], vert)
  id_ordenado <- id_colors2[order(id_order),]
  #Asignamos los atributos de comunidad y color al grafo
  V(graph)$carac <- id_ordenado[,2]
  V(graph)$color <- id_ordenado[,3]
  #Graficamos
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

#OJO!! Aquí necesitamos el graph que se genera en
#plot_louvain_by_comm()
#Obtener posición media de x y y de lay
#vertex_attr_names(graph) ayuda a ver qué atributos hay
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
#Número de nodos por comunidad
num_vert <- function(graph_by_comm) {
  # Obtener la tabla de frecuencias
  tabla_frecuencias <- table(V(graph_by_comm)$carac)
  #También sale con table(as.numeric(im[1,]))
  
  # Convertir la tabla de frecuencias en un data frame
  df_frecuencias <- as.data.frame(tabla_frecuencias)
  df_frecuencias <- as.matrix(df_frecuencias)
  # Renombrar las columnas
  colnames(df_frecuencias) <- c("Comm", "Frecuencia")
  ord <- order(as.numeric(df_frecuencias[,1]))
  freq_ordenado <- df_frecuencias[ord,]
  
  return(freq_ordenado)
}

#Densidad de conexiones
#Atributos: edge_attr_names(graph)
dens_con <- function(graph_by_comm) {
  #Obtenemos pares de individuos del grafo y nombramos
  conections <- get.edgelist(graph_by_comm)
  colnames(conections) <- c("ID1", "ID2")
  #Obtenemos lista de individuos y la comunidad a la que pertenecen
  indv_comm <- cbind(V(graph_by_comm)$name, V(graph_by_comm)$carac)
  colnames(indv_comm) <- c("ID1", "Comm")
  #Obtenemos las comunidades de ID1
  merge_1 <- merge(conections, indv_comm, by = "ID1")
  colnames(merge_1) <- c("ID1", "ID2","Comm1")
  #Cambiamos nombre para obtener comunidades para ID2
  colnames(indv_comm) <- c("ID2", "Comm")
  merge_2 <- merge(merge_1, indv_comm, by = "ID2")
  #Ordenamos y nombramos la matriz
  merged <- merge_2[,c(2,1,3,4)]
  colnames(merged) <- c("ID1", "ID2","Comm1","Comm2")
  #Ahora sólo nos quedamos con Comm1 y Comm2
  comms <- cbind(merged$Comm1,merged$Comm2)
  # Ordenar alfabéticamente los nombres de las entidades en cada fila
  relaciones_ordenadas <- t(apply(comms, 1, function(row) {
    if (row[1] < row[2]) {
      return(c(row[1], row[2]))
    } else {
      return(c(row[2], row[1]))
    }
  }))
  #Convertir la matriz resultante en un data frame
  relaciones_ordenadas <- as.data.frame(relaciones_ordenadas)
  #Obtenemos tabla de frecuencias
  ####AQUÍ ESTÁ EL ERROR!!! No era matriz cuadrada
  #freq2 <- xtabs( ~ relaciones_ordenadas[,1] + relaciones_ordenadas[,2])
  #Esto lo arregla:
  vector1 <- relaciones_ordenadas[,1]
  vector2 <- relaciones_ordenadas[,2]
  combined_levels <- unique(c(vector1, vector2))
  freq_table <- table(factor(vector1, levels = combined_levels), factor(vector2, levels = combined_levels))
  
  # Muestra la tabla de frecuencias
  #Hay conexiones dentro de las comunidades
  #Esas no las necesitamos y las hacemos 0
  for(i in 1:ncol(freq_table)){
    freq_table[i,i] = 0
  }
  #tabla_freq <- as.data.frame(freq)
  tabla_freq2 <- as.data.frame(freq_table)
  #En la tabla de frecuencia se integran todos contra todos
  #Así que quedan muchos ceros. Los quitamos.
  # nonzeros <- which(tabla_freq[,3] != 0)
  # tabla_freq_nonzeros <- tabla_freq[nonzeros,]
  # tabla_freq_nonzeros <- as.matrix(tabla_freq_nonzeros)
  # colnames(tabla_freq_nonzeros) <- c("from","to","weight")
  nonzeros2 <- which(tabla_freq2[,3] != 0)
  tabla_freq_nonzeros2 <- tabla_freq2[nonzeros2,]
  tabla_freq_nonzeros2 <- as.matrix(tabla_freq_nonzeros2)
  colnames(tabla_freq_nonzeros2) <- c("from","to","weight")
  
  return(tabla_freq_nonzeros2)
}

#Hacer el grafo y plotear
graph_comms <- function(graph_by_comm, tabla_freq_nonzeros, freq_ordenado) {
  #Obtenemos los nombres de las comundiades de menor a mayor
  vert <- sort(as.numeric(unique(V(graph_by_comm)$carac)))
  #Obtenemos los colores de las comunidades
  color <- unique(V(graph_by_comm)$color)
  #Obtenemos el orden de los colores para que coincidan
  #con la posición de las comunidades
  orden <- order(as.numeric(unique(V(graph_by_comm)$carac)))
  color <- color[orden]
  #Creamos un grafo con esta información
  graph_only_comms <-
    graph_from_data_frame(tabla_freq_nonzeros,
                          directed = F,
                          vertices = vert)
  #Asignamos atributos para graficar
  V(graph_only_comms)$carac <- vert
  V(graph_only_comms)$size <- as.numeric(freq_ordenado[, 2])
  V(graph_only_comms)$color <- color
  E(graph_only_comms)$weight <- as.numeric(E(graph_only_comms)$weight)
  E(graph_only_comms)$width <- E(graph_only_comms)$weight
  
  return(graph_only_comms)
}

# Función para contar números únicos por fila
contar_unicos <- function(fila) {
  longitud <- length(unique(fila))
  return(longitud)
}
#Medias de layout en 3D
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

# Para graficar el nuevo grafo
# En el archivo Pa_colores, no parece integrar la nueva normalización al grafo
# o a la función que grafica, por lo que no veríamos esos cambios en la gráfica.
# Vnorm = (rank(V(graph_only_comms_no1)$size) / length(V(graph_only_comms_no1)$size)) * 100
# Enorm = (rank(E(graph_only_comms_no1)$width) / length(E(graph_only_comms_no1)$width)) * 100
# Una forma de hacer lo anterior sería otorgando nuevos valores de la siguiente manera:
# V(graph)$size <- Vnorm
# E(graph)$width <- Enorm
# Pero se perderían los datos originales. Podríamos copiar el grafo antes y hacer la
# modificación en la copia.
# Pero me parece más sencillo si sólo modificamos la función "plot_by_density_no_legend"
# para que sea:
# vertex.size = (rank(V(graph_only_comms_no1)$size) / length(V(graph_only_comms_no1)$size)) * 100,
# edge.width = Enorm = (rank(E(graph_only_comms_no1)$width) / length(E(graph_only_comms_no1)$width)) * 100,
# Y así no tener que modificar o copiar el grafo.
# Como dice mi asesor de doctorado: Nos ahorramos tinta.
plot_by_density <- function(graph_only_comms_no1, means_table, path, plotname) {
  # La función ya no necesita introducir a Vnorm y a Enorm como parámetros
  #ordenar means_table
  comm_order <- match(as.numeric(V(graph_only_comms_no1)$carac),as.numeric(means_table[, 1]))
  means_table_ord_filt <- means_table[comm_order,]
  lay_mean <- as.matrix(means_table_ord_filt[,c(2,3)])
  #Graficamos y guardamos
  png(file.path(path, plotname))
  plot(
    graph_only_comms_no1,
    layout = lay_mean,
    #### Aquí inicia el cambio
    vertex.size = (rank(V(graph_only_comms_no1)$size) / length(V(graph_only_comms_no1)$size)) *
      50,
    edge.width = (rank(E(graph_only_comms_no1)$width) / length(E(graph_only_comms_no1)$width)) *
      20,
    #### Aquí termina el cambio
    edge.color = "gray",
    col = V(graph_only_comms_no1)$color,
  )
  dev.off()
}

prune_network <- function (net){
  while (sum(degree(net) == 1)> 0){
    net <- delete_vertices(net, names(which(degree(net) == 1)))
    net <- delete_vertices(net, names(which(degree(net) == 0)))
  }
  
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

#For IBD
#kind = "IBD"
#path = "C:\\Users\\Brenda Elizabeth L\\Documents\\sohail_Lab\\Last_call\\IBD"
#data = "IBD_cM_sum_NoDivision_segmentSize5_NoSelf_relatedrm"#Must have 3 columns
#info = "info_file_180424.txt"
#shiny_info = "info_file_180424.txt"
#max = 1000#To filter
#steps = 50#Number of lambdas to explore
#min_comms = 6#Tamaño mínimo de las comunidades
#r = 0#r is which layer (counting from bottom) will be perfectly sorted
#(choose on what looks prettiest. Try 0 as default.)
#Lambda = 0.01#Resolución para cl
#Prune = 1

#For PCA
#kind = "PCA"
#path = "C:\\Users\\Brenda Elizabeth L\\Documents\\sohail_Lab\\PCA_shiny"
#data = "PCA_ind_values_NoNull_V2.txt"
#info = "info_file.txt"
#max = 0 #No se usa
#steps = 50#Number of lambdas to explore
#min_comms = 20#Tamaño mínimo de las comunidades
#r = 0#r is which layer (counting from bottom) will be perfectly sorted
#(choose on what looks prettiest. Try 0 as default.)
#Lambda = 0.01#Resolución para cl

#For GRM
#kind = "GRM"
#path = "/data3/PGC/blopez/GRM_results"
#data = "hgdp_tgp_qced_autosomesMergedGRM_rare_NSingletons"
#info = "info_file.txt"
#max = 0 #No se usa
#steps = 50#Number of lambdas to explore
#min_comms = 20#Tamaño mínimo de las comunidades
#r = 0#r is which layer (counting from bottom) will be perfectly sorted
#(choose on what looks prettiest. Try 0 as default.)
#Lambda = 0.01 #Resolución para cl

#The permissible options for kind:
#'IBD' (Identical by descent),
#'PCA'(Principal Component Analysis) or
#'GRM' (Genetic Relationship Matrix).

args <- commandArgs(trailingOnly = TRUE)
kind = args[1]
path = args[2]
data = args[3]
info = args[4]
max = as.numeric(args[5])
steps = as.numeric(args[6])
Lambda = as.numeric(args[7]) #Resolución para cl
Prune = as.numeric(args[8])
min_comms = as.numeric(args[9])
r = 0 #r is which layer (counting from bottom) will be perfectly sorted

#shiny_info = args[8]
#min_comms = 6 #Tamaño mínimo de las comunidades
#r = 0 #r is which layer (counting from bottom) will be perfectly sorted

# Creating graph 
network <- input(kind, path, data, info, max, Prune)
graph <- network[[1]]
spop_color <- network[[2]]
lay <- layout_with_fr(graph)
tic()
cl_list <- get_lambda_results(graph, steps, path, "output/M_IBD_50.txt")
toc()

set.seed(2)
lay=layout_with_fr(graph)
svg(paste(path,"Graph_Simple.svg",sep="/"))
plot(graph, vertex.size = 3, layout=lay, vertex.label=NA)
dev.off()


# Plot: Network 
plot_louvain_by_spop(graph = graph, cl_list= cl_list, R = Lambda, lay = lay, 
                     spop_color = spop_color, name = "louvain_spop.svg")

plot_louvain_by_spop(graph = graph, cl_list= cl_list, 
                     R = as.character(logspace(-2,2,50)[5]), lay = lay, 
                     spop_color = spop_color, name = "louvain_spop_5.svg")

plot_louvain_by_spop(graph = graph, cl_list= cl_list, 
                     R = as.character(logspace(-2,2,50)[9]), lay = lay, 
                     spop_color = spop_color, name = "louvain_spop_9.svg")

# Plot: Heatmap 
overlap_hm(path, info, cl_list = cl_list, 
           R= Lambda, spop_color, name = "heatmap.png")

overlap_hm(path, info, cl_list = cl_list, 
           R= as.character(logspace(-2,2,50)[5]), spop_color, name = "heatmap_5.png")

overlap_hm(path, info, cl_list = cl_list, 
           R= as.character(logspace(-2,2,50)[9]), spop_color, name = "heatmap_9.png")

plots_metric(df_metrics(graph, "nid", steps = steps), metric_name = "NID", 
             path = path, name = "TEST2_NID")

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
cl <- cl_list[[as.character(Lambda)]]
ma_file = file.path(path, "output/mapbig_50.txt")
mapbig <- read.csv2(ma_file, sep = ",", header = F)

# Llamando a la funcion: shiny info
info_map(path= path, cl = cl, im = im, ord = ord, steps = steps, 
         outputname = "shinny_info.txt", info=info)

### Plot: Heatmap by comm
overlap_by_comm(path, info, graph, res_mat, res_indv , step= which(logspace(-2,2,50)==Lambda), "heatmap_white_comm.png")
overlap_by_comm(path, info, graph, res_mat, res_indv , step= 5, "heatmap_white_comm_5.png")
overlap_by_comm(path, info, graph, res_mat, res_indv , step= 9, "heatmap_white_comm_9.png")

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
#y_vals <- seq(-2, 2, length.out = 50)
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
yticklabels <- rev(-2:2)
# Crear el vector de posiciones para los ticks del eje y
yticks <- seq(1, nrow(im), length.out = 5)
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
for (x in 1:50){
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

# Guardamos todo en archivos
#file_name_distintive <- paste0(path, "\\", kind,"_","3Dplots_distintive",".rds")
#file_name_similar <- paste0(path, "\\", kind,"_","3Dplots_similar",".rds")

saveRDS(graph_distintive_list, file.path(path,"3Dplots_distintive.rds"))
saveRDS(graph_similar_list, file.path(path,"3Dplots_similar.rds"))

# Archivos
# load("your_file.RData")
# rglplot(graph_list[[3]], layout = as.matrix(graph_list[[3]]$coords), vertex.label = NA)
# Esto nos permite ver la gráfica
# rglwidget()

# Guardamos un .html
#htmlwidgets::saveWidget(rglwidget(width = 520, height = 520), 
#                        file = file.path(path, "plot_3D.html"),
#                        libdir = "libs",
#                        selfcontained = TRUE
#)
# Para borrar datos y que la siguiente gráfica no se empalme con la anterior
#clear3d()
