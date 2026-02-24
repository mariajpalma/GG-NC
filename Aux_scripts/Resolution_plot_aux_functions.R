### Auxiliar script
### This script aims to compare community detection results among different metrics and resolutions
### The x-axis order of the resolution plots can be defined by the super populations of the order of other metric
### We also added a color bar to show the super populations

library(dplyr)
library(png)
library(grid)
library(ggpubr)
library(ggplot2)

################################################################################
#                          res_plot_mash_up_slice                               #
################################################################################
# This function takes two resolutions (or steps) from the res_plots for comparison.
# It orders the individuals in metric B according to the order in metric A.
# If there are individuals in A that are not present in B, they are colored black in B
# at the position corresponding to the individual in A.
# If there are individuals in B that are not present in A, they are colored black in A
# and, according to their color in B, placed at the end of the slices.

resplot_mash_up_slice <- function(metricA,
                                  metricB,
                                  RA,
                                  RB,
                                  OA,
                                  OB,
                                  coloresA,
                                  coloresB,
                                  RA_lower_limit,
                                  RA_upper_limit,
                                  RB_lower_limit,
                                  RB_upper_limit,
                                  RA_steps,
                                  RB_steps,
                                  path,
                                  plotname,
                                  plotnameA,
                                  plotnameB,
                                  select_by,
                                  RA_val,
                                  RB_val) {
  #' @param metricA a character string specifying the name of the metric that set the individuals order.
  #' @param metricB a character string specifying the name of the metric to compare with.
  #' @param RA is the resolution matrix of metric A.
  #' @param RB is the resolution matrix of metric B.
  #' @param OA is the Community_detection/individual_ID_order.txt file with the individuals ID for metric A.
  #' @param OB is the Community_detection/individual_ID_order.txt file with the individuals ID for metric B.
  #' @param coloresA a vector with the colors for the metric A.
  #' @param coloresB a vector with the colors for the metric B.
  #' @param RA_lower_limit a numeric value indicating the lower bound of the resolution range for clustering of the metric A.
  #' @param RA_upper_limit a numeric value indicating the upper bound of the resolution range for clustering of the metric A.
  #' @param RB_lower_limit a numeric value indicating the lower bound of the resolution range for clustering of the metric B.
  #' @param RB_upper_limit a numeric value indicating the upper bound of the resolution range for clustering of the metric B.
  #' @param RA_steps is the row size of the resolucion matrix of metric A.
  #' @param RB_steps is the row size of the resolucion matrix of metric B.
  #' @param path a character string specifying the path to the directory containing the input files.
  #' @param plotname a character string specifying the file name for saving the final plot.
  #' @param plotnameA a character string specifying the file name for saving the metric A plot.
  #' @param plotnameB a character string specifying the file name for saving the metric B plot.
  #' @param select_by a character string specifying if the value to select a row in the resolution matrix is a "step" or a "resolution" value.
  #' @param RA_val a numeric value, "step" or "resolution" value, for the metric A.
  #' @param RB_val a numeric value, "step" or "resolution" value, for the metric B.
  # Alinear a los individuos de GRMc como en PCA
  RAB_columnas_comunes <- intersect(OA[, 1], OB[, 1])  # Columnas en ambos con orden como en OA
  RAB_columnas_faltantes <- setdiff(OA[, 1], OB[, 1])  # Columnas en OA pero no en OB
  RBA_columnas_faltantes <- setdiff(OB[, 1], OA[, 1])  # Columnas en OB pero no en OA
  # La opción drop = FALSE evita que R convierta matrices de una columna en vectores automáticamente.
  colnames(RA) <- OA[, 1]  # Asignar nombres a RB (necesario para el filtrado)
  colnames(RB) <- OB[, 1]  # Asignar nombres a RB (necesario para el filtrado)
  RA_ordenada <- RA
  # RA_ordenada <- RA[, RAB_columnas_comunes, drop = FALSE]  # Seleccionar columnas comunes
  RB_ordenada <- RB[, RAB_columnas_comunes, drop = FALSE]  # Seleccionar columnas comunes
  # Seleccionar un step, gamma o res
  if (select_by == "step") {
    RA_ordenada <- RA_ordenada[RA_val, ]
    RB_ordenada <- RB_ordenada[RB_val, ]
    RB_all <- RB[RB_val,]
  }
  if (select_by == "resolution") {
    RA_res <- seq(from = RA_lower_limit, to = RA_upper_limit, length.out = RA_steps)
    RB_res <- seq(from = RB_lower_limit, to = RB_upper_limit, length.out = RB_steps)
    val_A <- which(round(RA_res, 6) == RA_val)
    val_B <- which(round(RB_res, 6) == RB_val)
    RB_all <- RB[val_B,]
    if (length(val_A) == 1 && length(val_B) == 1) {
      RA_ordenada <- RA_ordenada[val_A, ]
      RB_ordenada <- RB_ordenada[val_B, ]
    } else{
      stop("There are several positions with this value, please try another value")
    }
  }
  # Elegir los colores para cada métrica
  # Obtenemos las comunidades de la res en orden de menor a mayor
  RA_comms <- unique(sort(as.vector(as.matrix(RA_ordenada))))
  RB_comms <- unique(sort(as.vector(as.matrix(RB_all))))
  # Obtener los colores que corresponden a esas comunidades.
  if (length(RAB_columnas_faltantes) > 0) {
    RB_color <- c("#9b9b9b", coloresB[RB_comms, 1])
  } else{
    RB_color <- coloresB[RB_comms, 1]
  }
  if (length(RBA_columnas_faltantes) > 0) {
    RA_color <- c("#9b9b9b", coloresA[RA_comms, 1])
  } else{
    RA_color <- coloresA[RA_comms, 1]
  }
  # Renumerar las comunidades de RA de 1:length(comms) para colorear correctamente
  RA_renamed <- matrix(as.numeric(factor(
    as.matrix(RA_ordenada),
    levels = RA_comms,
    labels = 1:length(RA_comms)
  )), ncol = ncol(RA_ordenada))
  colnames(RA_renamed) <- names(RA_ordenada)
  # Crear matriz para columnas faltantes de A en B
  # Caso 1, cuando todo elemento de B es también elemento de A.
  if (length(RAB_columnas_faltantes) > 0 && length(RBA_columnas_faltantes) == 0) {
    # Renumerar las comunidades de RB de 1:length(comms) para colorear correctamente
    RB_renamed <- matrix(as.numeric(factor(
      as.matrix(RB_ordenada),
      levels = RB_comms,
      labels = 1:length(RB_comms)
    )), ncol = ncol(RB_ordenada))
    colnames(RB_renamed) <- names(RB_ordenada)
    # Matriz de ceros para los faltantes
    matriz_faltantes <- matrix(0,
                               nrow = nrow(RB_renamed),
                               ncol = length(RAB_columnas_faltantes))
    colnames(matriz_faltantes) <- RAB_columnas_faltantes  # Nombrar columnas
    # Combinar con columnas faltantes
    RB_A <- cbind(RB_renamed, matriz_faltantes)
    #RA_B <- cbind(RA_renamed, matriz_faltantes)
    # Reordenar columnas según el orden de RA
    RB_A <- as.matrix(t(RB_A[, OA[, 1]]))
    RA_B <- RA_renamed
  }
  # Crear matriz para columnas faltantes de B en A
  # Caso 2, cuando todo elemento de B es también elemento de A.
  if (length(RAB_columnas_faltantes) == 0 && length(RBA_columnas_faltantes) > 0) {
    matriz_faltantes <- matrix(0,
                               nrow = nrow(RA_renamed),
                               ncol = length(RBA_columnas_faltantes))
    colnames(matriz_faltantes) <- RBA_columnas_faltantes  # Nombrar columnas
    RB_indiv_falt <- RB_all[, RBA_columnas_faltantes, drop = FALSE]
    RB_ordenada_falt <- cbind(RB_ordenada, RB_indiv_falt)
    # Renumerar las comunidades de RB de 1:length(comms) para colorear correctamente
    RB_renamed <- matrix(as.numeric(factor(
      as.matrix(RB_ordenada_falt),
      levels = RB_comms,
      labels = 1:length(RB_comms)
    )), ncol = ncol(RB_ordenada_falt))
    colnames(RB_renamed) <- names(RB_ordenada_falt)
    RA_B <- cbind(RA_renamed, matriz_faltantes)
    RB_A <- RB_renamed
  }
  # Caso 3, cuando A y B tienen elementos en común y no en común.
  if (length(RAB_columnas_faltantes) > 0 && length(RBA_columnas_faltantes) > 0) {
    RBA_matriz_faltantes <- matrix(0,
                                   nrow = nrow(RA_renamed),
                                   ncol = length(RBA_columnas_faltantes))
    colnames(RBA_matriz_faltantes) <- RBA_columnas_faltantes  # Nombrar columnas
    RA_B <- cbind(RA_renamed, RBA_matriz_faltantes)
    RAB_matriz_faltantes <- matrix(0,
                                   nrow = nrow(RA_renamed),
                                   ncol = length(RAB_columnas_faltantes))
    colnames(RAB_matriz_faltantes) <- RAB_columnas_faltantes  # Nombrar columnas
    RB_indiv_falt <- RB_all[, RBA_columnas_faltantes, drop = FALSE]
    RB_ordenada2 <- cbind(RB_ordenada, RB_indiv_falt)
    RB_renamed <- matrix(as.numeric(factor(
      as.matrix(RB_ordenada2),
      levels = RB_comms,
      labels = 1:length(RB_comms)
    )), ncol = ncol(RB_ordenada2))
    colnames(RB_renamed) <- names(RB_ordenada2)
    RB_ordenada3 <- cbind(RB_renamed[,1:ncol(RB_ordenada), drop = FALSE], RAB_matriz_faltantes)
    colnames(RB_ordenada3) <- c(colnames(RB_renamed[,1:ncol(RB_ordenada), drop = FALSE]), colnames(RAB_matriz_faltantes))
    RB_ordenada4 <- RB_ordenada3[, OA[,1], drop = FALSE]
    RB_A <- cbind(RB_ordenada4, RB_renamed[,(ncol(RB_ordenada)+1):ncol(RB_renamed), drop = FALSE])
  }
  # fix axis values
  x_vals <- 1:length(RA_B)
  y_vals <- 1
  # Graficamos y guardamos
  png(
    file = paste(path, plotnameA, sep = "/"),
    res = 300,
    width = 2000,
    height = 500
  )
  par(mar = c(0, 4, 2, 1)) # mar = c(inf, izq, sup, der)
  image(
    x = x_vals,
    y = y_vals,
    z = t(RA_B),
    useRaster = TRUE,
    col = RA_color,
    axes = F,
    xlab = " ",
    ylab = " "
  )
  if (select_by == "step") {
    mtext(
      side = 2,
      line = 1,
      text = paste0(metricA, "\n Step: ", RA_val),
      cex = 1,
    )
  }
  if (select_by == "resolution") {
    mtext(
      side = 2,
      line = 1,
      text = paste0(metricA, "\nRes: ", RA_val),
      cex = 1,
    )
  }
  dev.off()
  
  png(
    file = paste(path, plotnameB, sep = "/"),
    res = 300,
    width = 2000,
    height = 500
  )
  par(mar = c(2, 4, 0, 1)) # mar = c(inf, izq, sup, der)
  image(
    x = x_vals,
    y = y_vals,
    z = t(RB_A),
    useRaster = TRUE,
    col = RB_color,
    xaxs = "i",
    yaxt = "n",
    #axes = F,
    xlab = " ",
    ylab = " "
  )
  if (select_by == "step") {
    mtext(
      side = 2,
      line = 1,
      text = paste0(metricB, "\n Step: ", RB_val),
      cex = 1,
    )
  }
  if (select_by == "resolution") {
    mtext(
      side = 2,
      line = 1,
      text = paste0(metricB, "\nRes: ", RB_val),
      cex = 1,
    )
  }
  dev.off()
  
  plot1 <-
    create_plot_with_title(file.path(path, plotnameA), "")
  plot2 <-
    create_plot_with_title(file.path(path, plotnameB), " ")
  
  combined_plot <- ggarrange(
    plot1,
    plot2,
    labels = c("", ""),
    ncol = 1,
    nrow = 2
  )
  
  ggsave(
    file.path(path, plotname),
    plot = combined_plot,
    limitsize = FALSE,
    width = 12,
    height = 6,
    dpi = 300
  )
  
}

## Función accesoria para unir ambas slices

create_plot_with_title <- function(file, title) {
  img <- readPNG(file)
  ggplot() +
    annotation_custom(rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc")), -Inf, Inf, -Inf, Inf) +
    ggtitle(title) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 12))
}

# Esta función sirve y es para toda la resplot

################################################################################
#                          res_plot_mash_up                                    #
################################################################################
# This function uses the res_A_B matrix to plot it using the colors from B.
# res_A_B is a matrix obtained by reordering the individuals in B according to A.
# It contains the individuals' IDs as column names.
# If there are individuals in B that are not present in A, they are omitted.
# If there are individuals in A that are not present in B, they are colored black.
resplot_mash_up_all <- function(OA, OB, RB, coloresB, lower_limit, upper_limit, path, name){
  #' @param OA is the Community_detection/individual_ID_order.txt file with the individuals ID for metric A.
  #' @param OB is the Community_detection/individual_ID_order.txt file with the individuals ID for metric B.
  #' @param RB is the community detection matrix for metric B
  #' @param coloresB a vector with the colors for the metric B.
  #' @param path a character string specifying the path to the directory containing the input files.
  #' @param name a character string specifying the file name for saving the final plot.
  colnames(RB) <- OB[,1]  # Asignar nombres a RB (necesario para el filtrado)
  columnas_comunes <- intersect(OA[,1], OB[,1])  # Columnas en ambos: A, B, C, E
  columnas_faltantes <- setdiff(OA[,1], OB[,1])  # Columnas en OA pero no en OB: D
  RB_ordenada <- RB[, columnas_comunes, drop = FALSE]  # Seleccionar columnas comunes
  # Obtenemos las comunidades de la res en orden de menor a mayor
  comms <- unique(sort(as.vector(as.matrix(RB_ordenada))))
  # Obtener los colores que corresponden a esas comunidades.
  if(length(columnas_faltantes) > 0) {
    color <- c("#9b9b9b", coloresB[comms, 1])
  }else{ color <- coloresB[comms, 1]}
  # Renumerar las comunidades de 1:length(comms) para colorear correctamente
  RA_B_renamed <- matrix(
    as.numeric(factor(as.matrix(RB_ordenada), levels = comms, labels = 1:length(comms))),
    ncol = ncol(RB_ordenada)
  )
  colnames(RA_B_renamed) <- names(RB_ordenada)
  # Crear matriz para columnas faltantes (D)
  if (length(columnas_faltantes) > 0) {
    matriz_faltantes <- matrix(0,
                               nrow = nrow(RB),
                               ncol = length(columnas_faltantes))
    colnames(matriz_faltantes) <- columnas_faltantes  # Nombrar columnas
    # Combinar RB_ordenada con columnas faltantes
    RA_B <- cbind(RA_B_renamed, matriz_faltantes)
    # Reordenar columnas según el orden de RA
    RA_B <- RA_B[, OA[, 1]]
  }else{RA_B <- RA_B_renamed}
  # Preparamos la matriz para graficarla con image()
  im2 <- as.matrix(apply(RA_B, 2, rev))
  # fix axis values
  x_vals <- 1:ncol(RA_B)
  y_vals <- 1:nrow(RA_B)
  # Graficamos y guardamos
  png(file = paste(path, name, sep = "/"), res=300, width = 2000, height = 2000)
  image(
    x = x_vals,
    y = y_vals,
    z = t(im2),
    useRaster = TRUE,
    col = color,
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
  yticks <- seq(1, nrow(RA_B), length.out = length(yticklabels))
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

################################################################################
#                          resplot_superpop                                   #
################################################################################
# This function orders the x-axis of the whole matrix according to the pop information
# (for example: super populations) and by number ID
resplot_superpop <- function(OA, RB, OA_colors, pop_info, lower_limit, upper_limit, path, name){
  #' @param OA A matrix or data frame containing the reference ordering of individuals.
  #'   The first column (\code{OA[,1]}) must contain individual IDs, which will be
  #'   assigned as column names of \code{RB}.
  #'
  #' @param RB A numeric matrix of community memberships (e.g., resolution steps × individuals).
  #'   Rows typically represent resolution values and columns represent individuals.
  #'   Community labels are expected to be integers.
  #'
  #' @param OA_colors A matrix or data frame of colors indexed by community label.
  #'   Rows correspond to community IDs and the first column contains color values
  #'   (e.g., hex codes). Colors are selected according to the unique community
  #'   labels present in \code{RB}.
  #'
  #' @param pop_info A data frame containing individual-level metadata. It must include
  #'   at least the columns:
  #'   \itemize{
  #'     \item \code{indv}: individual IDs (matching column names of \code{RB}),
  #'     \item \code{genetic_region}: superpopulation or region label,
  #'     \item \code{color}: color associated with each superpopulation.
  #'   }
  #'   The data frame is used to reorder individuals and to draw the superpopulation
  #'   annotation bar.
  #' @param lower_limit Integer. Lower bound of the resolution values used to label
  #'   the y-axis.
  #' @param upper_limit Integer. Upper bound of the resolution values used to label
  #'   the y-axis.
  #' @param path Character. Output directory where the PNG file will be saved.
  #'
  #' @param name Character. File name (or relative path under \code{path}) for the
  #'   generated PNG image.
  colnames(RB) <- OA[,1]
  common_ids <- intersect(pop_info$indv, colnames(RB))
  pop_info <- pop_info %>% arrange(genetic_region, indv) %>% filter(indv %in% common_ids)
  RB_ordered <- RB[, pop_info$indv, drop = FALSE]
  comms <- unique(sort(as.vector(as.matrix(RB_ordered))))
  colors <- OA_colors[comms,1]
  im2 <- as.matrix(apply(RB_ordered, 2, rev))
  # fix axis values
  x_vals <- 1:ncol(RB_ordered)
  y_vals <- 1:nrow(RB_ordered)
  
  # Obtener colores y etiquetas de superpoblaciones en el orden de los individuos
  superpop_colors <- pop_info$color
  superpop_labels <- pop_info$genetic_region
  
  # Encontrar los límites donde cambia la superpoblación
  superpop_changes <- c(1, which(diff(as.numeric(factor(superpop_labels))) != 0) + 1, ncol(RB_ordered) + 1)
  
  # Configurar el gráfico
  png(file = paste(path, name, sep = "/"), res = 300, width = 2000, height = 2200)
  layout(matrix(c(1, 2), nrow = 2, ncol = 1), heights = c(0.9, 0.1))
  par(mar = c(2, 5, 2, 2))
  image(
    x = x_vals,
    y = y_vals,
    z = t(im2),
    useRaster = TRUE,
    col = colors,
    xaxs = "i",
    yaxt = "n",
    xaxt="n",
    xlab = "",
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
  yticks <- seq(1, nrow(RB_ordered), length.out = length(yticklabels))
  # Configurar los ticks y etiquetas del eje y
  axis(
    side = 2,
    at = yticks,
    labels = yticklabels,
    las = 1,
    cex.axis = 1.6
  )
  mtext("Resolution value", side = 2, line = 2, cex = 1.6)
  
  # Eje y con etiquetas invertidas
  yticklabels <- rev(lower_limit:upper_limit)
  yticks <- seq(1, nrow(RB_ordered), length.out = length(yticklabels))
  axis(2, at = yticks, labels = yticklabels, las = 1, cex.axis = 1.6)
  
  # Barra de superpoblaciones (abajo)
  par(mar = c(3, 5, 0, 2))
  plot(1, type="n", xlim = c(0.5, ncol(RB_ordered) + 0.5), ylim = c(0, 1),
       xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n", axes=F, ann=F)
  
  # Dibujar rectángulos para cada superpoblación
  for (i in 1:(length(superpop_changes) - 1)) {
    x1 <- superpop_changes[i]
    x2 <- superpop_changes[i + 1]
    current_superpop <- superpop_labels[superpop_changes[i]]
    current_color <- superpop_colors[superpop_changes[i]]
    rect(x1 - 0.5, 0, x2 - 0.5, 1, col = current_color, border = NA)
  }
  dev.off()
}


################################################################################
#                         relabel_order_colors                                 #
################################################################################
#' Relabel communities by first appearance and size, then plot a relabeled target matrix
#'
#' This function relabels community IDs in a reference membership matrix (\code{matrix})
#' so that labels are assigned in the order communities first appear across rows
#' (resolution steps), breaking ties by the community size at first appearance
#' (largest first). Community \code{1} is kept as label \code{1}; \code{0} is treated
#' as background and ignored for relabeling.
#'
#' The same relabeling-by-appearance procedure is applied independently to \code{target}
#' to obtain \code{relabeled_target}. A color palette is assigned to communities based
#' on the relabeling of \code{matrix}, and a PNG is written showing \code{relabeled_target}
#' colored according to that palette.
relabel_order_colors <- function(matrix, target, color_palette, path, name, upper_limit, lower_limit){
  #' @param matrix A numeric matrix of community memberships (e.g., a res_plot),
  #'   where rows typically correspond to resolution values and columns to individuals.
  #'   Values are expected to be integers with \code{0} meaning "no community/background".
  #' @param target A numeric matrix of community memberships to be relabeled and plotted.
  #'   Must have the same interpretation as \code{matrix}.
  #' @param color_palette Optional color palette to use for communities. If \code{NULL},
  #'   a palette is generated with \pkg{viridis} (if available) or \code{rainbow()} otherwise.
  #'   The palette is expected to provide at least one distinct color per community label
  #'   (excluding \code{0}). (See Details for accepted shapes.)
  #' @param path Character. Output directory where the PNG will be saved.
  #' @param name Character. File name (or relative path under \code{path}) for the PNG.
  #' @param upper_limit Integer. Upper bound of the resolution index/value used to label
  #'   the y-axis.
  #' @param lower_limit Integer. Lower bound of the resolution index/value used to label
  #'   the y-axis.
  matrix <- as.matrix(matrix)
  target <- as.matrix(target)
  nc_matrix <- max(matrix)
  nc_target <- max(target)
  all_communities <- sort(setdiff(unique(c(matrix)), 0))
  
  # Crear paleta de colores si no se proporciona
  if(is.null(color_palette)) {
    n_colors <- length(all_communities)
    if(requireNamespace("viridis", quietly = TRUE)) {
      color_palette <- viridis::viridis(n_colors)
    } else {
      color_palette <- rainbow(n_colors)
    }
  }
  
  # Verificar que hay suficientes colores
  if(length(color_palette[,]) < length(all_communities)) {
    stop("La paleta de colores no tiene suficientes colores para todas las comunidades")
  }
  
  # Asignar colores secuencialmente (comunidad 1 = color_palette[1], etc.)
  color_mapping <- data.frame(
    old_label = all_communities,
    color = color_palette[1:length(all_communities),],
    stringsAsFactors = FALSE
  )
  
  # Inicializar variables de seguimiento
  relabeled_matrix <- matrix
  next_label <- 2  # Comenzamos a reetiquetar desde 2 (1 se mantiene)
  
  # Diccionario para mapeos
  label_dict <- data.frame(
    old_label = integer(),
    new_label = integer(),
    color = character(),
    size_at_first_appearance = integer(),
    first_appearance_row = integer(),
    stringsAsFactors = FALSE
  )
  
  # Procesar comunidad 1 primero
  comm_1_rows <- which(matrix == 1, arr.ind = TRUE)
  if(length(comm_1_rows) > 0) {
    first_row_1 <- min(comm_1_rows[, 1])
    size_at_first_1 <- sum(matrix[first_row_1, ] == 1)
    
    label_dict <- rbind(label_dict, data.frame(
      old_label = 1,
      new_label = 1,
      color = color_mapping$color[color_mapping$old_label == 1],
      size_at_first_appearance = size_at_first_1,
      first_appearance_row = first_row_1
    ))
  }
  
  # Lista de comunidades ya procesadas (empezamos con 1)
  processed_communities <- c(1)
  
  # Procesar filas en orden
  for(i in 1:nrow(matrix)) {
    current_row <- matrix[i, ]
    unique_in_row <- unique(current_row)
    
    # Identificar comunidades no procesadas en esta fila (excluyendo 0 y 1)
    new_communities <- setdiff(unique_in_row, c(0, 1, processed_communities))
    
    if(length(new_communities) > 0) {
      # Calcular tamaños en primera aparición
      sizes <- sapply(new_communities, function(comm) sum(current_row == comm))
      
      # Ordenar comunidades por tamaño (mayor primero)
      new_communities <- new_communities[order(sizes, decreasing = TRUE)]
      
      # Procesar cada nueva comunidad
      for(comm in new_communities) {
        # Asignar nuevo etiqueta y color
        relabeled_matrix[matrix == comm] <- next_label
        
        # Registrar en el diccionario
        label_dict <- rbind(label_dict, data.frame(
          old_label = comm,
          new_label = next_label,
          color = color_mapping$color[color_mapping$old_label == comm],
          size_at_first_appearance = sum(current_row == comm),
          first_appearance_row = i
        ))
        
        # Actualizar variables de seguimiento
        processed_communities <- c(processed_communities, comm)
        next_label <- next_label + 1
      }
    }
  }
  
  # Ordenar diccionario por primera aparición
  label_dict <- label_dict[order(label_dict$new_label), ]
  
  
  all_communities <- sort(setdiff(unique(c(target)), 0))
  
  # Inicializar variables de seguimiento
  relabeled_target <- target
  next_label <- 2  # Comenzamos a reetiquetar desde 2 (1 se mantiene)
  
  # Diccionario para mapeos
  label_dict2 <- data.frame(
    old_label = integer(),
    new_label = integer(),
    size_at_first_appearance = integer(),
    first_appearance_row = integer(),
    stringsAsFactors = FALSE
  )
  
  # Procesar comunidad 1 primero
  comm_1_rows <- which(target == 1, arr.ind = TRUE)
  if(length(comm_1_rows) > 0) {
    first_row_1 <- min(comm_1_rows[, 1])
    size_at_first_1 <- sum(target[first_row_1, ] == 1)
    
    label_dict2 <- rbind(label_dict2, data.frame(
      old_label = 1,
      new_label = 1,
      size_at_first_appearance = size_at_first_1,
      first_appearance_row = first_row_1
    ))
  }
  
  # Lista de comunidades ya procesadas (empezamos con 1)
  processed_communities <- c(1)
  
  # Procesar filas en orden
  for(i in 1:nrow(target)) {
    current_row <- target[i, ]
    unique_in_row <- unique(current_row)
    
    # Identificar comunidades no procesadas en esta fila (excluyendo 0 y 1)
    new_communities <- setdiff(unique_in_row, c(0, 1, processed_communities))
    
    if(length(new_communities) > 0) {
      # Calcular tamaños en primera aparición
      sizes <- sapply(new_communities, function(comm) sum(current_row == comm))
      
      # Ordenar comunidades por tamaño (mayor primero)
      new_communities <- new_communities[order(sizes, decreasing = TRUE)]
      
      # Procesar cada nueva comunidad
      for(comm in new_communities) {
        # Asignar nuevo etiqueta y color
        relabeled_target[target == comm] <- next_label
        
        # Registrar en el diccionario
        label_dict2 <- rbind(label_dict2, data.frame(
          old_label = comm,
          new_label = next_label,
          size_at_first_appearance = sum(current_row == comm),
          first_appearance_row = i
        ))
        
        # Actualizar variables de seguimiento
        processed_communities <- c(processed_communities, comm)
        next_label <- next_label + 1
      }
    }
  }
  
  # Ordenar diccionario por primera aparición
  label_dict2 <- label_dict2[order(label_dict2$new_label), ]
  
  colors <- rep("gray", nc_target)
  colors[1:min(nc_matrix, nc_target)] <- label_dict$color
  im2 <- as.matrix(apply(relabeled_target, 2, rev))
  # fix axis values
  x_vals <- 1:ncol(relabeled_target)
  y_vals <- 1:nrow(relabeled_target)
  png(file = paste(path, name, sep = "/"), res=300, width = 2000, height = 2000)
  image(
    x = x_vals,
    y = y_vals,
    z = t(im2),
    useRaster = TRUE,
    col = colors,
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
  yticks <- seq(1, nrow(relabeled_target), length.out = length(yticklabels))
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

################################################################################
#                     Main Parameter Parsing Function                          #
################################################################################

run_plot_function <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    stop("Usage: Rscript plot_script.R <function_name> [parameters...]")
  }
  
  function_name <- args[1]
  parameters <- args[-1]
  
  # Parameter requirements for each function
  param_requirements <- list(
    "resplot_mash_up_all" = c("OA", "OB", "RB", "coloresB", 
                              "lower_limit", "upper_limit", "path", "name"),
    "resplot_superpop" = c("OA", "RB", "OA_colors", "pop_info",
                           "lower_limit", "upper_limit", "path", "name"),
    "resplot_mash_up_slice"= c("metricA","metricB","RA","RB","OA", "OB","coloresA","coloresB",
                               "RA_lower_limit","RA_upper_limit","RB_lower_limit",
                               "RB_upper_limit", "RA_steps","RB_steps","path","plotname","plotnameA",
                               "plotnameB","select_by","RA_val","RB_val"),
    "relabel_order_colors" = c("matrix", "target", "color_palette", "path", "name", 
                               "upper_limit", "lower_limit")
    
  )
  
  # Check if function exists
  if (!exists(function_name)) {
    stop(paste("Function", function_name, "does not exist"))
  }
  
  # Check parameter count
  required_params <- param_requirements[[function_name]]
  if (length(parameters) != length(required_params)) {
    stop(paste("Function", function_name, "requires", length(required_params), 
               "parameters:", paste(required_params, collapse = ", ")))
  }
  
  # Load parameters with type conversion
  load_parameter <- function(path) {
    if (endsWith(path, ".txt") || endsWith(path, ".csv")) {
      
      # parámetros que SÍ tienen header
      header_params <- c(
        "pop_info"
      )
      
      use_header <- param_name %in% header_params
      
      return(
        read.csv2(
          path,
          sep = ifelse(endsWith(path, ".csv"), ",", "\t"),
          header = use_header
        )
      )
    }
    return(path)
  }
  
  # Prepare parameter list
  params <- list()
  for (i in seq_along(required_params)) {
    param_name <- required_params[i]
    param_value <- parameters[i]
    
    # Special handling for numeric parameters
    if (param_name %in% c("lower_limit", "upper_limit")) {
      params[[param_name]] <- as.numeric(param_value)
    } else {
      params[[param_name]] <- load_parameter(param_value)
    }
  }
  
  # Call the appropriate function
  do.call(function_name, params)
}





################################################################################
#                                Execution                                    #
################################################################################

# Only run if executed as a script (not when sourced)
if (sys.nframe() == 0) {
  tryCatch({
    run_plot_function()
  }, error = function(e) {
    message("Error: ", e$message)
    quit(status = 1)
  })
}
