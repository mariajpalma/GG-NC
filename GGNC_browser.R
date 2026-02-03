# Interactive map 
# Version: 7

# Libraries: 
library(shiny)
library(shinytitle)
library(pracma)
library(leaflet)
library(plotly)
library(bslib)
library(igraph)
library(shinyalert)
library(shinyWidgets)
library(shinyBS)
library(rgl)
library(shinyscreenshot)
library("DescTools")
library("vembedr")
library("leaflet.minicharts")
library(shinybusy)

### Load functions
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
      cl <- cluster_leiden(graph, resolution = 1)
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
    if (!is_igraph(graph)) {
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

color_BinA <- function(step_A, step_B, use_leiden, mapbig_A, mapbig_B, lay_A, lay_B, ord_A, ord_B, im_A, im_B, net_A, net_B){
  # We color the networks with your own information
  net_A_step <- graph_louvain_by_comm(ord_A, im_A, step_A, net_A, mapbig_A, use_leiden= use_leiden)
  net_B_step <- graph_louvain_by_comm(ord_B, im_B, step_B, net_B, mapbig_B, use_leiden= use_leiden)
  
  # We execute the function that colors network A like network B
  red_A <- net_A_step
  red_B <- net_B_step
  lay <- lay_A
  
  # Get node names in net_ B
  nodos_B <- V(red_B)$name
  
  # Initialize colors for network A (default gray)
  colores_A_actualizado <- rep("#9b9b9b", length(V(red_A)))
  
  # Assign colors from network B to matching nodes in A
  for (i in 1:length(nodos_B)) {
    idx <- which(V(red_A)$name == nodos_B[i])
    if (length(idx) > 0) {
      colores_A_actualizado[idx] <- V(red_B)$color[i]
    }
  }
  
  # Create another network so as not to modify the ones we already have.
  red_C <- red_A
  
  # Change colors in red_A
  V(red_C)$color_actualizado <- colores_A_actualizado
  
  return(red_C)
}

### Reading files: data for the charts: IBD/PCA/GRM metrics
IBD <- read.csv("IBD.txt", sep = "\t", header =  T)
IBD5 <- read.csv("IBD5.txt", sep = "\t", header =  T)
PCA <- read.csv("PCA.txt", sep = "\t", header =  T)
GRM_rare <- read.csv("GRM_rare.txt", sep = "\t", header =  T)
GRM_common <- read.csv("GRM_common.txt", sep = "\t", header =  T)
T2D <- read.csv("T2D.txt", sep = "\t", header =  T)
skin <- read.csv("Skin.txt", sep = "\t", header =  T) 
alt <- read.csv("Alt.txt", sep = "\t", header =  T)

### Reading files: Colors for the pie charts 
IBD_colors <- scan("colors_IBD.txt", character(), quote = "")
IBD5_colors <- scan("colors_IBD5.txt", character(), quote = "")
PCA_colors <- scan("colors_PCA.txt", character(), quote = "")
GRM_rare_colors <- scan("colors_GRM_rare.txt", character(), quote = "")
GRM_common_colors <- scan("colors_GRM_common.txt", character(), quote = "")
T2D_colors <- scan("colors_T2D.txt", character(), quote = "")
skin_colors <- scan("colors_Skin.txt", character(), quote = "")
alt_colors <- scan("colors_Alt.txt", character(), quote = "")

## Loading similar colors for each metric 
IBD_similar_colors <- scan("colors_similar_IBD.txt", character(), quote = "")
IBD5_similar_colors <- scan("colors_similar_IBD5.txt", character(), quote = "")
PCA_similar_colors <- scan("colors_similar_PCA.txt", character(), quote = "")
GRM_similar_rare_colors <- scan("colors_similar_GRM_rare.txt", character(), quote = "")
GRM_similar_common_colors <- scan("colors_similar_GRM_common.txt", character(), quote = "")
T2D_similar_colors <- scan("colors_similar_T2D.txt", character(), quote = "")
skin_similar_colors <- scan("colors_similar_Skin.txt", character(), quote = "")
alt_similar_colors <- scan("colors_similar_Alt.txt", character(), quote = "")

### Reading resolution plot info for each metric 
IBD_resolution <- read.csv2("resolution_IBD.txt", sep = ",", header = F)
IBD5_resolution <- read.csv2("resolution_IBD5.txt", sep = ",", header = F)
PCA_resolution <- read.csv2("resolution_PCA.txt", sep = ",", header = F)
GRM_rare_resolution <- read.csv2("resolution_GRM_rare.txt", sep = ",", header = F)
GRM_common_resolution <- read.csv2("resolution_GRM_common.txt", sep = ",", header = F)
T2D_resolution <- read.csv2("resolution_T2D.txt", sep = ",", header = F)
skin_resolution <- read.csv2("resolution_Skin.txt", sep = ",", header = F)
alt_resolution <- read.csv2("resolution_Alt.txt", sep = ",", header = F)

## Load information for the networks 
# Distintive colors 
IBD_dis.network <- readRDS("IBD_3Dplots_distintive.rds")
IBD5_dis.network <- readRDS("IBD5_3Dplots_distintive.rds")
PCA_dis.network <- readRDS("PCA_3Dplots_distintive.rds")
GRM_common_dis.network <- readRDS("GRMc_3Dplots_distintive.rds")
GRM_rare_dis.network <- readRDS("GRMr_3Dplots_distintive.rds")
T2D_dis.network <- readRDS("T2D_3Dplots_distintive.rds")
skin_dis.network <- readRDS("Skin_3Dplots_distintive.rds")
alt_dis.network <- readRDS("Alt_3Dplots_distintive.rds")

# Similar colors
IBD_sm.network <- readRDS("IBD_3Dplots_similar.rds")
IBD5_sm.network <- readRDS("IBD5_3Dplots_similar.rds")
PCA_sm.network <- readRDS("PCA_3Dplots_similar.rds")
GRM_common_sm.network <- readRDS("GRMc_3Dplots_similar.rds")
GRM_rare_sm.network <- readRDS("GRMr_3Dplots_similar.rds")
T2D_sm.network <- readRDS("T2D_3Dplots_similar.rds")
skin_sm.network <- readRDS("Skin_3Dplots_similar.rds")
alt_sm.network <- readRDS("Alt_3Dplots_similar.rds")

# graphml files 
IBD_graphml <- read_graph("graph_by_comm_IBD_v3.graphml", format = "graphml")
IBD5_graphml <- read_graph("graph_by_comm_IBD5.graphml", format = "graphml")
PCA_graphml <- read_graph("graph_by_comm_PCA_v3.graphml", format = "graphml")
GRM_common_graphml <- read_graph("graph_by_comm_GRM_common_v3.graphml", format = "graphml")
GRM_rare_graphml <- read_graph("graph_by_comm_GRM_rare_v3.graphml", format = "graphml")

# Layout
IBD_layout<- readRDS("IBD_layout.rds")
IBD5_layout<- readRDS("IBD5_layout.rds")
PCA_layout <- readRDS("PCA_layout.rds")
GRM_common_layout <- readRDS("GRM_common_layout.rds")
GRM_rare_layout <- readRDS("GRM_rare_layout.rds")

# Index files 
IBD_index <- read.csv2("IBD_index.txt", sep = ",", header = F)
IBD5_index <- read.csv2("IBD5_index.txt", sep = ",", header = F)
PCA_index <- read.csv2("PCA_index.txt", sep = ",", header = F)
GRM_common_index <- read.csv2("GRM_common_index.txt", sep = ",", header = F)
GRM_rare_index <- read.csv2("GRM_rare_index.txt", sep = ",", header = F)

# Lambda values 
PCAlv <- logspace(-3,0,50)
IBDlv <- logspace(-3,1,50)
IBD5lv <- logspace(-5,0,50)
GRMclv <- logspace(-5,-1,50)
GRMrlv <- logspace(-6,-2,50)

# By default we load the data from IBD metrics/distintive colors: 
datos_ciudades <- IBD
colors_name <- IBD_colors
resolution_plot <- IBD_resolution
network_data <- IBD_dis.network

# For network panel, by default NetworkA = IBD
graphml.netA <- IBD_graphml
layout.netA <- IBD_layout
index_file.netA <- IBD_index
netA_colors <- IBD_colors
resolutionA <- IBD_resolution
lower_limit <- -3
upper_limit <- 1 

# For network panel, by default NetworkA = GRM_common
graphml.netB <- GRM_common_graphml
layout.netB <- GRM_common_layout
index_file.netB <- GRM_common_index
netB_colors <- GRM_common_colors
resolutionB <- GRM_common_resolution

ui <- fluidPage(
  tags$head(tags$meta(property = "og:image", content = "prev.png"),),
  title = "GG-NC",
  navbarPage("Global Genomic Network Communities (GG-NC) Browser",
    tabPanel("Map",
      sidebarLayout(
        sidebarPanel(
          selectInput("metrica", "Select a genetic metric:", 
                    choices= c("Identity By Descent > 2 cM", 
                               "Identity By Descent > 5 cM",
                                "Principal Component Analysis",
                                "Genetic Relationship Matrix based on rare variants", 
                                "Genetic Relationship Matrix based on common variants",
                                "PCA type 2 diabetes",
                                "PCA altitude adaptation",
                                "PCA skin pigmentation"
                               )),
        # Add tooltip
        bsTooltip("metrica", "Choose from the provided genetic metrics: Identity by Descent (IBD), Principal Component Analysis (PCA), or Genetic Relationship Matrix (GRM) for rare or common variants. The selected metric is utilized to construct a network, and community detection algorithms are subsequently applied.",
                    "right", options = list(container = "body")),
        selectInput("ubicacion", "Select a cohort:",
                    choices = unique(datos_ciudades$Pop)),
        # Add tooltip
        bsTooltip("ubicacion", "Choose a cohort from either the 1000 Genomes Project or the Human Genome Diversity Project, and the corresponding pie chart will be presented.",
                  "right", options = list(container = "body")),
        tags$head(tags$style(HTML('.irs-from, .irs-to, .irs-min, .irs-max, .irs-single {
                                    visibility: hidden !important;
        }'))),
        sliderInput("slider_id", label = "Select a resolution value:", min = 1, max = 50, value = 25,step = 1, ticks = FALSE), 
        bsTooltip("slider_id", "Select one of the 50 resolution values. Lower values lean towards larger communities, whereas higher values lean towards smaller communities",
                  "right", options = list(container = "body")),
        # Value_text for lambda value
        verbatimTextOutput("value_text"),
        verbatimTextOutput("community.num"),
        br(),
      
        imageOutput("resolution_plot"), br(),
        screenshotButton(label= "Download plot", id= "resolution_plot",filename = "Resolution_plot"), br(),
        tags$head(tags$style(".butt{background-color:white;} .butt{color: black;} .butt{border: none}")),
        br(), br(), br(), br(), br(), br(), br(),
        useShinyalert(force = TRUE),
        
        tags$div(
          title = "Opting for this feature will present the resolution plot such that more 
          similar communities (in a genetic sense) 
          will be represented with more similar colors, promoting visual coherence.",
          shinyWidgets::switchInput(inputId = "mySwitch",      
                                    value = FALSE,
                                    label = "Similar colors")
        ),
        actionButton("help", "Help", class = "h"),
        tags$head(tags$style(".h{background-color:white;} .h{color: black;} .h{border: none}"))
        ),
        mainPanel(
        leafletOutput("mapa"),
        tags$head(tags$style(".leaflet-popup-content-wrapper{box-shadow:none !important; }")),
        tags$head(tags$style(".leaflet-popup-tip{box-shadow: none;}")),
          
        fluidRow(br(),
          column(width = 12, align = "right",screenshotButton(label= "Download map", id= "mapa", filename = "GG-NC"))),
          
        fluidRow(div(style = "display:flex; flex-wrap:wrap; justify-content:space-between; align-items:center; gap:30px;", #width:100%;
          div(style = "flex:1 1 60%; min-width:300px;",plotlyOutput("diagramaPie", height = "400px")),
          div(style = "flex:1 1 35%; min-width:300px; text-align:center;", 
                div(style = "margin:auto 0;",
                    rglwidgetOutput("network", width = "400px", height = "400px")
                )
              )
            ),
            br(), br(), br(), br(), br(), br(), br(), br(), br() 
          ),

          tags$footer(style = "
                            width: 100%;
                            background-color: white;
                            padding: 15px 0;
                            margin-top: 30px;
                            border-top: 1px solid #ddd;
                            text-align: center;
                            position: relative;
                            bottom: 0;
                          ",
            
            # --- Primera fila: enlace al laboratorio ---
            div(style = "margin-bottom: 8px;",
              a(href = 'https://www.sohaillab.com/contact', 'Sohail-Lab', target = "_blank",
                style = 'color:#2c3e50; font-weight:bold; text-decoration:none;')
            ),
            
            # --- Segunda fila: autores ---
            div(style = "font-size:10px; margin-bottom: 8px;",
              "Developed by María J. Palma-Martínez, Yuridia S. Posadas-García, Brenda E. López-Ángeles, Claudia Quiroz-López, Diego Ortega-Del Vecchyo, Anna C. F. Lewis, Kevin A. Bird, Tina Lasisi, Arslan A. Zaidi, and Mashaal Sohail"
            ),
            
            # --- Tercera fila: logos alineados ---
            div(
              style = "display:flex; flex-wrap:wrap; justify-content:center; align-items:center; gap:15px; margin-bottom: 10px;",
              img(height = 70, width = 160,
                  src = "https://www.ciatec.mx/images/paginas/daafc7198552967c65d2947c39e6cad1.jpg"),
              img(height = 70, width = 80,
                  src = "https://oferta.unam.mx/assets/img/dummies/ccg1.png"),
              img(height = 40, width = 35,
                  src = "https://i0.wp.com/www.atmosfera.unam.mx/wp-content/uploads/2019/06/unam-escudo-azul.png?ssl=1")
            ),
            
            # --- Cuarta fila: financiamiento ---
            div(
              style = "font-size:10px; margin-bottom: 8px;",
              a(
                href = "https://www.ciatec.mx/pagina/6/14-informes-secihti",
                "Funding provided by SECIHTI",
                target = "_blank",
                style = "color:#2c3e50; text-decoration:none;"
              )
            ),
            
            # --- Quinta fila: licencia CC ---
            div(a(href = "https://creativecommons.org/licenses/by/4.0/?ref=chooser-v1",
                target = "_blank",
                img(height = 20, width = 20,
                src = "https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1")
              )
            ),
            
            div(style = "font-size:10px; margin-bottom: 8px;",
                "This app is designed for use on a desktop or laptop computer."
            )
          ) # footer
        )
      )
    ), # Tabpanel: networks 
    tabPanel("Network of individuals",
            sidebarLayout(
              sidebarPanel(
                  # Select a metric 
                  selectInput("metrica_network", "Network of individuals based on:", 
                              choices= c("Identity By Descent > 2 cM", 
                                         "Identity By Descent > 5 cM",
                                         "Principal Component Analysis",
                                         "Genetic Relationship Matrix (rare variants)", 
                                         "Genetic Relationship Matrix (common variants)")),
                  # Add tooltip
                  shinyBS::bsTooltip("metrica_network", "Choose the genetic metric to display the individuals network (left panel).",
                            "right", options = list(container = "body")),
                  # if wanted 
                  selectInput("color_metric", "Color nodes based on:", 
                              choices= c("Identity By Descent > 2 cM", 
                                        "Identity By Descent > 5 cM",
                                         "Principal Component Analysis",
                                         "Genetic Relationship Matrix (rare variants)", 
                                         "Genetic Relationship Matrix (common variants)")),
                  # Add tooltip
                  bsTooltip("color_metric", "Choose a genetic metric to color the previously selected network of individuals (left). The network’s community composition at each selected resolution will be shown based on this metric. Individuals absent from this metric will appear in gray. The network of individuals from this selected metric will appear in the right panel.",
                            "right", options = list(container = "body")),
                  
                  sliderInput("slider_id_network", label = "Select a resolution value:", min = 1, max = 50, value = 25,step = 1, ticks = FALSE), 
                  bsTooltip("slider_id_network", "This slider lets you explore the community composition of the second selected metric across different resolution values, projected onto a network generated from the first selected metric (left panel).",
                            "right", options = list(container = "body")),
                  # Value_text for lambda value
                  verbatimTextOutput("value_text.network")
                  
              ), # Main panel
              mainPanel(
                fluidRow(
                  column(6,
                         uiOutput("colorby.net"),
                         plotOutput("networkA")
                  ),
                  column(6,
                         uiOutput("display.net"),
                         plotOutput("networkB")
                  )
                )
              )
            )
    ),
    # tabpanel: map
    tabPanel("Tutorial",
             mainPanel(
               br(),
               column( 12, align= "right", embed_url("https://youtu.be/2bdod1RuRVk")),
               br(),
               column( 12, align= "right", embed_url("https://youtu.be/mVegRKbWSAI"))
            )
             ),
    tabPanel("Customize", 
             sidebarLayout(
               sidebarPanel(
                 textInput("lambdavalue.customize", label= "Insert a resolution value index", value = "1", width = NULL),
                 verbatimTextOutput("value.customize"),
                 verbatimTextOutput("possible.customize"),
                 verbatimTextOutput("community.customize"),
                 fileInput("shiny.info", "Choose a shiny file", accept = ".txt"), 
                 tags$head(tags$style(HTML('.irs-from, .irs-to, .irs-min, .irs-max, .irs-single {
                                    visibility: hidden !important;}'))),
                 fileInput("shiny.col", "Choose a shiny colors file", accept = ".txt"),
                 fileInput("resolution.data", "Choose a resolution plot file", accept = ".txt"),
                 imageOutput("resolution.plot.customize"),
                 actionButton("help.custumize", "Help", class = "c")
               ),
               mainPanel(
                 leafletOutput("mapa.customize")
               )
             )
    )
  )
)

server <- function(input, output, session) {
  
  ##### Main Panel
  # Help message 
  observeEvent(input$help, {
    shinyalert(
      title = "Help",
      text = HTML("
Welcome to the Global Genomic Network Communities (GG-NC) Browser!<br><br>

How do we figure out who's genetically similar to who? Think of each person's genome like its own unique city. The viewpoint you choose and what you decide to focus on can completely change how you see its similarity to another city. Are you looking from a helicopter focusing on skyscrapers, or are you on the street counting parks? That's the kind of choice the Global Genomic Network Communities (GG-NC) Browser offers. You can select different metrics—like Identity by Descent (IBD), Principal Component Analysis (PCA), or Genetic Relationship Matrix (GRM)—to change your 'viewpoint' and focus on different 'landmarks' to compare two genomes. Dive in and see how the story changes when you start looking at things differently!<br><br>

To begin you can pick a genetic metric and a resolution value in the left panel. So far you can choose between three metrics: Identity by Descent (IBD), Principal Component Analysis (PCA), and Genetic Relationship Matrix (GRM). The fifty resolution values are in the logarithmic space (base 10) from -2 to 2. A network is computed using the genetic metric chosen, and communities are detected on that network at the chosen resolution using Louvain’s algorithm.<br><br>

Any pie chart shows the individuals in that cohort that belong to the network communities detected at that resolution. You can play around and see how the communities change across these values. On the bottom left, you will see a graph showing the communities detected across our range of resolution values, colored in the same way as the pie charts. The individuals in this graph are ordered according to the communities detected, and can also be observed on the world map. As a note, any community with less than 6 members is colored in white.<br><br>

The available cohorts come from the 1000 Genomes Project and the Human Genome Diversity Project. To display a cohort’s name, hover over a pie chart on the map. If you wish to take a closer look into a particular cohort, select its name in the left panel. You can also zoom in on the map for a closer view of a particular geographic area.<br><br>

This is still a project in progress and modifications will be included in the near future in new versions. Questions and comments are very welcome!<br><br>

Contact:<br>
mashaal at ccg dot unam dot mx<br>
mpalma at ccg dot unam dot mx<br><br>

If you use the maps in other works, include the source link:<br>
<a href='http://sohaillab.ccg.unam.mx/gg-nc/' target='_blank'>http://sohaillab.ccg.unam.mx/gg-nc//</a><br><br>

Find our pre-print <a href='https://www.biorxiv.org/content/10.1101/2024.12.11.627824v1' target='_blank'>here</a>.<br><br>

Version 0.2 (beta)
"),
      type = "info",
      html = TRUE
    )
  })
  
  # Obvserve event when selecting a population.
  observeEvent(input$ubicacion, {
    pop_selected = datos_ciudades[datos_ciudades$Pop == input$ubicacion,]
    proxy <- leafletProxy('mapa')
    
    leafletProxy("mapa") %>%
      addPopups(
        lng=pop_selected$Longitud, 
        lat=pop_selected$Latitud,
        pop_selected$Pop,
        options = popupOptions(closeButton = FALSE)
      )
  })
 
 ##############################################################################################
  # --- MAPA INICIAL (solo una vez) ---
  output$mapa <- renderLeaflet({
    leaflet() %>%
      addTiles()
  })
  
  
  # --- ACTUALIZACIÓN CUANDO CAMBIA CUALQUIER INPUT ---
  observeEvent({
    input$slider_id
    input$metrica
    input$mySwitch
  }, {
    
    # ---- 1. Seleccionar dataset según métrica ----
    
    if(input$metrica == "Identity By Descent > 5 cM"){
      datos_ciudades <- IBD5
      resolution_plot <- IBD5_resolution
      lower_limit <- -5
      upper_limit <- 0
      colors_name <- if (input$mySwitch) IBD5_similar_colors else IBD5_colors
    }
    
    if(input$metrica == "Identity By Descent > 2 cM"){
      datos_ciudades <- IBD
      resolution_plot <- IBD_resolution
      lower_limit <- -3
      upper_limit <- 1
      colors_name <- if (input$mySwitch) IBD_similar_colors else IBD_colors
    }
    
    if(input$metrica=="Principal Component Analysis"){
      datos_ciudades <- PCA
      resolution_plot <- PCA_resolution
      lower_limit <- -3
      upper_limit <- 0
      colors_name <- if (input$mySwitch) PCA_similar_colors else PCA_colors
    }
    
    if(input$metrica=="Genetic Relationship Matrix based on rare variants"){
      datos_ciudades <- GRM_rare
      resolution_plot <- GRM_rare_resolution
      lower_limit <- -6
      upper_limit <- -2
      colors_name <- if (input$mySwitch) GRM_similar_rare_colors else GRM_rare_colors
    }
    
    if(input$metrica=="Genetic Relationship Matrix based on common variants"){
      datos_ciudades <- GRM_common
      resolution_plot <- GRM_common_resolution
      lower_limit <- -5
      upper_limit <- -1
      colors_name <- if (input$mySwitch) GRM_similar_common_colors else GRM_common_colors
    }
    
    if(input$metrica=="PCA type 2 diabetes"){
      datos_ciudades <- T2D
      resolution_plot <- T2D_resolution
      lower_limit <- -3
      upper_limit <- 0
      colors_name <- if (input$mySwitch) T2D_similar_colors else T2D_colors
    }
    
    if(input$metrica=="PCA altitude adaptation"){
      datos_ciudades <- alt
      resolution_plot <- alt_resolution
      lower_limit <- -2
      upper_limit <- 0
      colors_name <- if (input$mySwitch) alt_similar_colors else alt_colors
    }
    
    if(input$metrica=="PCA skin pigmentation"){
      datos_ciudades <- skin
      resolution_plot <- skin_resolution
      lower_limit <- -3
      upper_limit <- 0
      colors_name <- if (input$mySwitch) skin_similar_colors else skin_colors
    }
    
    
    # ---- 2. Texto de resolución ----
    lambda <- sort(unique(datos_ciudades$Lambda))
    log_lambda <- lambda[input$slider_id]
    
    output$value_text <- renderText({
      paste("Resolution:", round(log10(log_lambda), 3))
    })
    
    
    # ---- 3. Plot de resolución ----
    output$resolution_plot <- renderPlot({
      
      im2 <- apply(resolution_plot, 2, rev)
      x_vals <- 1:ncol(resolution_plot)
      y_vals <- 1:nrow(resolution_plot)
      
      image(
        x = x_vals,
        y = y_vals,
        z = t(im2),
        useRaster = TRUE,
        col = colors_name,
        xaxs = "i",
        yaxt = "n",
        xlab = "Individuals",
        ylab = " "
      )
      
      mtext(
        side = 2, line = 2,
        text = "Resolution value"
      )
      
      yticklabels <- rev(lower_limit:upper_limit)
      yticks <- seq(1, nrow(resolution_plot), length.out = length(yticklabels))
      
      axis(2, at = yticks, labels = yticklabels)
      
      # Línea seleccionada
      selected_value <- (50:1)[input$slider_id]
      abline(h = selected_value, lwd=3, lty=2)
    })
    
    
    # ---- 4. Filtrar datos para el lambda seleccionado ----
    data <- datos_ciudades[datos_ciudades$Lambda == log_lambda,]
    
    c1 <- which(colnames(datos_ciudades)=="C1")
    lcolumn <- ncol(datos_ciudades)
    
    communities.info <- data[,c1:lcolumn]
    number.com <- length(which(colSums(!is.na(communities.info)) > 0))
    
    output$community.num <- renderText({
      paste("Total of detected communities:", number.com)
    })
    
    show_modal_spinner(text = HTML("<br>Loading...<br>Higher resolution values may take longer to load."))
    
    # ---- 5. ACTUALIZAR EL MAPA SIN RECREARLO ----
    leafletProxy("mapa") %>%
      clearMinicharts() %>%
      addMinicharts(
        lat = data$Latitud,
        lng = data$Longitud,
        type = "pie",
        colorPalette = colors_name,
        legend = FALSE,
        popup = popupArgs(html = paste0(data$Pop)),
        #chartdata = data[,c1:lcolumn]
        chartdata = round(data[,c1:lcolumn], 3)
      )
    
    session$onFlushed(function() {
      remove_modal_spinner()
    }, once = TRUE)
  })

###########################################################################################33  
 
  # Add tooltip for the resolution plot
  addPopover(session,"resolution_plot", "Resolution Plot", 
             content = "This plot displays the identified communities at various resolution values. 
             The y-axis represents the resolution values, 
               while the x-axis denotes the index of individuals. The colors in this plot correspond 
               to those displayed in the accompanying map. Communities that were less tan six members or were dected in only one 
             resolution are colored in white.", 
             placement = "right", options = list(container = "body"))
  
  # Add tooltip for the map
  addPopover(session,"mapa", "Map", content = "In this map, cohorts are spatially linked to their respective sampling locations. 
             Each pie chart corresponds to a cohort sourced from either the 1000 Genomes Project or the Human Genome Diversity Project. 
             The distinctive colors identified communities at the selected resolution value.", 
             placement = "bottom", options = list(container = "body"))
  
  output$diagramaPie <- renderPlotly({
    
    ## Select data set
    # For IBD
    if (input$mySwitch) {
      colors_name <- IBD_similar_colors
    }
    
    if(input$metrica=="Identity By Descent > 5 cM"){
      datos_ciudades <- IBD5
      if (input$mySwitch) {
        colors_name <- IBD5_similar_colors
      } else {
        colors_name <- IBD5_colors
      }
    }
    
    if(input$metrica=="Principal Component Analysis"){
      datos_ciudades <- PCA
      if (input$mySwitch) {
        colors_name <- PCA_similar_colors
      } else {
        colors_name <- PCA_colors
      }
    }
    
    if(input$metrica=="Genetic Relationship Matrix based on rare variants"){
      datos_ciudades <- GRM_rare
      if (input$mySwitch) {
        colors_name <- GRM_similar_rare_colors
      } else {
        colors_name <- GRM_rare_colors
      }
    }
    
    if(input$metrica=="Genetic Relationship Matrix based on common variants"){
      datos_ciudades <- GRM_common
      if (input$mySwitch) {
        colors_name <- GRM_similar_common_colors
      } else {
        colors_name <- GRM_common_colors
      }
    }
    
    if(input$metrica=="PCA type 2 diabetes"){
      datos_ciudades <- T2D
      if (input$mySwitch) {
        colors_name <- T2D_similar_colors
      } else {
        colors_name <- T2D_colors
      }
    }
    
    if(input$metrica=="PCA skin pigmentation"){
      datos_ciudades <- skin
      if (input$mySwitch) {
        colors_name <- skin_similar_colors
      } else {
        colors_name <- skin_colors
      }
    }
    
    if(input$metrica=="PCA altitude adaptation"){
      datos_ciudades <- alt
      if (input$mySwitch) {
        colors_name <- alt_similar_colors
      } else {
        colors_name <- alt_colors
      }
    }
    
    ubicacion_seleccionada <- input$ubicacion
    selected_value <- input$slider_id
    lambda <- sort(unique(datos_ciudades$Lambda))
    log_lambda <- lambda[selected_value]
    
    data <- datos_ciudades[datos_ciudades$Pop == ubicacion_seleccionada &
                             datos_ciudades$Lambda == log_lambda, ]
    
    # Number of columns 
    lcolumn <- length(datos_ciudades)
    
    # Possition of the first community
    c1 <- which(colnames(datos_ciudades)=="C1")
    
    # Selecting the data for all communities
    comunidades <- data[,c1:lcolumn]
    
    # Getting rid of the NA values 
    plotdata <- comunidades[!is.na(comunidades)]
    
    # Selecting colors 
    get_colors <- colors_name[!is.na(comunidades)]
    
    # Generating community names: labels  
    cvector <- 1:length(comunidades)
    clabels <- cvector[!is.na(comunidades)]
    clabels <- paste("Network community", clabels)
    
    plot_ly(labels= clabels,
            values = plotdata,
            sort = FALSE,
            marker = list(colors = get_colors, line = list(color = "black", width = 0.3)),
            type = 'pie', textfont = list(size = 15), width = 595, height = 595,
            textinfo= "none",
            hovertemplate = "%{label}<br>%{percent}<extra></extra>") %>% 
      layout(title = paste("<br>Proportion of individuals in each network community <br>Cohort:<br>",ubicacion_seleccionada, "\nData set:", data$Project,
                           "<br>Resolution value: ", round(log10(log_lambda),3)), legend = list(x = 1.1, y = 0.5))
  })
  
  # Add tooltip for the pie chart 
  addPopover(session, "diagramaPie", "Pie chart", content = "The pie chart illustrates the distribution of individuals within a 
             chosen cohort across various network communities.",
             placement = "right", options = list(container = "body"))
  
  # Network 
  output$network <- renderRglwidget({
    # Lambda selected 
    x <- input$slider_id
    
    # Select data set for the networks 
    # For IBD
    if (input$mySwitch) {
      network_data <- IBD_sm.network
    }
    
    if(input$metrica=="Identity By Descent > 5 cM"){
      if (input$mySwitch) {
        network_data <- IBD5_sm.network
      } else {
        network_data <- IBD5_dis.network
      }
    }
    
    if(input$metrica=="Principal Component Analysis"){
      if (input$mySwitch) {
        network_data <- PCA_sm.network
      } else {
        network_data <- PCA_dis.network
      }
    }
    
    if(input$metrica=="Genetic Relationship Matrix based on rare variants"){
      if (input$mySwitch) {
        network_data <- GRM_rare_sm.network
      } else {
        network_data <- GRM_rare_dis.network
      }
    }
    
    if(input$metrica=="Genetic Relationship Matrix based on common variants"){
      if (input$mySwitch) {
        network_data <- GRM_common_sm.network
      } else {
        network_data <- GRM_common_dis.network
      }
    }
    
    if(input$metrica=="Genetic Relationship Matrix based on common variants"){
      if (input$mySwitch) {
        network_data <- GRM_common_sm.network
      } else {
        network_data <- GRM_common_dis.network
      }
    }
    
    if(input$metrica=="PCA type 2 diabetes"){
      if (input$mySwitch) {
        network_data <- T2D_sm.network
      } else {
        network_data <- T2D_dis.network
      }
    }
    
    if(input$metrica=="PCA skin pigmentation"){
      if (input$mySwitch) {
        network_data <- skin_sm.network
      } else {
        network_data <- skin_dis.network
      }
    }
    
    if(input$metrica=="PCA altitude adaptation"){
      if (input$mySwitch) {
        network_data <- alt_sm.network
      } else {
        network_data <- alt_dis.network
      }
    }
    
    # Plot the network
    rglplot(network_data[[x]], layout = as.matrix(network_data[[x]]$coords), vertex.label = NA)
    rglwidget()
  })
  
  # Tooltip for the network
  addPopover(session,"network", "Community network", content = "This network is constructed by averaging the positions of individuals within each community in the network of individuals. The node sizes are proportional to the respective community sizes. 
              Edges indicate the density of connections between these communities. Communities in white are not shown here.",
              placement = "left", options = list(container = "body"))
  
  # Clear the existing plot when
  #   similar colors button on/off
  #   lambda value changes 
  #   metric changes 
  observeEvent(input$mySwitch, {
    clear3d()
  })
  
  observeEvent(input$slider_id, {
    clear3d()
  })
  
  observeEvent(input$metrica, {
    clear3d()
  })
  
  ##### Tutorial panel 
  
  ##### Network panel 
  
  # Print lambda value 
  output$value_text.network <- renderText({
    if (input$color_metric == "Principal Component Analysis") {
      paste("Resolution:", round(log10(PCAlv[input$slider_id_network]), 3))
    } else if (input$color_metric == "Genetic Relationship Matrix (rare variants)") {
      paste("Resolution:", round(log10(GRMrlv[input$slider_id_network]), 3))
    } else if (input$color_metric == "Genetic Relationship Matrix (common variants)") {
      paste("Resolution:", round(log10(GRMclv[input$slider_id_network]), 3))
    } else if (input$color_metric == "Identity By Descent > 2 cM") {
      paste("Resolution:", round(log10(IBDlv[input$slider_id_network]), 3))
    } 
    #else if (input$color_metric == "Identity By Descent > 5 cM") {
    #  paste("Resolution:", round(log10(IBD5lv[input$slider_id_network]), 3))
    #}
  })
  
  output$colorby.net <- renderUI({ 
  if(input$metrica_network == input$color_metric){
      h3(paste(input$metrica_network), "Network")
  } else {
      h3(paste(input$metrica_network, 
               "Network colored by communities detected on", input$color_metric, "Network"))
  }
  })
  
  output$display.net <- renderUI({ 
    if(input$metrica_network == input$color_metric){
      return(NULL)
    } else {
      h3(paste(input$color_metric, "Network"))
    }
  })
  
  output$networkA <- renderPlot({
    show_modal_spinner(text = HTML("<br>Loading..."))
    selected_value <- input$slider_id_network
    
    # Network A
     if(input$metrica_network== "Identity By Descent > 5 cM"){
        graphml.netA <- IBD5_graphml
        layout.netA <- IBD5_layout
        index_file.netA <- IBD5_index
        netA_colors <- IBD5_colors
        resolutionA <- IBD5_resolution
      }
      if(input$metrica_network== "Principal Component Analysis"){
        graphml.netA <- PCA_graphml
        layout.netA <- PCA_layout
        index_file.netA <- PCA_index
        netA_colors <- PCA_colors
        resolutionA <- PCA_resolution
      }
      if(input$metrica_network== "Genetic Relationship Matrix (rare variants)"){
        graphml.netA <- GRM_rare_graphml
        layout.netA <- GRM_rare_layout
        index_file.netA <- GRM_rare_index
        netA_colors <- GRM_rare_colors
        resolutionA <- GRM_rare_resolution
      }
      if(input$metrica_network== "Genetic Relationship Matrix (common variants)"){
        graphml.netA <- GRM_common_graphml
        layout.netA <- GRM_common_layout
        index_file.netA <- GRM_common_index
        netA_colors <- GRM_common_colors
        resolutionA <- GRM_common_resolution
      }
    
    if(input$metrica_network == input$color_metric){
      netA <- graph_louvain_by_comm(index_file.netA, resolutionA, selected_value, graphml.netA, as.data.frame(netA_colors), TRUE)
    } else {
      # Color by other metric 
      if(input$color_metric == "Identity By Descent > 5 cM"){
        graphml.netB <- IBD5_graphml
        layout.netB <- IBD5_layout
        index_file.netB <- IBD5_index
        netB_colors <- IBD5_colors
        resolutionB <- IBD5_resolution
      }
      if(input$color_metric == "Principal Component Analysis"){
        graphml.netB <- PCA_graphml
        layout.netB <- PCA_layout
        index_file.netB <- PCA_index
        netB_colors <- PCA_colors
        resolutionB <- PCA_resolution
      }
      if(input$color_metric== "Genetic Relationship Matrix (rare variants)"){
        graphml.netB <- GRM_rare_graphml
        layout.netB <- GRM_rare_layout
        index_file.netB <- GRM_rare_index
        netB_colors <- GRM_rare_colors
        resolutionB <- GRM_rare_resolution
      }
      if(input$color_metric== "Identity By Descent > 2 cM"){
        graphml.netB <- IBD_graphml
        layout.netB <- IBD_layout
        index_file.netB <- IBD_index
        netB_colors <- IBD_colors
        resolutionB <- IBD_resolution
      }
      netA <- color_BinA(selected_value, 
                         selected_value, 
                         TRUE, 
                         mapbig_A=as.data.frame(netA_colors), 
                         mapbig_B=as.data.frame(netB_colors), 
                         lay_A= layout.netA, 
                         lay_B= layout.netB, 
                         ord_A= index_file.netA, 
                         ord_B= index_file.netB, 
                         im_A= resolutionA, 
                         im_B= resolutionB, 
                         net_A= graphml.netA, 
                         net_B= graphml.netB)
    }
    
    vertex_colors <- if ("color_actualizado" %in% vertex_attr_names(netA)) {
      V(netA)$color_actualizado
    } else {
      V(netA)$color
    }
    
    plot(netA, 
         layout = layout.netA,
         vertex.size = 6, 
         vertex.label.dist=1,
         vertex.label = NA, 
         vertex.color = vertex_colors,
         edge.color = "gray"
    )
    
    session$onFlushed(function() {
      remove_modal_spinner()
    }, once = TRUE)
    
  })
  
  output$networkB <- renderPlot({
    selected_value <- input$slider_id_network
    
    if (input$metrica_network == input$color_metric) return(NULL)
    
    if(input$color_metric== "Identity By Descent > 5 cM"){
      graphml.netB <- IBD5_graphml
      layout.netB <- IBD5_layout
      index_file.netB <- IBD5_index
      netB_colors <- IBD5_colors
      resolutionB <- IBD5_resolution
    }
    if(input$color_metric== "Principal Component Analysis"){
        graphml.netB <- PCA_graphml
        layout.netB <- PCA_layout
        index_file.netB <- PCA_index
        netB_colors <- PCA_colors
        resolutionB <- PCA_resolution
      }
    if(input$color_metric== "Genetic Relationship Matrix (rare variants)"){
        graphml.netB <- GRM_rare_graphml
        layout.netB <- GRM_rare_layout
        index_file.netB <- GRM_rare_index
        netB_colors <- GRM_rare_colors
        resolutionB <- GRM_rare_resolution
      }
    if(input$color_metric== "Identity By Descent > 2 cM"){
        graphml.netB <- IBD_graphml
        layout.netB <- IBD_layout
        index_file.netB <- IBD_index
        netB_colors <- IBD_colors
        resolutionB <- IBD_resolution
      }

    netB <- graph_louvain_by_comm(index_file.netB, resolutionB, selected_value, graphml.netB, as.data.frame(netB_colors), TRUE)
    plot(netB, 
         layout = layout.netB,
         vertex.size = 6, 
         vertex.label.dist=1,
         vertex.label = NA, 
         vertex.color = V(netB)$color,
         edge.color = "gray",
    )
  })
  
  ##### Customize panel 

  output$mapa.customize <- renderLeaflet({
    leaflet(options = leafletOptions(zoomControl = TRUE)) %>%
      addTiles() %>%  # Add default OpenStreetMap map tiles
      setView(lng = 0, lat = 0, zoom = 2) 
  })
  first_render <- TRUE
  observeEvent(input$lambdavalue.customize, {
    if (is.null(input$shiny.info) || is.null(input$shiny.col) || is.null(input$resolution.data)){
      return(NULL)
    }
    
    # Archivo shiny_info
    observeEvent(input$shiny_info, {
      #req(input$shiny_info)
      
      if (input$shiny_info$size > 4*1024*1024) {
        showNotification(
          paste("File", input$shiny_info$name, "exceeds maximum allowed size"),
          type = "error",
          duration = 8
        )
        return(NULL)
      }
    }, ignoreInit = TRUE)
    
    # Archivo shiny.col
    observeEvent(input$shiny.col, {
      #req(input$shiny.col)

      if (input$shiny.col$size > 3*1024) {
        showNotification(
          paste("File", input$shiny.col$name, "exceeds maximum allowed size"),
          type = "error",
          duration = 8
        )
        return(NULL)
      }
    }, ignoreInit = TRUE)
    
    # Archivo resolution.data
    observeEvent(input$resolution.data, {
      #req(input$resolution.dat)
  
      if (input$resolution.data$size > 1000*1024) {
        showNotification(
          paste("File", input$resolution.data$name, "exceeds maximum allowed size"),
          type = "error",
          duration = 8
        )
        return(NULL)  # ← Detiene la ejecución aquí
      }
    }, ignoreInit = TRUE)
    
    # Lambda
    observeEvent(input$lambdavalue.customize, {
      req(input$lambdavalue.customize)
      lambda_val <- as.numeric(input$lambdavalue.customize)
      
      # Solo pasa si no es entero
      if (is.na(lambda_val) || lambda_val %% 1 != 0) {
        showNotification(
          "Please introduce integers only",
          type = "error",
          duration = 8,
          id = "lambda_notification"
        )
        return(NULL)
      }
    }, ignoreInit = TRUE)
    
    
    file <- input$shiny.info
    datos_ciudades <- read.csv(file$datapath, sep = "\t", header = T)
    
    file2 <- input$shiny.col
    colors_name <- scan(file2$datapath, character(), quote = "")
    
    file3 <- input$resolution.data
    resolution_plot <- read.csv2(file3$datapath, sep = ",", header = F)
    
    selected_value <- as.numeric(input$lambdavalue.customize)
    lambda <- sort(unique(datos_ciudades$Lambda))
    log_lambda <- lambda[selected_value]
    
    data <- datos_ciudades[datos_ciudades$Lambda == log_lambda,]
    
    # Get the number of columns
    lcolumn <- length(datos_ciudades)
    
    # Obtaining the position of the first community
    if ("C1" %in% colnames(datos_ciudades)){
      c1 <- which(colnames(datos_ciudades)=="C1")
    } else {
      c1 <- which(colnames(datos_ciudades)=="C2")
    }
    
    # Printing number of communities 
    communities.info <- data[, c1:lcolumn]
    number.com <- length(colnames(communities.info)[apply(communities.info, 2, function(col) any(!is.na(col)))])
    
    output$community.customize <- renderText({
      paste("Total of detected communities:", number.com)
    })
    
    # Use leafletProxy to update the existing map
    leafletProxy("mapa.customize") %>%
      clearMarkers() %>%
      addMinicharts(
        lat = data$Latitud,
        lng = data$Longitud,
        type = "pie",
        colorPalette = colors_name,
        legend = FALSE,
        popup = popupArgs(html = paste0(data$Pop)),
        chartdata = data[, c1:lcolumn]
      )
    
    # Update other UI elements
    output$value.customize <- renderText({
      paste("Resolution:", round(log10(log_lambda), 3))
    })
    
    output$possible.customize <- renderText({
      paste("Resolution value range indexes: 1 -", length(lambda))
    })
    
    output$resolution.plot.customize <- renderPlot({
      im2 <- apply(resolution_plot, 2, rev)
      x_vals <- 1:ncol(resolution_plot)
      y_vals <- 1:nrow(resolution_plot)
      
      image(
        x = x_vals,
        y = y_vals,
        z = t(im2),
        useRaster = TRUE,
        col = colors_name,
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
        cex = 1
      )
      yticklabels <- rev(round(log10(min(lambda)), 3):round(log10(max(lambda)), 3))
      yticks <- seq(1, nrow(resolution_plot), length.out = length(yticklabels))
      
      axis(
        side = 2,
        at = yticks,
        labels = yticklabels,
        las = 1
      )
      
      abline(h = (length(lambda):1)[selected_value], col = "black", lwd = 3, lty = 2)
    })
  })
  
  # Help message 
  observeEvent(input$help.custumize, {
    shinyalert("Help", "We will be providing custom software to analyze your own genetic dataset. This will allow you to generate files you can upload on this page to visualize your own dataset in the Global Genomic Network Communities (GG-NC) browser.
The custom code will be available shortly at: 
\n https://github.com/mariajpalma/GG-NC \n

Please stay tuned!
", 
               type = "info")})
}

# Calling shiny app 
shinyApp(ui, server)