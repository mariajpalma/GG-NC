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

### Reading files: data for the charts: IBD/PCA/GRM metrics
IBD <- read.csv("IBD.txt", sep = "\t", header =  T)
PCA <- read.csv("PCA.txt", sep = "\t", header =  T)
GRM_rare <- read.csv("GRM_rare.txt", sep = "\t", header =  T)
GRM_common <- read.csv("GRM_common.txt", sep = "\t", header =  T)

### Reading files: Colors for the pie charts 
IBD_colors <- scan("colors_IBD.txt", character(), quote = "")
PCA_colors <- scan("colors_PCA.txt", character(), quote = "")
GRM_rare_colors <- scan("colors_GRM_rare.txt", character(), quote = "")
GRM_common_colors <- scan("colors_GRM_common.txt", character(), quote = "")

## Loading similar colors for each metric 
IBD_similar_colors <- scan("colors_similar_IBD.txt", character(), quote = "")
PCA_similar_colors <- scan("colors_similar_PCA.txt", character(), quote = "")
GRM_similar_rare_colors <- scan("colors_similar_GRM_rare.txt", character(), quote = "")
GRM_similar_common_colors <- scan("colors_similar_GRM_common.txt", character(), quote = "")

### Reading resolution plot info for each metric 
IBD_resolution <- read.csv2("resolution_IBD.txt", sep = ",", header = F)
PCA_resolution <- read.csv2("resolution_PCA.txt", sep = ",", header = F)
GRM_rare_resolution <- read.csv2("resolution_GRM_rare.txt", sep = ",", header = F)
GRM_common_resolution <- read.csv2("resolution_GRM_common.txt", sep = ",", header = F)

## Load information for the networks 
# Distintive colors 
IBD_dis.network <- readRDS("IBD_3Dplots_distintive.rds")
PCA_dis.network <- readRDS("PCA_3Dplots_distintive.rds")
GRM_common_dis.network <- readRDS("GRM_common_3Dplots_distintive.rds")
GRM_rare_dis.network <- readRDS("GRM_rare_3Dplots_distintive.rds")

# Similar colors
IBD_sm.network <- readRDS("IBD_3Dplots_similar.rds")
PCA_sm.network <- readRDS("PCA_3Dplots_similar.rds")
GRM_common_sm.network <- readRDS("GRM_common_3Dplots_similar.rds")
GRM_rare_sm.network <- readRDS("GRM_rare_3Dplots_similar.rds")

# By default we load the data from IBD metrics/distintive colors: 
datos_ciudades <- IBD
colors_name <- IBD_colors
resolution_plot <- IBD_resolution
network_data <- IBD_dis.network

ui <- fluidPage(
  tags$head(tags$meta(property = "og:image", content = "prev.png"),),
  title = "GG-NC",
  navbarPage("Global Genomic Network Communities (GG-NC) Browser",
    tabPanel("Map",
      sidebarLayout(
        sidebarPanel(
          selectInput("metrica", "Select a genetic metric:", 
                    choices= c("Identity By Descent", 
                                "Principal Component Analysis",
                                "Genetic Relationship Matrix based on rare variants", 
                                "Genetic Relationship Matrix based on common variants")),
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
        fluidRow(
          br(),
          column(7, offset= 5, align = "right",screenshotButton(label= "Download map", id= "mapa",filename = "GG-NC")),
          column(7, plotlyOutput("diagramaPie")),
          column(4, rglwidgetOutput("network")), 
          br(), br(), br(), br(),br(), br(), br(),
          column(8, offset=2, a(href = "https://www.sohaillab.com/contact", "Sohail-Lab")),
          column(5,align="center", 
               div(style= "font-size:9px;","Developed by María J. Palma, Yuridia Posadas, Claudia Quiroz, Brenda E. López, Anna Lewis, Tina Lasisi, Kevin A. Bird, Arslan Zaidi and Mashaal Sohail")),
          column(1, offset = 2, img(height = 50, width = 70, src = "https://conahcyt.mx/wp-content/uploads/2021/10/logo_conacyt_con_sintagma_azul_completo.svg")),
          column(1, offset =1, img(height = 40, width = 80, src = "https://covid19.ccg.unam.mx/images/ccg.png")),
          column(1, offset = 1, img(height = 40, width = 35, src = "https://i0.wp.com/www.atmosfera.unam.mx/wp-content/uploads/2019/06/unam-escudo-azul.png?ssl=1")),
          column(4, offset = 7, align="center",div(style= "font-size:10px",a(href= "https://conahcyt.mx/convocatorias/convocatoria-de-ciencia-basica-y-o-ciencia-de-frontera-modalidad-paradigmas-y-controversias-de-la-ciencia-2022/", "Funding provided by Conahcyt"))),
          column(10, offset=2,a(href= "https://creativecommons.org/licenses/by/4.0/?ref=chooser-v1",img(height = 20, width = 20, src = "https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1"))),
          column(8)
          )
        )
      )
    ), # tabpanel: map
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
    shinyalert("Help", "Welcome to the Global Genomic Network Communities (GG-NC) Browser!
    
How do we figure out who's genetically similar to who? Think of each person's genome like its own unique city. The viewpoint you choose and what you decide to focus on can completely change how you see its similarity to another city. Are you looking from a helicopter focusing on skyscrapers, or are you on the street counting parks? That's the kind of choice the Global Genomic Network Communities (GG-NC) Browser offers. You can select different metrics—like Identity by Descent (IBD), Principal Component Analysis (PCA), or Genetic Relationship Matrix (GRM)—to change your 'viewpoint' and focus on different 'landmarks' to compare two genomes. Dive in and see how the story changes when you start looking at things differently!

To begin you can pick a genetic metric and a resolution value in the left panel. So far you can choose between three metrics: Identity by Descent (IBD), Principal Component Analysis (PCA), and Genetic Relationship Matrix (GRM).  The fifty resolution values are in the logarithmic space (base 10) from -2 to 2. A network is computed using the genetic metric chosen, and communities are detected on that network at the chosen resolution using Louvain’s algorithm.

Any pie chart shows the individuals in that cohort that belong to the network communities detected at that resolution. You can play around and see how the communities change across these values. On the bottom left, you will see a graph showing the communities detected across our range of resolution values, colored in the same way as the pie charts. The individuals in this graph are ordered according to the communities detected, and can also be observed on the world map. As a note, any community with less than 6 members is colored in white.

The available cohorts come from the 1000 Genomes Project and the Human Genome Diversity Project. To display a cohort’s name, hover over a pie chart on the map.  If you wish to take a closer look into a particular cohort, select its name in the left panel. You can also zoom in on the map for a closer view of a particular geographic area.

This is still a project in progress and modifications will be included in the near future in new versions. Questions and comments are very welcome!


Contact: 
mashaal at ccg dot unam dot mx
mpalma at ccg dot unam dot mx\n
If you use the maps in other works, include the source link:
https://sohail-lab.shinyapps.io/GG-NC/ \n
We plan to have a pre-print posted shortly.\n
Version\n0.2 (beta) 
", 
               type = "info")})
  
  output$value_text <- renderText({
    selected_value <- input$slider_id
    lambda <- sort(unique(datos_ciudades$Lambda))
    log_lambda <- lambda[selected_value]
    paste("Resolution:", round(log10(log_lambda),3))
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
  
  output$mapa <- renderLeaflet({
    # Select the data, the colors and the resolution plot according to the metrics selected by the user
    # For IBD
    if (input$mySwitch) {
      colors_name <- IBD_similar_colors
    }
    
    if(input$metrica=="Principal Component Analysis"){
      datos_ciudades <- PCA
      resolution_plot <- PCA_resolution
      if (input$mySwitch) {
      colors_name <- PCA_similar_colors
      } else {
      colors_name <- PCA_colors
      }
    }
    
    if(input$metrica=="Genetic Relationship Matrix based on rare variants"){
      datos_ciudades <- GRM_rare
      resolution_plot <- GRM_rare_resolution
      if (input$mySwitch) {
        colors_name <- GRM_similar_rare_colors
      } else {
      colors_name <- GRM_rare_colors
      }
    }
    
    if(input$metrica=="Genetic Relationship Matrix based on common variants"){
      datos_ciudades <- GRM_common
      resolution_plot <- GRM_common_resolution
      if (input$mySwitch) {
        colors_name <- GRM_similar_common_colors
      } else {
      colors_name <- GRM_common_colors
      }
    }
    
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
      yticks <- seq(1, nrow(resolution_plot), length.out = 5)
      
      axis(
        side = 2,
        at = yticks,
        labels = yticklabels,
        las = 1
      )
      
      # Lambda value 
      selected_value <- (50:1)[input$slider_id]
      abline(h = selected_value, col = "black", lwd=3, lty=2)
    })
    
    # Lambda value 
    selected_value <- input$slider_id
    lambda <- sort(unique(datos_ciudades$Lambda))
    log_lambda <- lambda[selected_value]
    
    # Data for the lambda selected 
    data <- datos_ciudades[datos_ciudades$Lambda == log_lambda,]
    
    # Get the number of columns
    lcolumn <- length(datos_ciudades)
    
    # Obtaining the possition of the first community
    c1 <- which(colnames(datos_ciudades)=="C1")
    
    # Printing number of communities 
    communities.info <- data[,c1:lcolumn]
    number.com <- length(colnames(communities.info)[apply(communities.info, 2, function(col) any(!is.na(col)))])
    # Delete white comunnity 
    number.com <- number.com 
    
    output$community.num <- renderText({
      paste("Total of detected comunities:", number.com)
    })
    
    leaflet() %>%
      addTiles() %>%
      
      # Add pie charts
      addMinicharts(
        lat = data$Latitud,
        lng = data$Longitud,
        type = "pie",
        colorPalette = colors_name,
        legend = FALSE,
        popup = popupArgs(html = paste0(data$Pop)),
        chartdata = data[,c1:lcolumn])
  })
  
  # Add tooltip for the resolution plot
  addPopover(session,"resolution_plot", "Resolution Plot", 
             content = "This plot displays the identified communities at various resolution values. 
               The resolution values span a logarithmic space from -2 to 2. The y-axis represents the resolution values, 
               while the x-axis denotes the index of individuals. The colors in this plot correspond 
               to those displayed in the accompanying map. Communities that were less tan six members or were dected in only one 
             resolution are colored in white", 
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