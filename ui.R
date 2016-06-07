
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(title = "Integrated Genome View",
                  
                  # Application title
                  titlePanel("Tracks for Pfizer Drug Treatments"),
                  #   tags$a("Help", href = "http://45.56.125.191/help.html", target = "_blank"),
                  # Sidebar with a slider input for number of bins
                  sidebarLayout(
                    sidebarPanel(
                      tags$head(
                        tags$style(type="text/css", "select { max-width: 140px; }"),
                        tags$style(type="text/css", ".span4 { max-width: 210px; }"),
                        tags$style(type="text/css", ".well { max-width: 280px; }")
                      ),
                      h3("Location (hg38)"),
                      textInput("geneSymbol", "updates only when valid gene symbol entered", value = "PGR"),
                      radioButtons("featureType", label = "", choices = c("promoter", "gene body"), selected = "promoter"),
                      conditionalPanel(
                        condition = "input.featureType == 'promoter'",
                        selectInput("promoterWidth", label = "Promoter Width", choices = c(500, 1000, 2000, 5000, 10000, 25000, 50000, 100000), selected = 5000)
                      ),
                      #       actionButton(inputId = "gotoSymbol", label = "Go!"),
                      textInput("chrPos", "chr:start-end may be entered manually", value = ""),
                      h3("States"),
                      checkboxInput(inputId = "interpretStates", label = "Interpret States", value = T),
                      h3("Cell Line"),
                      radioButtons("cellType", label = "", choices = c("MCF10A", "MCF7"), selected = "MCF7"),
                      h3("Drug Treatments"),
                      radioButtons("drugTreatment", label = "", choices = c("ctrl", "e2", "bza", "e2bza", "gc10", "gc10bza"), selected = "e2"),
                      
                      #       h3("Shortcuts"),
                      #       actionButton("gotoRunx1", "Runx1"),
                      #       actionButton("gotoRunx2", "Runx2"),
                      #       actionButton("gotoRunx3", "Runx3"),
                      h3("Navigation"),
                      fluidRow(actionButton("zoomOut", "Zoom Out"),
                               actionButton("zoomIn", "Zoom In")),
                      fluidRow(actionButton("shiftLeft", "Shift Left"),
                               actionButton("shiftRight", "Shift Right")),
                      
                      h3("Export"),
                      downloadButton("dlImage", label = "Download Current View"), 
                      fileInput(inputId = "upGenes", label = "Upload Genes", accept = c(".txt")),
                      uiOutput(outputId = "availableGeneLists"),
                      downloadButton("dlGeneLists", label = "Download Selected Gene Lists")
                    ),
                    
                    # Show a plot of the generated distribution
                    mainPanel(
                      uiOutput(outputId = "plotUI")
                      #       plotOutput("profilePlots", height = 800, width = 1000,
                      #                  brush = brushOpts(id = 'plot_brush', fill = rgb(0,0,1), stroke = "black", opacity = .2, clip = T, 
                      #                                    direction = "x", delay = 600,
                      #                                    delayType = 'debounce', resetOnNew = T),
                      #                  dblclick = dblclickOpts(id = "plot_dblClick"))
                    )
                  )
)
)

