library(png)
source("functions_integrated_tracks_figure.R")
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
#options(shiny.reactlog = T)

library(shiny)

shinyServer(function(input, output, session) {
  
  output$profilePlots = renderPlot({
    req(input$chrPos, input$cellType, input$drugTreatment)
    print("profilePlots")
    print(values$plot_number)
    pos = isolate(input$chrPos)
    print(paste(input$cellType, input$drugTreatment, pos))
    if(pos != "none"){
      figure_track_plots(cell = input$cellType, drug = input$drugTreatment, ucsc_rng = pos, add_ref_img = T)
      
    }else{
      #par(mai = c(.887232,  2.81696, 1.963, .321621)) #measured margins from figure
      
      plot(0:1, 0:1); text(.5, .5, "waiting for position")
    }
    layout(1)
    par(mai = c(.887232,  2.81696, 1.963, .321621))
  })
  output$plotUI = renderUI({
    print(values$ui_number)
    values$plot_number = isolate(values$plot_number) + 1
    plotOutput("profilePlots", height = 800, width = 1000,
               brush = brushOpts(id = 'plot_brush', fill = rgb(0,0,1), stroke = "black", opacity = .2, clip = T, 
                                 direction = "x", delay = 600,
                                 delayType = 'debounce', resetOnNew = T),
               dblclick = dblclickOpts(id = "plot_dblClick"))
  })
  
  observe({#updates chrPos if gene symbol is valid
    print("observer geneSymbol")
    req(input$geneSymbol, input$featureType)
    met = NULL
    gs = toupper(input$geneSymbol)
    if(length(intersect(all_symbols, gs)) < 1){
      return()#do nothing if no match
    }
    if(!is.null(input$featureType)){
      print("start tests")
      if(input$featureType == "promoter"){
        possible = symbol2promoter_ucsc(gs, ext = as.integer(input$promoterWidth))
      }else if(input$featureType == "gene body"){
        possible = symbol2ucsc(gs)
      }else{
        warning("invalid featureType set")
      }
      #       print(possible)
      if(possible != "NA:NA-NA"){
        updateTextInput(session, inputId = "chrPos", value = possible)
      }
    }
  })
  
  observe({#check chrPos updates plotNumber each time a new plot should be drawn
    req(input$chrPos)
    print("observe chrPos")
    
    pos = input$chrPos
    
    if(pos != "none"){
      if(grepl(",", pos)){
        pos = gsub(",", "", pos)
        updateTextInput(session, inputId = "chrPos", value = pos)  
      }
      
      pos = strsplit(pos, "[:-]")[[1]]
      chrm = pos[1]; start = as.integer(pos[2]); end = as.integer(pos[3])
      if(end > start){
        values$ui_number = isolate(values$ui_number) + 1
      }
    }
  })
  values = reactiveValues(plot_number = 0, ui_number = 0)
  make_filename = function(ext){
    paste0(input$cellType, "_", input$geneSymbol, ".", ext)
  }
  png_name = function(){
    make_filename("png")
  }
  ref_name = function(){
    make_filename("ref.png")
  }
  
  output$dlImage = downloadHandler(filename = ref_name,
                                   content = function(file){
                                     img_width = session$clientData$output_profilePlots_width
                                     img_height = session$clientData$output_profilePlots_height
                                     png(file, img_width, img_height)
                                     figure_bw_plots(cell = input$cellType, drug = input$drugTreatment, ucsc_rng = symbol2ucsc(input$geneSymbol))
                                     dev.off()
                                   }
  )
  
  observeEvent(input$plot_dblClick, {
    dbl = input$plot_dblClick
    disp_start = input$plot_brush$xmin
    disp_end = input$plot_brush$xmax
    print(c(disp_start, disp_end))
    chrm = strsplit(input$chrPos, "[:-]")[[1]][1]
    newPos = paste0(chrm, ":", as.integer(disp_start), "-", as.integer(disp_end))
    updateTextInput(session, inputId = "chrPos", value = newPos)
  })
  
  observeEvent(input$zoomOut, handlerExpr = {
    chr = strsplit(input$chrPos, "[:-]")[[1]][1]
    chr_start = as.integer(strsplit(input$chrPos, "[:-]")[[1]][2])
    chr_end = as.integer(strsplit(input$chrPos, "[:-]")[[1]][3])
    len = chr_end - chr_start
    cen = chr_start + len / 2
    len = len * 10
    new_start = round(cen - len / 2)
    new_end = round(cen + len / 2)
    updateTextInput(session, inputId = "chrPos", value = 
                      paste0(chr, ":", new_start, "-", new_end))
  })
  observeEvent(input$zoomIn, handlerExpr = {
    chr = strsplit(input$chrPos, "[:-]")[[1]][1]
    chr_start = as.integer(strsplit(input$chrPos, "[:-]")[[1]][2])
    chr_end = as.integer(strsplit(input$chrPos, "[:-]")[[1]][3])
    len = chr_end - chr_start
    cen = chr_start + len / 2
    len = len / 10
    new_start = round(cen - len / 2)
    new_end = round(cen + len / 2)
    updateTextInput(session, inputId = "chrPos", value = 
                      paste0(chr, ":", new_start, "-", new_end))
  })
  observeEvent(input$shiftLeft, handlerExpr = {
    chr = strsplit(input$chrPos, "[:-]")[[1]][1]
    chr_start = as.integer(strsplit(input$chrPos, "[:-]")[[1]][2])
    chr_end = as.integer(strsplit(input$chrPos, "[:-]")[[1]][3])
    len = chr_end - chr_start
    cen = chr_start + len / 2
    new_start = round(chr_start - len / 3)
    new_end = round(chr_end - len / 3)
    updateTextInput(session, inputId = "chrPos", value = 
                      paste0(chr, ":", new_start, "-", new_end))
  })
  observeEvent(input$shiftRight, handlerExpr = {
    chr = strsplit(input$chrPos, "[:-]")[[1]][1]
    chr_start = as.integer(strsplit(input$chrPos, "[:-]")[[1]][2])
    chr_end = as.integer(strsplit(input$chrPos, "[:-]")[[1]][3])
    len = chr_end - chr_start
    cen = chr_start + len / 2
    new_start = round(chr_start + len / 3)
    new_end = round(chr_end + len / 3)
    updateTextInput(session, inputId = "chrPos", value = 
                      paste0(chr, ":", new_start, "-", new_end))
  })
  
})
