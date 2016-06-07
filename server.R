library(shiny)
library(png)
source("functions_integrated_tracks_figure.R")
na_pos = "NA:NA-NA"
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
#options(shiny.reactlog = T)
#options(shiny.error = browser)


shinyServer(function(input, output, session) {
  
  output$profilePlots = renderPlot({
    req(input$cellType, input$drugTreatment)
    print("profilePlots")
    print(paste("plot", values$plot_number))
    pos = values$setPos
    print(paste(input$cellType, input$drugTreatment, pos))
    if(values$plot_number > -1){
      plot_title = paste0(values$setGeneSymbol, "\n", input$cellType, " ", input$drugTreatment)
      figure_track_plots(cell = input$cellType, drug = input$drugTreatment, ucsc_rng = pos, add_ref_img = T, plot_title = plot_title, interpret_states = input$interpretStates)
    }else{
      plot(0:1, 0:1); text(.5, .5, "waiting for input")
    }
  })
  output$plotUI = renderUI({
    print(paste("ui", values$ui_number))
    values$plot_number = isolate(values$plot_number) + 1
    plotOutput("profilePlots", height = 800, width = 1000,
               brush = brushOpts(id = 'plot_brush', fill = rgb(0,0,1), stroke = "black", opacity = .2, clip = T, 
                                 direction = "x", delay = 600,
                                 delayType = 'debounce', resetOnNew = T),
               dblclick = dblclickOpts(id = "plot_dblClick"))
  })
  
  getPossible = function(gs, type, width){
    possible = NA
    if(type == "promoter"){
      possible = symbol2promoter_ucsc(gs, ext = as.integer(width))
    }else if(type == "gene body"){
      possible = symbol2ucsc(gs)
    }else{
      warning("invalid featureType set")
    }
    return(possible)
  }
  
  observeEvent(input$geneSymbol, {#when geneSymbols changes, check if valid and update values$setGeneSymbol
    print("observer input$geneSymbol")
    req(input$geneSymbol)
    gs = toupper(input$geneSymbol)
    if(length(intersect(all_symbols, gs)) < 1){
      return()#do nothing if no match
    }
    updateTextInput(session, inputId = "geneSymbol", value = gs)
    if(values$plot_number > 1){
      values$geneLists$history = union(isolate(values$geneLists$history), gs)
    }
    values$setGeneSymbol = gs
  })
  
  observeEvent(c(values$setGeneSymbol, input$featureType, input$promoterWidth), {#when values$SetgeneSymbol changes, update input$chrPos and values$chrPos
    #input$geneSymbol, {#updates chrPos if gene symbol is valid or promoter/feature changes
    print("observer values$setGeneSymbol")
    req(values$setGeneSymbol, input$featureType, input$promoterWidth)
    possible = getPossible(gs = values$setGeneSymbol,  type = input$featureType, width = input$promoterWidth)
    if(possible != na_pos && possible != values$lastPos){
      print(paste('updating possition with', possible))
      #       values$setPos = possible
      updateTextInput(session, inputId = "chrPos", value = possible)
    }
  })
  
  observeEvent(input$chrPos, { #check chrPos updates plotNumber each time a new plot should be drawn
    req(input$chrPos)
    print("observe chrPos")
    
    input_pos = input$chrPos #should not call when chrPos is being set here or elsewhere
    
    if(input_pos != na_pos){
      if(grepl(",", input_pos)){
        print("cleanup chrPos")
        input_pos = gsub(",", "", input_pos)
        updateTextInput(session, inputId = "chrPos", value = input_pos)  
      }
      
      split_pos = strsplit(input_pos, "[:-]")[[1]]
      chrm = split_pos[1]; 
      start_pos = as.integer(split_pos[2]); 
      end_pos = as.integer(split_pos[3])
      if(end_pos > start_pos && values$setPos != input_pos){#Pos must be valid and new
        print("force plot update with new chrPos")
        values$lastPos = values$setPos
        values$setPos = input_pos
        #         values$ui_number = isolate(values$ui_number) + 1
      }
    }
  })
  
  values = reactiveValues(
    plot_number = 0, 
    ui_number = 0, 
    geneLists = list(history = character()),#geneLists available for export
    lastPos = na_pos,
    setPos = na_pos,#position tested for validity
    setGeneSymbol = "PGR"#symbol tested for validiting
  )
  
  
  make_filename = reactive({
    bname = paste(input$geneSymbol, sub(":", "-", input$chrPos), sep = "_")
    desc = paste(input$cellType, input$drugTreatment, input$featureType, sep = "_")
    if(input$featureType == "promoter") desc = paste(desc, input$promoterWidth, sep = "")
    bname = paste(bname, desc, sep = ".")
    fname = paste0(bname, ".pdf")
    return(fname)
  })
  
  output$dlImage = downloadHandler(filename = make_filename,
                                   content = function(file){
                                     img_width = session$clientData$output_profilePlots_width
                                     img_height = session$clientData$output_profilePlots_height
                                     pdf(file, width = img_width / 72, height = img_height / 72)
                                     pos = getPossible(gs = input$geneSymbol, type = input$featureType, width = input$promoterWidth)
                                     plot_title = paste0(input$geneSymbol, "\n", input$cellType, " ", input$drugTreatment)
                                     figure_track_plots(cell = input$cellType, drug = input$drugTreatment, ucsc_rng = input$chrPos, add_ref_img = T, plot_title = plot_title)
                                     dev.off()
                                   })
  
  output$dlGeneLists = downloadHandler(
    
    filename = ifelse(length(input$selectedGeneLists) == 1, 
                      paste0(input$selectedGeneLists, ".zip"),
                      "batch.zip"),
    content = function(file){
      print("dlGeneLists")
      out_dir = paste("tmp", as.integer(Sys.time()), sep = "_")
      dir.create(out_dir)
      img_width = session$clientData$output_profilePlots_width
      img_height = session$clientData$output_profilePlots_height
      
      for(gl_name in input$selectedGeneLists){
        bname = strsplit(gl_name, "\\.")[[1]][1]
        desc = paste(input$cellType, input$drugTreatment, input$featureType, sep = "_")
        if(input$featureType == "promoter") desc = paste(desc, input$promoterWidth, sep = "")
        bname = paste(bname, desc, sep = ".")
        fname = paste0(out_dir, "/" , bname, ".pdf")
        pdf(fname, width = img_width / 72, height = img_height / 72)
        glist = values$geneLists[[gl_name]]
        glist = intersect(glist, all_symbols)
        for(gs in glist){
          print(gs)
          pos = getPossible(gs = gs, type = input$featureType, width = input$promoterWidth)
          plot_title = paste0(gs, "\n", input$cellType, " ", input$drugTreatment)
          figure_track_plots(cell = input$cellType, drug = input$drugTreatment, ucsc_rng = pos, add_ref_img = T, plot_title = plot_title)
        }
        dev.off()
      }
      zipfiles = dir(out_dir, full.names = T)
      if(length(zipfiles) == 0){
        write("", file = "no_files_to_download")  
        zipfiles = "no_files_to_download"
      }
      print(zipfiles)
      zip(zipfile = file, files = zipfiles, flags = "-j")
    })
  
  observeEvent(input$upGenes, {#parse and record uploaded genes
    print(input$upGenes)
    genes = read.table(input$upGenes$datapath, stringsAsFactors = F, header = F, quote = "")[,1]
    values$geneLists[[input$upGenes$name]] = toupper(genes)
  })
  
  output$availableGeneLists = renderUI({
    selectInput(inputId = "selectedGeneLists", label = "Gene Lists to Export", choices = names(values$geneLists), multiple = T)
  })
  
  #   source("server_plot_reactivity.R")
  observeEvent(input$plot_dblClick, {
    print("dbl_click")
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
