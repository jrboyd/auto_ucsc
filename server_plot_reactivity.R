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