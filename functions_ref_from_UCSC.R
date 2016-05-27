library(RCurl)
library(png)

#retrieve a formatted ucsc png
fetch_ref_image = function(chrm, start, end){
  base_url = "http://genome.ucsc.edu/cgi-bin/hgRenderTracks?position=CHR:START-END&hideTracks=1&refGene=pack&pubs=pack%27http://genome.ucsc.edu/cgi-bin/hgRenderTracks?&refGene=pack&pubs=pack%27&hgsid=496004633_di6d4J35W6lXIBfFnZtvxmUGqlM1"
  sym_url = base_url
  sym_url = sub("CHR", chrm, sym_url)
  sym_url = sub("START", start, sym_url)
  sym_url = sub("END", end, sym_url)
  my_image <-  readPNG(getURLContent(sym_url))
  return(my_image)
}

plot_ref = function(my_image, start, end){
  img_height = nrow(my_image)
  lab_width = 117
  img_width = ncol(my_image)
#   par(mai = rep(0, 4))
  plot(x = c(start, end), y = c(0,1), axes = F, type = "n", xlab = "", ylab = "")
#   axis(1, at = c(start, end))
  w = end - start
  xleft = start -  w / (800 - 117) * 117
  rasterImage(my_image, xleft = xleft, ybottom = 0, xright = end, ytop = 1)
  return(img_height)
}
