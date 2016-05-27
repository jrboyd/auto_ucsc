library(RColorBrewer)
library(png)
files = dir("state_ngsplots_FE/", pattern = "avgprof.RData", recursive = T, full.names = T)
keep = grepl("posState", files)
files = files[keep]
o = order(as.integer(sapply(strsplit(files, "[\\./_]"), function(x)x[10])))
files = files[o]

plot_FE_states = function(f){
  
  mark_colors = c("black", brewer.pal(5, "Dark2"))#[c(1,1,2,3)]
  
  name = basename(dirname(f))
  name = gsub("_", " ", name)
  load(f)
  dat = regcovMat
  y_rng = range(regcovMat)
  bp_rng = -2100:2100
  x_rng = range(bp_rng)
  # layout(rbind(1:2), widths = c(3,1))
  par(mai = rep(.5,4))
  plot(0, type = "n", xlim = x_rng, ylim = y_rng, axes = F, xlab = "", ylab = "")
  title(name)
  
  xs = c(-2100 + 2000 * 0:39/40,
         -100 + 200 * 0:19/20,
         100 + 2000 * 0:40/40)
  for(i in 1:ncol(dat)){
    lines(xs, dat[,i], col = mark_colors[i], lwd = 3)
  }
  lines(rep(-100, 2), y_rng)
  lines(rep(100, 2), y_rng)
  axis(1, at = c(-2100, 0, 2100), c("-2000", "200 bp window", "2000"))
  axis(2)
  lines(x_rng, rep(0,2), lty = 2, col = rgb(.3,.3,.3))
#   plot(0, type = "n", xlab = "", ylab = "", axes = F)
  if(grepl("H3K4ME3", f)) legend("topright", horiz = F, fill = mark_colors, legend = colnames(dat), bty = "n", xpd = NA)
}
# pdf("MCF7_ctrl_states.pdf")
for(f in files){
  name = paste0(basename(dirname(f)), ".png")
  # name = gsub("_", " ", name)
  png(name, res = 150, units = "in", width = 5, height = 5)
  plot_FE_states(f)
  dev.off()
}
# dev.off()

plot0 = function(width = 1, height = 1){
  fudge = 0.037037037037
  plot(c(0+fudge*width, width-fudge * width), c(0+fudge*height, height-fudge * height), type = 'n', xlab = '', ylab = '', axes = F)
}

pdf("all_posStates.pdf", width = 8, height = 38)
layout(1)
par(mai = rep(0,4))
plot0(width = 8, height = 38)
for(i in 1:length(files)){
  name = paste0(basename(dirname(files[i])), ".png")
  x = (i - 1) %% 4 * 2
  y = 36 - floor( (i - 1) / 4) * 2
  # print(paste(x, y))
  rasterImage(readPNG(name, native = FALSE), xleft = x, xright = x + 2, ytop = y + 2, ybottom = y, interpolate = FALSE)  
}
dev.off()

files = dir("state_ngsplots_FE/", pattern = "hm1.txt", recursive = T, full.names = T)
keep = grepl("posState", files)
files = files[keep]
keep = grepl("H3K4ME3", files)
files = files[keep]
o = order(as.integer(sapply(strsplit(files, "[\\./_]"), function(x)x[10])))
files = files[o]
for(f in files){
  name = basename(dirname(f))
  name = gsub("_", " ", name)
  name = paste(strsplit(name, " ")[[1]][5:6], collapse = " ")
  nr = nrow(read.table(f))
  print(paste(name, "count", nr))
}

files = dir("state_ngsplots_FE/", pattern = "hm1.txt", recursive = T, full.names = T)
keep = !grepl("posState", files)
files = files[keep]
# keep = grepl("H3K4ME3", files)
# files = files[keep]
o = order(as.integer(sapply(strsplit(files, "[\\./_]"), function(x)x[9])))
files = files[o]
for(f in files){
  name = basename(dirname(f))
  name = gsub("_", " ", name)
  name = paste(strsplit(name, " ")[[1]][4:5], collapse = " ")
  nr = nrow(read.table(f))
  print(paste(name, "count", nr))
}
