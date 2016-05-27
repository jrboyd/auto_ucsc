library(RColorBrewer )
# already_seen = list(); names(already_seen) = character()

all_bigwigs = dir("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/drug_treated/hg38//pooled", full.names = T)
bigwigs = all_bigwigs[0:5 * 4 + 1]
mark_colors = brewer.pal(5, "Set1")[c(1:3,5)]
names(mark_colors) = c("H3K27ME3", "H3K4AC", "H3K4ME3", "H3K27AC")

ucsc_rng = "chr10:17,229,220-17,231,407"


plot_bw = function(bw, chrm, start, end, bw_color = "black"){
  attrib = strsplit(basename(bw), "_")[[1]]
  arg_this = arg_template
  arg_this = sub("BIGWIG", bw, arg_this)
  arg_this = sub("CHRM", chrm, arg_this)
  arg_this = sub("START", start, arg_this)
  arg_this = sub("END", end, arg_this)
  
  out_raw = system2(command = "bigWigToWig", stderr = F, stdout = T,
                    args =  arg_this) #query bigwig range
  keep = !grepl("^#", out_raw) #find and remove commented lines
  out_raw = out_raw[keep]
  out_split = strsplit(out_raw[-1], "\t")
  nr = length(out_split)
  nc = length(out_split[[1]])
  out_mat = matrix(unlist(out_split), nrow = nr, ncol = nc, byrow = T)
  
  starts = as.integer(out_mat[,2])
  ends = as.integer(out_mat[,3])
  vals = as.numeric(out_mat[,4])
  #zipper starts and ends and double in place vals
  xs = sort(c(starts,ends))
  ys = vals[sort(rep(1:length(vals),2))]
  
  MIN = 0
  MAX = 1
  ys = ifelse(is.infinite(ys), MIN, ys)
  ys = ifelse(ys < MIN, MIN, ys)
  
  
  plot(xs,ys, type = "n", ylim = c(MIN,MAX), xlim = c(as.integer(start), as.integer(end)), bty = "L", axes = F)
  at = as.numeric(c(start, axTicks(1), end))
  lab = F
  if(bw == bigwigs[length(bigwigs)]) lab = T
  axis(side = 1, labels = lab, at = at)
  #   lines(rep(as.integer(start), 2), c(MIN,MAX))
  #   lines(rep(as.integer(end), 2), c(MIN,MAX))
  #each positive island between zeroes needs to be plotted separately
  if(max(ys) == 0) next#case of all 0
  
  
  
  toplot_xs = xs
  toplot_ys = ys
  done = F
  while(max(toplot_ys) > 0){#case ends on 0
    print(length(toplot_xs))
    s = (1:length(toplot_ys))[toplot_ys > 0][1]
    toplot_xs = toplot_xs[s:length(toplot_xs)]
    toplot_ys = toplot_ys[s:length(toplot_ys)]
    if(min(toplot_ys) > 0){
      e = length(toplot_ys)#case ends above 0
      
    }else{
      e = (1:length(toplot_ys))[toplot_ys == 0][1] - 1#should be most common
    }
    
    block_xs = toplot_xs[1:e]
    block_ys = toplot_ys[1:e]
    block_xs = c(min(block_xs), block_xs)#these close the polygon
    block_ys = c(MIN, block_ys)
    block_xs = c(block_xs, max(block_xs))
    block_ys = c(block_ys, MIN)
    block_xs = c(block_xs, min(block_xs))
    block_ys = c(block_ys, MIN)
    
    polygon(block_xs, block_ys, col = bw_color, border = F)
    if(e == length(toplot_ys)){
      break
    }else{
      toplot_xs = toplot_xs[(e+1):length(toplot_xs)]
      toplot_ys = toplot_ys[(e+1):length(toplot_ys)]
    }
  }
}


figure_bw_plots = function(bigwigs, 
                     mark_colors = c("H3K27ME3" = "#E41A1C", 
                                     "H3K4AC" = "#377EB8",   
                                     "H3K4ME3" ="#4DAF4A",  
                                     "H3K27AC" = "#FF7F00"), 
                     ucsc_rng = "chr10:17,229,220-17,231,407"){
  ucsc_rng = gsub(",", "", ucsc_rng)
  
  arg_template = "BIGWIG -chrom=CHRM -start=START -end=END stdout" #ALLCAPS will be replaced by args
  
  chrm = strsplit(ucsc_rng, ":")[[1]][1]
  start = strsplit(ucsc_rng, "[:-]")[[1]][2]
  end = strsplit(ucsc_rng, "[:-]")[[1]][3]
  #need chrm start and end defined before proceeding
  
  layout(c(0:length(bigwigs),0))
  par(mai = c(0,0,0,0))
  for(bw in bigwigs){
    #   if(!any(grepl(ucsc_rng, names(already_seen)))){
    #     newList = list()
    #     names(newList) = character()
    #     already_seen[[ucsc_rng]] = newList
    #   }
    #   if(!any(grepl(ucsc_rng, names(already_seen)))){
    #     newList = list()
    #     names(newList) = character()
    #     already_seen[[ucsc_rng]] = newList
    #   }
    attrib = strsplit(basename(bw), "_")[[1]]
    cell = attrib[1]
    drug = attrib[2]
    mark = attrib[3]
    
    col_this = mark_colors[mark]
    
    plot_bw(bw = bw, chrm = chrm, start = start, end = end, bw_color = col_this)
  }
}
