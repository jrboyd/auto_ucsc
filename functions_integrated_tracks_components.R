
# already_seen = list(); names(already_seen) = character()

all_bigwigs = dir("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/drug_treated/hg38//pooled", full.names = T)
mcf7_rna_bigwigs = dir("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/hg38/RNA_seq/Breast_MCF7_Pfizer", full.names = T)
mcf10a_rna_bigwigs = dir("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/hg38/RNA_seq/Breast_MCF10A_Pfizer", full.names = T)
rna_bigwigs = c(mcf10a_rna_bigwigs, mcf7_rna_bigwigs)
drugs_list = list(bza = character(), gc10 = character(), gc10bza = character(), e2 = character(), e2bza = character(), ctrl = character())
cells_list = list(MCF10A = drugs_list, MCF7 = drugs_list)
rna_bw = list(positive = cells_list, negative = cells_list)
for(strand in names(rna_bw)){
  for(cell in names(rna_bw[[strand]])){
    for(drug in names(rna_bw[[strand]][[cell]])){
      d = paste0("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/hg38/RNA_seq/Breast_", cell, "_Pfizer")
      files = dir(d, pattern = strand, full.names = T)
      key = drug
      if(key == "ctrl") key = "Control"
      key = paste0(cell, "_", key, "_", "R")
      keep = grepl(key, files, ignore.case = T)
      if(sum(keep) == 0){
        key = sub("e2bza", "bza_e2", key)
        key = sub("gc10bza", "bza_gc10", key)
        keep = grepl(key, files, ignore.case = T)
      }
      rna_bw[[strand]][[cell]][[drug]] = files[keep]
    }
  }
}

# lapply(rna_bw, function(x){
#   sapply(x, function(y){
#     sapply(y, length)
#   })
# })

drugs = sapply(strsplit(basename(rna_bigwigs), "_"), function(x)x[2])
unique(drugs)
bigwigs = all_bigwigs[0:5 * 4 + 1]
bigwigs = all_bigwigs[1:4]
mark_colors = brewer.pal(5, "Set1")[c(1:3,5)]
names(mark_colors) = c("H3K27ME3", "H3K4AC", "H3K4ME3", "H3K27AC")

sum_bw_files = function(bw_files, chrm, start, end){
  bw_vals = lapply(bw_files, function(x){
    out_mat = get_bw_matrix(x, chrm, start, end)
    starts = as.integer(out_mat[,2])
    ends = as.integer(out_mat[,3])
    vals = as.integer(out_mat[,4])
    lens = ends - starts
    rle = Rle(values = c(0, vals), lengths = c(min(starts), lens))
    return(rle)
  })
  sum_prof = bw_vals[[1]]
  hidden = lapply(bw_vals[2:length(bw_vals)], function(x){
    suppressWarnings(sum_prof <<- sum_prof + x)
  })
  return(sum_prof)
}

get_plot_coords_from_Rle = function(prof_sum){
  ends = cumsum(runLength(prof_sum))
  starts = c(0, ends[-length(ends)])
  vals = runValue(prof_sum)
  xs = sort(c(starts,ends))
  ys = vals[sort(rep(1:length(vals),2))]
  return(list(xs = xs, ys = ys))
}

fetch_rna_profile = function(cell = c("MCF10A", "MCF7"), drug = c("ctrl", "bza", "gc10bza", "gc10", "e2", "e2bza"),
                             chrm, start, end){
  drug = drug[1]
  cell = cell[1]
  pos_files = rna_bw$positive[[cell]][[drug]]
  neg_files = rna_bw$negative[[cell]][[drug]]
  pos_sum = sum_bw_files(pos_files, chrm, start, end)
  neg_sum = sum_bw_files(neg_files, chrm, start, end)
  total_sum = suppressWarnings(pos_sum + neg_sum)
  
  pos_coords = get_plot_coords_from_Rle(pos_sum)
  neg_coords = get_plot_coords_from_Rle(neg_sum)
  total_coords = get_plot_coords_from_Rle(total_sum)
  out = list(pos_coords = pos_coords,
             neg_coords = neg_coords,
             total_coords = total_coords)
  return(out)
}

get_bw_matrix = function(bw, chrm, start, end){
  arg_template = "BIGWIG -chrom=CHRM -start=START -end=END stdout" #ALLCAPS will be replaced by args
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
  if(length(out_split) > 0){
    nr = length(out_split)
    nc = length(out_split[[1]])
    out_mat = matrix(unlist(out_split), nrow = nr, ncol = nc, byrow = T)
  }else{
    out_mat = matrix("", nrow = 1, ncol = 4)
    out_mat[,1] = chrm
    out_mat[,2] = start
    out_mat[,3] = end
    out_mat[,4] = 0
  }
  return(out_mat)
}

plot_prof = function(xs, ys, start, end, ylab = "",
                     log_transform = F,
                     prof_color = "black", new_plot = F,  
                     add_axis_labels = F, add_y_scale = T, 
                     ymin = 0, ymax = 1, plot_y = 0, plot_h = 1){
  if(!new_plot){add_y_scale = F; add_axis_labels = F}#when not a new plot axes should not be messed with
  if(log_transform){ys = log10(ys); ys = ifelse(ys < 0, 0, ys); ymax = log10(ymax)}
  if(new_plot) plot(xs,ys, type = "n", ylim = c(ymin, ymax), xlim = c(as.integer(start), as.integer(end)), bty = "L", axes = F, xlab = "", ylab = ylab, log = )
  at = as.numeric(c(start, axTicks(1), end))
  axis(side = 1, labels = add_axis_labels, at = at)
  
  if(add_y_scale){
    def_at = c(0,.5,1)
    def_lab = c(0, "", 1)
    if(ymax > 1){
      def_at = c(ymin, ymax)
      def_lab = c(ymin, ifelse(log_transform, 10^ymax, ymax))
    }
    axis(side = 2, at = def_at, labels = def_lab)
  }
  #each positive island between zeroes needs to be plotted separately
  if(max(ys) == 0) return()#case of all 0
  toplot_xs = xs
  toplot_ys = ys
  done = F
  while(max(toplot_ys) > 0){#case ends on 0
    #     print(length(toplot_xs))
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
    block_ys = c(ymin, block_ys)
    block_xs = c(block_xs, max(block_xs))
    block_ys = c(block_ys, ymin)
    block_xs = c(block_xs, min(block_xs))
    block_ys = c(block_ys, ymin)
    
    norm_ys = (block_ys / ymax * plot_h)
    norm_ys = ifelse(norm_ys > plot_h, plot_h, norm_ys)
    
    polygon(block_xs, norm_ys  + plot_y, col = prof_color, border = F)
    if(e == length(toplot_ys)){
      break
    }else{
      toplot_xs = toplot_xs[(e+1):length(toplot_xs)]
      toplot_ys = toplot_ys[(e+1):length(toplot_ys)]
    }
  }
}

plot_bw = function(bw, chrm, start, end, bw_color = "black", plot_y = 0, plot_h = .75){
  attrib = strsplit(basename(bw), "_")[[1]]
  out_mat = get_bw_matrix(bw, chrm, start, end)
  #   if(is.null(out_mat)){
  #     xs = c(start, end)
  #     ys = c(0, 0)
  #   }else{
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
  #   }
  plot_prof(xs = xs, ys = ys, start = start, end = end, prof_color = bw_color, new_plot = F, add_axis_labels = F, plot_y = plot_y, plot_h = plot_h)
}



#cell = "MCF7"; drug = "e2"; ucsc_rng = symbol2ucsc("PGR"); add_ref_img = T; mark_colors = c("H3K27ME3" = "#E41A1C", "H3K4AC" = "#377EB8","H3K4ME3" ="#4DAF4A","H3K27AC" = "#FF7F00"); add_ref_img = T; plot_title = ucsc_rng

# figure_bw_plots(cell = "MCF7", drug = "e2", ucsc_rng = symbol2ucsc("PGR"), add_ref_img = T)
