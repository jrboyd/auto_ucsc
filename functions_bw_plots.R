library(RColorBrewer )
library(GenomicRanges)
source("functions_state_plots.R")
source("functions_ref_from_UCSC.R")
source("fetch_states.R")
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

lapply(rna_bw, function(x){
  sapply(x, function(y){
    sapply(y, length)
  })
})

drugs = sapply(strsplit(basename(rna_bigwigs), "_"), function(x)x[2])
unique(drugs)
bigwigs = all_bigwigs[0:5 * 4 + 1]
bigwigs = all_bigwigs[1:4]
mark_colors = brewer.pal(5, "Set1")[c(1:3,5)]
names(mark_colors) = c("H3K27ME3", "H3K4AC", "H3K4ME3", "H3K27AC")

ucsc_rng = "chr10:17,229,220-17,231,407"

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
                     prof_color = "black", new_plot = T,  
                     add_axis_labels = F, add_y_scale = T, 
                     ymin = 0, ymax = 1){
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
    
    polygon(block_xs, block_ys, col = prof_color, border = F)
    if(e == length(toplot_ys)){
      break
    }else{
      toplot_xs = toplot_xs[(e+1):length(toplot_xs)]
      toplot_ys = toplot_ys[(e+1):length(toplot_ys)]
    }
  }
}

plot_bw = function(bw, chrm, start, end, bw_color = "black"){
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
  plot_prof(xs = xs, ys = ys, start = start, end = end, prof_color = bw_color, new_plot = T, add_axis_labels = F)
}


figure_bw_plots = function(cell = "MCF10A", drug = "e2", 
                           #                            bigwigs, 
                           mark_colors = c("H3K27ME3" = "#E41A1C", 
                                           "H3K4AC" = "#377EB8",   
                                           "H3K4ME3" ="#4DAF4A",  
                                           "H3K27AC" = "#FF7F00"), 
                           ucsc_rng = "chr10:17,229,220-17,231,407",
                           add_ref_img = F,
                           plot_title = ucsc_rng){
  mark_drug = paste(cell, drug, sep = "_")
  bigwigs = all_bigwigs[grepl(paste0(mark_drug, "_"), all_bigwigs)]
  ucsc_rng = gsub(",", "", ucsc_rng)
  
  
  
  chrm = strsplit(ucsc_rng, ":")[[1]][1]
  start = as.integer(strsplit(ucsc_rng, "[:-]")[[1]][2])
  end = as.integer(strsplit(ucsc_rng, "[:-]")[[1]][3])
  if(add_ref_img){
    ref_img = fetch_ref_image(chrm, start, end)
  }else{
    ref_img = NA
  }
  #   ref_img = ifelse(add_ref_img, fetch_ref_image(chrm, start, end), NA)
  #need chrm start and end defined before proceeding
  added = 4
  if(!is.na(ref_img[1])) added = 5 #add area in layout for ref_img
  l_mat = cbind(c(0, 1:(length(bigwigs) + added), 0))
  l_mat = cbind(rep(1, nrow(l_mat)), ifelse(l_mat > 0, l_mat + 1, 0))
  nr = nrow(l_mat)
  l_mat[nr - 2, 1] = l_mat[nr - 2, 2]
  l_mat[nr - 2, 2] = l_mat[nr - 2, 2] + 1
  l_mat[nr - 1, 2] = l_mat[nr - 1, 2] + 1
  l_mat[nr - 1, 1] = l_mat[nr - 1, 2]
  l_mat[nr - 1, 2] = l_mat[nr - 1, 2] + 1
  l_hei = rep(1, nrow(l_mat)) #default layout heights are all equal
  if(!is.na(ref_img[1])){#need to calculate ref_img height based on ref_img
    dev_w = dev.size()[1]
    dev_h = dev.size()[2]
    plot_w = dev_w * 5 / 6 #based on default widths = c(1,5)
    img_w = ncol(ref_img) - 117 #117 may not hold if image width changes from 800
    img_h = nrow(ref_img)
    plot_h = plot_w / img_w * img_h
    img_r = plot_h / dev_h
    other_r = (1 - img_r) / (nrow(l_mat) - 1)
    l_hei = c(rep(other_r, 2), img_r, rep(other_r, nrow(l_mat) - 3))
  }
  layout(l_mat, widths = c(1,5), heights = l_hei)
  par(cex = 1.2)
  par(mai = rep(0, 4), xpd = T)
  plot(0:1, 0:1, type = "n", xlab = "", ylab = "", axes = F) #left side of plot is legend
  text(.5, .75, paste(cell, drug))
  legend("center", legend = c("RNA+", "RNA-", names(mark_colors)), fill = c("gray", "black", mark_colors), bty = "n")
  #second plot is title
  plot(c(start, end), 0:1, type = "n", xlab = "", ylab = "", axes = F) #second plot is title and scale
  rng = end - start
  text(start + rng / 2, .75, plot_title, cex = 2)
  scale = 10^round(log10(rng / 5))
  s = start + rng / 2 - scale / 2
  e = start + rng / 2 + scale / 2
  lines(c(s, e), rep(.25, 2))
  lines(rep(s, 2), c(.15, .35))
  lines(rep(e, 2), c(.15, .35))
  text(e + rng / 50, .25, paste(scale, "bp"), adj = 0)
  #end title/scale area
  
  #plot ref_img if supplied
  if(!is.na(ref_img[1])){
    par(mai = c(0,.8,0,0), xpd = T)
    plot_ref(my_image = ref_img, start = start, end = end)
  }
  #end plotting ref_img
  #plot rna profiles
  par(mai = c(.2,.8,.2,0), xpd = T)
  rna_profs = fetch_rna_profile(cell = cell, drug = drug, chrm = chrm, start = start, end = end)
  #   plot(0, bw = bw, chrm = chrm, start = start, end = end, bw_color = col_this)
  plot_prof(xs = rna_profs$total_coords$xs,
            ys = rna_profs$total_coords$ys, 
            #             ylab = round(quantile(rna_profs$total_coords$ys, .95), 0),
            log_transform = T,
            start = start, end = end, 
            prof_color = "black", 
            new_plot = T, 
            ymax = quantile(rna_profs$total_coords$ys, .95), 
            add_y_scale = T)
  plot_prof(xs = rna_profs$pos_coords$xs, 
            ys = rna_profs$pos_coords$ys, 
            log_transform = T,
            start = start, end = end, 
            prof_color = "gray", 
            new_plot = F)
  for(mark in names(mark_colors)){
    bw = bigwigs[grepl(mark, bigwigs)]
    col_this = mark_colors[mark]
    plot_bw(bw = bw, chrm = chrm, start = start, end = end, bw_color = col_this)
  }
  par(xpd = NA)
  k = states.IDEAS$chrm == chrm
  kept.states.IDEAS = states.IDEAS[k,]
  k = kept.states.IDEAS$end > start & kept.states.IDEAS$start < end
  kept.states.IDEAS = kept.states.IDEAS[k,]
  plot(0:1, 0:1, type = "n", xlab = "", ylab = "", axes = F)
  text(.5, .5, "IDEAS states")
  plot_state_numbers(starts = kept.states.IDEAS$start, ends = kept.states.IDEAS$end, states = kept.states.IDEAS[[mark_drug]], xlim = c(start, end))
  
  k = states.cHMM$chrm == chrm
  kept.states.cHMM = states.cHMM[k,]
  k = kept.states.cHMM$end > start & kept.states.cHMM$start < end
  kept.states.cHMM = kept.states.cHMM[k,]
  plot(0:1, 0:1, type = "n", xlab = "", ylab = "", axes = F)
  text(.5, .5, "chromHMM states")
  plot_state_numbers(starts = kept.states.cHMM$start, ends = kept.states.cHMM$end, states = kept.states.cHMM[[mark_drug]], xlim = c(start, end), hide_numbers = F)
}

#cell = "MCF7"; drug = "e2"; ucsc_rng = symbol2ucsc("PGR"); add_ref_img = T; mark_colors = c("H3K27ME3" = "#E41A1C", "H3K4AC" = "#377EB8","H3K4ME3" ="#4DAF4A","H3K27AC" = "#FF7F00"); add_ref_img = T; plot_title = ucsc_rng

# figure_bw_plots(cell = "MCF7", drug = "e2", ucsc_rng = symbol2ucsc("PGR"), add_ref_img = T)
