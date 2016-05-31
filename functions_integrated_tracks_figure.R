library(RColorBrewer )
library(GenomicRanges)
source("functions_state_plots.R")
source("functions_ref_from_UCSC.R")
source("fetch_states.R")
source("functions_integrated_tracks_components.R")

all_symbols = toupper(all_symbols)


figure_track_plots = function(cell = "MCF10A", drug = "e2", marks = names(mark_colors),
                              #                            chip_bigwigs, 
                              mark_colors = c("H3K27ME3" = "#E41A1C", 
                                              "H3K4AC" = "#377EB8",   
                                              "H3K4ME3" ="#4DAF4A",  
                                              "H3K27AC" = "#FF7F00"), 
                              ucsc_rng = "chr10:17,219,220-17,241,407",
                              add_ref_img = T,
                              plot_title = paste(cell, drug)){
  #cell = "MCF10A"; drug = "e2"; marks = names(mark_colors); mark_colors = c("H3K27ME3" = "#E41A1C", "H3K4AC" = "#377EB8", "H3K4ME3" ="#4DAF4A", "H3K27AC" = "#FF7F00"); ucsc_rng = "chr10:17,219,220-17,241,407"; add_ref_img = T; plot_title = paste(cell, drug)
  bot_mai = 2
  left_mai = 3
  top_mai = 0
  right_mai = .5
  cell_drug_mark_keys = paste(cell, drug, marks, sep = "_")
  cell_drug_keys = paste(cell, drug, sep = "_")
  chip_bigwigs = sapply(cell_drug_mark_keys, function(x){
    all_bigwigs[grepl(x, all_bigwigs)]
  })
  ucsc_rng = gsub(",", "", ucsc_rng)
  
  nplots = length(cell) * length(drug) * (length(marks) + 2 + 1) + add_ref_img
  
  chrm = strsplit(ucsc_rng, ":")[[1]][1]
  start = as.integer(strsplit(ucsc_rng, "[:-]")[[1]][2])
  end = as.integer(strsplit(ucsc_rng, "[:-]")[[1]][3])
  if(add_ref_img){
    ref_img = fetch_ref_image(chrm, start, end)
  }else{
    ref_img = NA
  }
  if(!is.na(ref_img[1])){#need to calculate ref_img height based on ref_img
    plot_Win = dev.size()[1] - left_mai - right_mai #plot w in screen inches
    plot_Hin = dev.size()[2] - top_mai - bot_mai #plot h in screen inches
    img_Wpx = ncol(ref_img) - 117 #image dimension in px
    img_Hpx = nrow(ref_img)
    img_Hin = img_Hpx * plot_Win / img_Wpx #convert images px to Hin it should take
    img_Hplot = (img_Hin * nplots) / (plot_Hin - img_Hin)
  }
  par(mai = c(bot_mai, left_mai, top_mai, right_mai))
  par(xpd = NA)
  plot_h = nplots - 1 + img_Hplot
  plot(c(start, end), c(0, plot_h), type = "n", xlab = "", ylab = "", axes = F) #left side of plot is legend
  par(usr = c(start, end, 0, plot_h))
  at = axisTicks(c(start, end), log = F)
  for(a in at){
    lines(rep(a, 2), c(0, plot_h - img_Hplot), lty = 2)
  }
  text(x = start - (end - start)*.25, y = plot_h * 3.5 / 5, labels = plot_title, adj = 0)
  legend(x = start - (end - start)*.3, y = plot_h * 3 / 5,  legend = c("RNA+", "RNA-", names(mark_colors)), fill = c("gray", "black", mark_colors), bty = "n", xpd = NA)
  #   lines(rep(par("usr")[1], 2), par("usr")[3:4], lty = 2)
  #   lines(rep(par("usr")[2], 2), par("usr")[3:4], lty = 2)
  #plot ref_img if supplied
  if(!is.na(ref_img[1])){
    #     par(mai = c(0,.8,0,0), xpd = T)
    plot_ref(my_image = ref_img, xmin = start, xmax = end, ymin = nplots - 1, ymax = nplots - 1 + img_Hplot)
  }
  #end plotting ref_img
  par(xpd = F)
  #plot rna profiles
  plot_heights = .7
  rna_profs = fetch_rna_profile(cell = cell, drug = drug, chrm = chrm, start = start, end = end)
  YMAX = quantile(rna_profs$total_coords$ys, .95)
  axis(side = 2, at = c(6, 6 + plot_heights), labels = c(0, round(YMAX)))
  plot_prof(xs = rna_profs$total_coords$xs,
            ys = rna_profs$total_coords$ys, 
            log_transform = T,
            start = start, end = end, 
            prof_color = "black", 
            new_plot = F, 
            plot_y = 6, plot_h = plot_heights,
            ymax = YMAX)
  plot_prof(xs = rna_profs$pos_coords$xs, 
            ys = rna_profs$pos_coords$ys, 
            log_transform = T,
            start = start, end = end, 
            prof_color = "gray", 
            new_plot = F,
            plot_y = 6, plot_h = plot_heights,
            ymax = YMAX)
  #plot chip profiles
  i = length(chip_bigwigs) + 1
  for(name in names(chip_bigwigs)){
    bw = chip_bigwigs[name]
    mark = strsplit(name, "_")[[1]][3]
    col_this = mark_colors[mark]
    plot_bw(bw = bw, chrm = chrm, start = start, end = end, bw_color = col_this, plot_y = i, plot_h = plot_heights)
    axis(side = 2, at = c(i, i + plot_heights), labels = c(0,1))
    i = i - 1
  }
  #plot states
  k = states.IDEAS$chrm == chrm
  kept.states.IDEAS = states.IDEAS[k,]
  k = kept.states.IDEAS$end > start & kept.states.IDEAS$start < end
  kept.states.IDEAS = kept.states.IDEAS[k,]
  axis(side = 2, at = c(1.5), labels = "IDEAS states", las = 2)
  plot_state_numbers(starts = kept.states.IDEAS$start, ends = kept.states.IDEAS$end, states = kept.states.IDEAS[[cell_drug_keys]], xlim = c(start, end), ylim = c(1, 2))
  
  k = states.cHMM$chrm == chrm
  kept.states.cHMM = states.cHMM[k,]
  k = kept.states.cHMM$end > start & kept.states.cHMM$start < end
  kept.states.cHMM = kept.states.cHMM[k,]
  axis(side = 2, at = c(0.5), labels = "chromHMM states", las = 2)
  plot_state_numbers(starts = kept.states.cHMM$start, ends = kept.states.cHMM$end, states = kept.states.cHMM[[cell_drug_keys]], xlim = c(start, end), ylim = c(0, 1))
  
  axis(side = 1, at = at, labels = at)
  axis(side = 1, at = c(start, end), labels = paste(c(start, end), "    "), las = 2)
}