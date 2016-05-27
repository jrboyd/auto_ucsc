if(!exists("symbol2ucsc")) source("fetch_states.R")

#load IDEAS state data
if(!exists("states.IDEAS")){
  if(!file.exists("states.IDEAS.save")){
    states.IDEAS = read.table(file = "data_IDEAS/mcf10a_and_mcf7_ideas_full.state", stringsAsFactors = F)
    states.IDEAS[,5:ncol(states.IDEAS)] = states.IDEAS[,5:ncol(states.IDEAS)] + 1
    states.IDEAS[,4] = states.IDEAS[,4] + 200
    drugs = c("bza", "ctrl", "e2bza", "e2", "gc10bza", "gc10")
    cnames = sapply(c("MCF10A", "MCF7"), function(x)paste(x, drugs, sep = "_"))[1:12]
    colnames(states.IDEAS) = c("row", "chrm", "start", "end", cnames)
    save(states.IDEAS, file = "states.IDEAS.save")
  }else{
    load("states.IDEAS.save")
  }
}

#load chromHMM state data
if(!exists("states.cHMM")){
  if(!file.exists("states.cHMM.save")){
    states.cHMM = read.table(file = "data_chromHMM//chromHMM_states_combined.body", stringsAsFactors = F)
    cnames = read.table(file = "data_chromHMM//chromHMM_states_combined.header", stringsAsFactors = F)
    cnames = sub("_30_999_segments.expanded.bed", "", cnames)
    cnames = sub("MCF10A", "MCF10A_", cnames)
    cnames = sub("MCF7", "MCF7_", cnames)
    colnames(states.cHMM) = cnames
    save(states.cHMM, file = "states.cHMM.save")
  }else{
    load("states.cHMM.save")
  }
}

# symbol2range("RUNX2")
# local =  fetch_IDEAS_by_symbol("RUNX2")
# 
# starts = local[,3]
# ends = local[,4]
# states = local[,5]

plot_state_numbers = function(starts, ends, states, xlim, ylim, hide_numbers = T){
  uniq_states = 1:35
  state_colors = rainbow(length(uniq_states))
  names(state_colors) = uniq_states
  ranges = matrix(0, ncol = 3, nrow = 0)
  i = 1
  len = length(states)
  while(i < len){
    curr_state = states[i]
    s = starts[i]
    e = ends[i]
    i = i + 1
    next_state = states[i]
    while(curr_state == next_state){
      i = i + 1
      if(i > len) break
      next_state = states[i]
    }
    e = ends[i - 1]
    ranges = rbind(ranges, c(start = s, end = e, state = curr_state))
  }
  #   plot(0, xlim = xlim, ylim = c(-1,1), type = "n", xlab = "", ylab = "", bty = "n", axes = F)
  at = c(xlim[1], axisTicks(xlim, log = F), xlim[2])
  axis(side = 1, at = at, labels = ifelse(rep(hide_numbers, length(at)), rep("", length(at)), at))
  #   axis(side = 1, labels = !hide_numbers, outer = 
  yrange = ylim[2] - ylim[1]
  apply(ranges, 1, function(rng){
    width = rng[2] - rng[1]
    rect(rng[1], ylim[1] + yrange * 2 / 3, rng[2], ylim[2] - yrange * 2 / 3, col = state_colors[rng[3]])
    #     lines(rng[1:2], rep(0,2), lwd = 3, col = state_colors[rng[3]])
    x = mean(c(rng[1], rng[2]))
    text(x = x, y = ylim[2] - yrange * 1 / 3 + .1, labels = rng[3], adj = c(.5,0))
  })
}