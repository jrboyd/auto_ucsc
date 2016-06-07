if(!exists("symbol2ucsc")) source("fetch_states.R")

#load IDEAS state data
if(!exists("states.IDEAS")){
  if(!file.exists("../data_auto_ucsc/states.IDEAS.save")){
    states.IDEAS = read.table(file = "../data_auto_ucsc/data_IDEAS/mcf10a_and_mcf7_ideas_full.state", stringsAsFactors = F)
    states.IDEAS[,5:ncol(states.IDEAS)] = states.IDEAS[,5:ncol(states.IDEAS)] + 1
    states.IDEAS[,4] = states.IDEAS[,4] + 200
    drugs = c("bza", "ctrl", "e2bza", "e2", "gc10bza", "gc10")
    cnames = sapply(c("MCF10A", "MCF7"), function(x)paste(x, drugs, sep = "_"))[1:12]
    colnames(states.IDEAS) = c("row", "chrm", "start", "end", cnames)
    save(states.IDEAS, file = "../data_auto_ucsc/states.IDEAS.save")
  }else{
    load("../data_auto_ucsc/states.IDEAS.save")
  }
}

#load chromHMM state data
if(!exists("states.cHMM")){
  if(!file.exists("../data_auto_ucsc/states.cHMM.save")){
    states.cHMM = read.table(file = "../data_auto_ucsc/data_chromHMM//chromHMM_states_combined.body", stringsAsFactors = F)
    cnames = read.table(file = "../data_auto_ucsc/data_chromHMM//chromHMM_states_combined.header", stringsAsFactors = F)
    cnames = sub("_30_999_segments.expanded.bed", "", cnames)
    cnames = sub("MCF10A", "MCF10A_", cnames)
    cnames = sub("MCF7", "MCF7_", cnames)
    colnames(states.cHMM) = cnames
    save(states.cHMM, file = "../data_auto_ucsc/states.cHMM.save")
  }else{
    load("../data_auto_ucsc/states.cHMM.save")
  }
}

#convert state number to interpretable descriptions
IDEAS2desc = read.table("../data_auto_ucsc/IDEAS_states.txt", sep = "\t")
IDEAS2desc = apply(IDEAS2desc, 1, function(x)ifelse(x[2] == "", x[1], x[2]))
names(IDEAS2desc) = as.character(1:length(IDEAS2desc))

chromHMM2desc = read.table("../data_auto_ucsc/chromHMM_states.txt", sep = "\t")
chromHMM2desc = apply(chromHMM2desc, 1, function(x)ifelse(x[2] == "", x[1], x[2]))
names(chromHMM2desc) = as.character(1:length(chromHMM2desc))

plot_state_descriptions = function(starts, ends, states, xlim, ylim, type = c("IDEAS", "chromHMM")[1], hide_numbers = T){
  if(type == "IDEAS"){
    plot_state_numbers(starts, ends, IDEAS2desc[states], xlim, ylim, hide_numbers)
  }else if(type == "chromHMM"){
    plot_state_numbers(starts, ends, chromHMM2desc[states], xlim, ylim, hide_numbers)
  }else{
    stop("unrecognized type")
  }
}

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
  alt = 1
  apply(ranges, 1, function(rng){
    s = as.integer(rng[1])
    e = as.integer(rng[2])
    txt = rng[3]
    width = e - s
    rect(s, ylim[1] + yrange * 2 / 3, e, ylim[2] - yrange * 2 / 3, col = state_colors[txt])
    #     lines(rng[1:2], rep(0,2), lwd = 3, col = state_colors[txt])
    x = mean(c(s, e))
    text(x = x, y = ylim[1] + yrange / 2 + (.05) * alt, labels = txt, adj = c(.5,ifelse(alt == 1, 0, 1)), cex = .6)
    if(alt == 1){
      alt <<- -1
    }else{
      alt <<- 1
    }
  })
}