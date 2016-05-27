source("plot.R")
source("fetch_states.R")
readPara("mcf10a_and_mcf7_ideas_full", 1, 4)
rP$mat = matrix(rP$m, ncol = 4, byrow = T)
colnames(rP$mat) = c("H3K27ac", "H3K27me3", "H3K4ac", "H3K4me3")
rownames(rP$mat) = 1:nrow(rP$mat)

source("H:/R_workspace/jrb_R_scripts/heatmap.3-split.R")
source("H:/R_workspace/jrb_R_scripts/heatmap.3-kmeans_wrapper.R")

plot_mat = rP$mat
MAX = quantile(plot_mat, .95)
plot_mat = ifelse(plot_mat > MAX, MAX, plot_mat)

cr = colorRamp(c("gray", "blue", "red"))
vals = 0:200 /100

cr = colorRamp(c("white", "blue"))
vals = 0:100 /100


plot_colors = rgb(cr(ifelse(vals >1 , 1, vals)) / 255)

heatmap.2(plot_mat, trace = "n", col = plot_colors, scale = "r")

read.table(file = "mcf10a_and_mcf7_ideas_full.state", nrows = 10)

states = read.table(file = "mcf10a_and_mcf7_ideas_full.state", stringsAsFactors = F)
states[,5:ncol(states)] = states[,5:ncol(states)] + 1
states[,4] = states[,4] + 200
drugs = c("bza", "ctrl", "e2bza", "e2", "gc10bza", "gc10")
cnames = sapply(c("MCF10A", "MCF7"), function(x)paste(x, drugs, sep = "_"))[1:12]
colnames(states) = c("row", "chrm", "start", "end", cnames)
