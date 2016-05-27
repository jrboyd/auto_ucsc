source("plot.R")
clust = read.table("mcf7_ideas_full.cluster", stringsAsFactors = F)
if(!exists("states")) states = read.table("mcf7_ideas_full.state", stringsAsFactors = F)
rP = readPara("mcf7_ideas_full", 1, 4)
mat = rP$m
mat = matrix(mat, nrow = 4, byrow = F)
mat = t(mat)
marks = c("27ac", "27me3", "4ac", "k4me3")
colnames(mat) = marks

head(states)
o = c(2,4,1,3,6,5)
cell_lines = c("bza", "ctrl", "e2bza", "e2", "gc10bza", "gc10")
cell_lines = cell_lines[o]

states[,5:10] = states[,(5:10)[o]] + 1#states are now equal to indexes in parameters
states[,4] = states[,4] + 199

colnames(states) = c("index", "chrm", "start", "end", cell_lines, "position_state")

if(F){#write bed files for ngsplotting
  for(col_num in 5:11){
    rng = sort(unique(states[,col_num]))
    for(state_num in rng){
      keep = states[,col_num] == state_num
      st = states[keep,]
      MAX = 50000
      if(nrow(st) > MAX){#cap size of output to randomly selected MAX
        r = runif(nrow(st))
        k = order(r)[1:MAX]
        st = st[k,]
      }
      fname = paste("MCF7", colnames(st)[col_num], "state", state_num, "ngsplot.bed", sep = "_")
      print(paste(fname, nrow(st)))
      write.table(st[,2:4], file = fname, sep = "\t", quote = F, col.names = F, row.names = F)
    }
  }
}

for(i in 1:nrow(mat)){
  print(("-------------"))
  print(paste("---STATE is", i))
  keep = apply(states[,5:10], 1, function(x)all(x == i))
  st = states[keep,]
  print("---MEAN")
  print(mat[i,])
  print(paste("---RANDOM 5 of", sum(keep)))
  print(st[order(runif(nrow(st)))[1:5],])
}



state_names = character(nrow(mat))
state_names[1] = "weak K4ac and K27me3"
state_names[2] = "K4me3 and K4ac, moderate K27ac"
state_names[3] = "weak K27ac"
state_names[4] = "NA"
state_names[5] = "K4ac and K27ac"
state_names[6] = "moderate K27ac and weak others"
state_names[7] = "strong K4me3, K4ac, and K27ac"
state_names[8] = "strong K4me3, low K4ac and K27ac"
state_names[9] = "blank, very weak all 4 marks"
state_names[10] = "strong K4me3, low K4ac and K27ac, K27me3 nearby"
state_names[11] = "strong K4me3, K4ac and K27ac almost as strong"
state_names[12] = "blank or K4me3 and K4ac"
state_names[13] = "strong K4me3, moderate K4ac and K27ac"
state_names[14] = "strong K4me3, moderate K4ac and K27ac"
state_names[15] = "moderate K4me3, K4ac, and K27ac"
state_names[16] = "blank"
state_names[34] = "blank"

pos_states[12] = "blank"
pos_states[10] = "blank - intermittent weak"
pos_states[17] = "active with K4me3/ac dominant"
pos_states[3] = "active with K27ac dominant"
pos_states[4] = "active with all 3 marks"
pos_states[15] = "promoter with K4me3 and weak acetyl"