load("ensg_ref.save")
dupe = (duplicated(ensg_ref$gene_name))
sym_ref = ensg_ref[!dupe,]
rownames(sym_ref) = sym_ref$gene_name
all_symbols = sort(unique(sym_ref$gene_name))
#states.IDEAS = read.table("mcf7_ideas_full.state", stringsAsFactors = F)

symbol2range = function(sym){#convert symbol to range (chrm, start, end, strand)
  sym = toupper(sym)
  range = sym_ref[sym,3:6]
  return(range)
}

symbol2ucsc = function(sym){#convert symbol chrm:start-end
  range = symbol2range(sym)
  range = paste0(range$chrm, ":", range$start, "-", range$end)
  return(range)
}

symbol2promoter_ucsc = function(sym, ext = 2000){
  range = symbol2range(sym)
  tss = range$start
  if(range$strand == "-") tss = range$end
  s = tss - ext
  e = tss + ext
  
  return(paste0(range$chrm, ":", s, "-", e))
}

fetch_IDEAS_by_symbol = function(sym, flank = 0){#assume states.IDEAS has been loaded
  range = symbol2range(sym)
  st = fetch_IDEAS_by_range(chrm = range$chrm, start = range$start - flank, end = range$end + flank, strand = range$strand)
  return(st)
}


fetch_IDEAS_by_range = function(chrm, start, end, strand = "+"){
  range = data.frame(chrm = chrm, start = start, end = end, strand = strand)
  st = states.IDEAS[states.IDEAS[,2] == range$chrm, ]
  keep = st[,3] > range$start & st[,4] < range$end
  st = st[keep,]
  if(strand == "-"){
    st = st[nrow(st):1,]
  }
  return(st)
}