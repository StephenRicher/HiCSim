#!/usr/bin/env Rscript

# Set .libPaths to last path (the Conda environemnt) to avoid conflicts
paths = .libPaths()
lastPath = paths[length(paths)]
.libPaths(lastPath)

library(hicrep)

args = commandArgs(trailingOnly=TRUE)

bin = as.integer(args[1])
start = as.integer(args[2])
end = as.integer(args[3])
matrix1 = args[4]
matrix2 = args[5]

h_hat = NA
scc = NA

m1 = read.table(matrix1)
m2 = read.table(matrix2)

# Get sequencing depth of experimental input matrix
expSum = sum(m1[,-c(1:3)])

# Adjust simulated matrix to match sequencing depth
m2 <- depth.adj(m2, expSum, bin, out = 0)

write(paste("sample1", "sample2", "maxInteractionDistance", "h", "scc", sep="\t"), stdout())

# Set max interaction as half capture region size, or 1,000,000bp
maxInteraction = as.integer((end - start)/2)
for (interactionDistance in seq(bin * 3, maxInteraction, by=(bin * 2))) {
  # Calculate optimal smoothing parameter for a given region and bin
  h_hat = tryCatch({
    htrain(m1, m2, resol = bin, max = interactionDistance, range = 0:20)
  }, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
  
  if (is.null(h_hat)) {
    print("Unable to calculate smoothing parameter - SCC not computed")
    h_hat = NA
    scc = NA
  } else {
    
    # Format HiC matrix pairs and smooth
    pre_hic <- prep(m1, m2, resol = bin, h = h_hat, max = interactionDistance)
    
    # Calculate SCC
    scc = get.scc(pre_hic, resol = bin, max = interactionDistance)$scc
  }
  
  write(paste(matrix1, matrix2, interactionDistance, h_hat, scc, sep="\t"), stdout())
}
