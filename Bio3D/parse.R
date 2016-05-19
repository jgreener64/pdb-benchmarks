# Benchmark the parsing of a PDB file given as an argument

library(bio3d)
library(microbenchmark)

pdb_filepath <- commandArgs(trailingOnly=TRUE)[1]
runs <- 3

parsepdb <- function() {
    struc <- read.pdb(pdb_filepath, multi=TRUE)
}

bench <- microbenchmark(parsepdb(), times=runs)

cat("Average time per run: ", mean(bench$time) / 10^9, "\n", sep="")
