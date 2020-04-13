# Benchmark the parsing of a PDB file given as an argument

library(bio3d)
library(microbenchmark)

pdb_filepath <- commandArgs(trailingOnly=TRUE)[1]

bench <- microbenchmark(read.pdb(pdb_filepath, multi=TRUE), times=1)

cat(bench$time / 10^9, "\n", sep="")
