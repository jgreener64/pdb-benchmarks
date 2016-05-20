# Benchmark the parsing of a PDB file given as an argument

library(Rpdb)
library(microbenchmark)

pdb_filepath <- commandArgs(trailingOnly=TRUE)[1]

bench <- microbenchmark(read.pdb(pdb_filepath), times=1)

cat(bench$time / 10^9, "\n", sep="")
