# Benchmark the counting of alanine residues in a PDB file

library(bio3d)
library(microbenchmark)

pdb_filepath <- "data/1AKE.pdb"
struc <- read.pdb(pdb_filepath, multi=TRUE)

count <- function() {
    resnums <- struc$atom$resno[struc$atom$resid=="ALA"]
    chains <- struc$atom$chain[struc$atom$resid=="ALA"]
    resids <- paste(resnums, chains, sep="")
    return(length(unique(resids)))
}

bench <- microbenchmark(count(), times=1)

cat(bench$time / 10^9, "\n", sep="")
