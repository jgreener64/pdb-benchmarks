# Benchmark the counting of alanine residues in a PDB file

library(Rpdb)
library(microbenchmark)

pdb_filepath <- "data/1AKE.pdb"
struc <- read.pdb(pdb_filepath)

count <- function() {
    resnums <- struc$atoms$resid[struc$atoms$resname=="ALA"]
    chains <- struc$atoms$chainid[struc$atoms$resname=="ALA"]
    resids <- paste(resnums, chains, sep="")
    return(length(unique(resids)))
}

bench <- microbenchmark(count(), times=1)

cat(bench$time / 10^9, "\n", sep="")
