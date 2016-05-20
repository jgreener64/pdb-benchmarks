# Benchmark the counting of alanine residues in a PDB file

library(bio3d)
library(microbenchmark)

pdb_filepath <- "pdbs/1AKE.pdb"
runs <- 1000

struc <- read.pdb(pdb_filepath, multi=TRUE)

count <- function() {
    resnums <- struc$atom$resno[struc$atom$resid=="ALA"]
    chains <- struc$atom$chain[struc$atom$resid=="ALA"]
    resids <- paste(resnums, chains, sep="")
    return(length(unique(resids)))
}

bench <- microbenchmark(count(), times=runs)

cat("Average time per run: ", mean(bench$time) / 10^9, "\n", sep="")
