# Benchmark the calculation of a distance in a PDB file
# The distance is the closest distance between any atoms of residues 50 and 60
#   of chain A in 1AKE

library(Rpdb)
library(microbenchmark)

pdb_filepath <- "data/1AKE.pdb"
struc <- read.pdb(pdb_filepath)

distance <- function() {
    is.res50 <- struc$atoms$resid == 50 & struc$atoms$chainid == "A"
    is.res60 <- struc$atoms$resid == 60 & struc$atoms$chainid == "A"
    d <- distances(struc, is.res50, is.res60)
    return(min(norm(d, type="xyz")))
}

bench <- microbenchmark(distance(), times=1)

cat(bench$time / 10^9, "\n", sep="")
