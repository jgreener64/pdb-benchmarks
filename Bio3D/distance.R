# Benchmark the calculation of a distance in a PDB file
# The distance is the closest distance between any atoms of residues 50 and 60
#   of chain A in 1AKE

library(bio3d)
library(microbenchmark)

pdb_filepath <- "pdbs/1AKE.pdb"
struc <- read.pdb(pdb_filepath)

distance <- function() {
    coords <- matrix(struc$xyz, length(struc$xyz) / 3, 3, byrow=TRUE)
    is_res50 <- which(struc$atom$resno == 50 & struc$atom$chain == "A")
    is_res60 <- which(struc$atom$resno == 60 & struc$atom$chain == "A")
    return(min(dist.xyz(coords[is_res50,], coords[is_res60,])))
}

bench <- microbenchmark(distance(), times=1)

cat(bench$time / 10^9, "\n", sep="")
