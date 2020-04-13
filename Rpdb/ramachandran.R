# Benchmark the calculation of Ramachandran phi/psi angles from a PDB file

library(Rpdb)
library(microbenchmark)

pdb_filepath <- "data/1AKE.pdb"
struc <- read.pdb(pdb_filepath)

ramachandran <- function() {
    is_n <- which(struc$atoms$elename=="N")
    is_ca <- which(struc$atoms$elename=="CA")
    is_c <- which(struc$atoms$elename=="C")
    res_count <- length(is_ca)
    phi_angles <- dihedral(struc, is_c[1:res_count-2], is_n[2:res_count-1], is_ca[2:res_count-1], is_c[2:res_count-1])
    psi_angles <- dihedral(struc, is_n[2:res_count-1], is_ca[2:res_count-1], is_c[2:res_count-1], is_n[3:res_count])
    return(phi_angles, psi_angles)
}

bench <- microbenchmark(ramachandran(), times=1)

cat(bench$time / 10^9, "\n", sep="")
