# Benchmark the counting of alanine residues in a PDB file

import time
import MDAnalysis as mda

pdb_filepath = "data/1AKE.pdb"
u = mda.Universe(pdb_filepath)

def count():
    return (u.residues.resnames == "ALA").sum()

start = time.time()
count()
end = time.time()

print(end - start)
