# Benchmark the counting of alanine residues in a PDB file

import time
import MDAnalysis as mda

pdb_filepath = "pdbs/1AKE.pdb"
u = mda.Universe(pdb_filepath)

def count():
    return (u.residues.resnames == "ALA").sum()

start = time.time()
count()
elapsed = time.time() - start

print elapsed
