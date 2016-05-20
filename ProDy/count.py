# Benchmark the counting of alanine residues in a PDB file

import time
from prody import *

pdb_filepath = "pdbs/1AKE.pdb"
struc = parsePDB(pdb_filepath)

def count():
    count = 0
    for res in struc.getHierView().iterResidues():
        if res.getResname() == "ALA":
            count += 1
    return count

start = time.time()
count()
elapsed = time.time() - start

print elapsed
