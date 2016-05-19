# Benchmark the counting of alanine residues in a PDB file

import time
from prody import *

pdb_filepath = "pdbs/1AKE.pdb"
runs = 1000

struc = parsePDB(pdb_filepath)
times = []

def count():
    hv = struc.getHierView()
    count = 0
    for res in hv.iterResidues():
        if res.getResname() == "ALA":
            count += 1
    return count

for i in range(runs):
    start = time.time()
    c = count()
    elapsed = time.time() - start
    times.append(elapsed)

print "Average time per run:", sum(times) / runs
