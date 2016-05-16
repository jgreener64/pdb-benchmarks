# Benchmark the counting of alanine residues in a PDB file

import time
from Bio.PDB import PDBParser

pdb_filepath = "pdbs/1AKE.pdb"
runs = 100

parser = PDBParser()
struc = parser.get_structure("", pdb_filepath)
times = []

def count():
    count = 0
    for res in struc.get_residues():
        if res.get_resname() == "ALA":
            count += 1
    return count

for i in range(runs):
    start = time.time()
    c = count()
    elapsed = time.time() - start
    times.append(elapsed)

print "Average time per run:", sum(times) / runs
