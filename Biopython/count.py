# Benchmark the counting of alanine residues in a PDB file

import time
from Bio.PDB import PDBParser

pdb_filepath = "pdbs/1AKE.pdb"
runs = 1000

parser = PDBParser()
struc = parser.get_structure('test', pdb_filepath)

def multicount(n):
    for i in range(n):
        count = 0
        for res in struc.get_residues():
            if res.get_resname() == "ALA":
                count += 1

start = time.time()
multicount(runs)
elapsed = time.time() - start

print "Average time per run:", elapsed / runs
