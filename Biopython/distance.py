# Benchmark the calculation of a distance in a PDB file
# The distance is the closest distance between any atoms of residues 50 and 60
#   in 1AKE

import time
from Bio.PDB import PDBParser

pdb_filepath = "pdbs/1AKE.pdb"
runs = 10000

parser = PDBParser()
struc = parser.get_structure('test', pdb_filepath)

def multidistance(n):
    for i in range(n):
        min_dist = float("inf")
        for atom_a in struc[0]['A'][50]:
            for atom_b in struc[0]['A'][60]:
                if atom_a - atom_b < min_dist:
                    min_dist = atom_a - atom_b

start = time.time()
multidistance(runs)
elapsed = time.time() - start

print "Average time per run:", elapsed / runs
