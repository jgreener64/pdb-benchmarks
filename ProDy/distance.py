# Benchmark the calculation of a distance in a PDB file
# The distance is the closest distance between any atoms of residues 50 and 60
#   of chain A in 1AKE

import time
from prody import *

pdb_filepath = "data/1AKE.pdb"
struc = parsePDB(pdb_filepath)

def distance():
    min_dist = float("inf")
    for atom_a in struc['A', 50]:
        for atom_b in struc['A', 60]:
            if calcDistance(atom_a, atom_b) < min_dist:
                min_dist = calcDistance(atom_a, atom_b)
    return min_dist

start = time.time()
distance()
end = time.time()

print(end - start)
