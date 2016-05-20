# Benchmark the calculation of a distance in a PDB file
# The distance is the closest distance between any atoms of residues 50 and 60
#   of chain A in 1AKE

import time
from Bio.PDB import PDBParser

pdb_filepath = "pdbs/1AKE.pdb"
parser = PDBParser()
struc = parser.get_structure("", pdb_filepath)

def distance():
    min_dist = float("inf")
    for atom_a in struc[0]['A'][50]:
        for atom_b in struc[0]['A'][60]:
            if atom_a - atom_b < min_dist:
                min_dist = atom_a - atom_b
    return min_dist

start = time.time()
distance()
elapsed = time.time() - start

print elapsed
