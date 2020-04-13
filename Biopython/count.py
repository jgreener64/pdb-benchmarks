# Benchmark the counting of alanine residues in a PDB file

import time
from Bio.PDB import PDBParser

pdb_filepath = "data/1AKE.pdb"
parser = PDBParser()
struc = parser.get_structure("", pdb_filepath)

def count():
    count = 0
    for res in struc.get_residues():
        if res.get_resname() == "ALA":
            count += 1
    return count

start = time.time()
count()
end = time.time()

print(end - start)
