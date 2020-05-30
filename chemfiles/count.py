# Benchmark the counting of alanine residues in a PDB file

import time
from chemfiles import Trajectory, Selection


def count(frame):
    selection = Selection("resname ALA")
    return len(selection.evaluate(frame))


pdb_filepath = "data/1AKE.pdb"
frame = Trajectory(pdb_filepath).read()

start = time.time()
count(frame)
end = time.time()

print(end - start)
