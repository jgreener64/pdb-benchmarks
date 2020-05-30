# Benchmark the calculation of a distance in a PDB file
# The distance is the closest distance between any atoms of residues 50 and 60
#   of chain A in 1AKE

import time
from chemfiles import Trajectory, Selection


def distance(frame):
    # FIXME: this should use Selection("resid 50 and [chainname] A") which will
    # be available in chemfiles 0.10 (the next release)
    r50 = Selection("resid 50 and index < 1000").evaluate(frame)
    r60 = Selection("resid 60 and index < 1000").evaluate(frame)

    min = float('inf')
    for i in r50:
        for j in r60:
            r = frame.distance(i, j)
            if r < min:
                min = r

    return min


pdb_filepath = "data/1AKE.pdb"
frame = Trajectory(pdb_filepath).read()

start = time.time()
distance(frame)
end = time.time()

print(end - start)
