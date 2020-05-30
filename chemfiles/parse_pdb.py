# Benchmark the parsing of a PDB file given as an argument

import sys
import time
from chemfiles import Trajectory

pdb_filepath = sys.argv[1]

start = time.time()
trajectory = Trajectory(pdb_filepath)
for frame in trajectory:
    pass
end = time.time()

print(end - start)
