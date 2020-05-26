# Benchmark the parsing of a PDB file given as an argument

import sys
import time
from chemfiles import Trajectory

pdb_filepath = sys.argv[1]

start = time.time()
Trajectory(pdb_filepath).read()
end = time.time()

print(end - start)
