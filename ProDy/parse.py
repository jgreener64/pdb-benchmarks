# Benchmark the parsing of a PDB file given as an argument

import sys
import time
from prody import *

pdb_filepath = sys.argv[1]
runs = 100

times = []

for i in range(runs):
    start = time.time()
    struc = parsePDB(pdb_filepath)
    elapsed = time.time() - start
    times.append(elapsed)

print "Average time per run:", sum(times) / runs
