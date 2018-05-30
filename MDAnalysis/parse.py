# Benchmark the parsing of a PDB file given as an argument

import sys
import time
import MDAnalysis as mda

pdb_filepath = sys.argv[1]

start = time.time()
mda.Universe(pdb_filepath)
end = time.time()

print(end - start)
