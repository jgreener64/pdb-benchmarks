import sys
import time
import MDAnalysis as mda

pdb_filepath = sys.argv[1]

start = time.time()
mda.Universe(pdb_filepath)
elapsed = time.time() - start

print elapsed
