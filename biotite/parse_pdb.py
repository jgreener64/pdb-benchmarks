# Benchmark the parsing of a PDB file given as an argument

import sys
import time
import biotite.structure.io.pdb as pdb

pdb_filepath = sys.argv[1]

start = time.time()
file = pdb.PDBFile()
file.read(pdb_filepath)
file.get_structure()
end = time.time()

print(end - start)
