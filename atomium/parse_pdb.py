# Benchmark the parsing of a PDB file given as an argument

import sys
import time
import atomium

pdb_filepath = sys.argv[1]

start = time.time()
atomium.open(pdb_filepath)
end = time.time()

print(end - start)
