# Benchmark the parsing of a PDB file given as an argument

import sys
import time
from prody import *

pdb_filepath = sys.argv[1]

start = time.time()
parsePDB(pdb_filepath)
end = time.time()

print(end - start)
