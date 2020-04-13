# Benchmark the parsing of a PDB file given as an argument

import sys
import time
from Bio.PDB import PDBParser

pdb_filepath = sys.argv[1]
parser = PDBParser()

start = time.time()
parser.get_structure("", pdb_filepath)
end = time.time()

print(end - start)
