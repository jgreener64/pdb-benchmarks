# Benchmark the parsing of a PDB file given as an argument

import sys
import time
from Bio.PDB import PDBParser

pdb_filepath = sys.argv[1]
runs = 100

parser = PDBParser()

def multiparse(n):
    for i in range(n):
        struc = parser.get_structure('test', pdb_filepath)

start = time.time()
multiparse(runs)
elapsed = time.time() - start

print "Average time per run:", elapsed / runs
