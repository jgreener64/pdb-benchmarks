# Benchmark the parsing of a MMTF file given as an argument

import sys
import time
from Bio.PDB.mmtf import MMTFParser

mmtf_filepath = sys.argv[1]

start = time.time()
MMTFParser.get_structure(mmtf_filepath)
end = time.time()

print(end - start)
