# Benchmark the parsing of a mmCIF file given as an argument

import sys
import time
from Bio.PDB import MMCIFParser

mmcif_filepath = sys.argv[1]
parser = MMCIFParser()

start = time.time()
parser.get_structure("", mmcif_filepath)
end = time.time()

print(end - start)
