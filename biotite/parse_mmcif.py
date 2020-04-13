# Benchmark the parsing of a mmCIF file given as an argument

import sys
import time
import biotite.structure.io.pdbx as pdbx

mmcif_filepath = sys.argv[1]

start = time.time()
file = pdbx.PDBxFile()
file.read(mmcif_filepath)
pdbx.get_structure(file)
end = time.time()

print(end - start)
