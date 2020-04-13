# Benchmark the parsing of a MMTF file given as an argument

import sys
import time
import biotite.structure.io.mmtf as mmtf

mmtf_filepath = sys.argv[1]

start = time.time()
file = mmtf.MMTFFile()
file.read(mmtf_filepath)
mmtf.get_structure(file)
end = time.time()

print(end - start)
