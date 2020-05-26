# Benchmark the parsing of a mmCIF file given as an argument

import sys
import time
from chemfiles import Trajectory

mmcif_filepath = sys.argv[1]

start = time.time()
Trajectory(mmcif_filepath, 'r', 'mmCIF').read()
end = time.time()

print(end - start)
