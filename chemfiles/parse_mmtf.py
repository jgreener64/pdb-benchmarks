# Benchmark the parsing of a MMTF file given as an argument

import sys
import time
from chemfiles import Trajectory

mmtf_filepath = sys.argv[1]

start = time.time()
Trajectory(mmtf_filepath).read()
end = time.time()

print(end - start)
