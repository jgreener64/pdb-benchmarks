# Benchmark the parsing of a mmCIF file given as an argument

import sys
import time
from chemfiles import Trajectory

mmcif_filepath = sys.argv[1]

start = time.time()
trajectory = Trajectory(mmcif_filepath, 'r', 'mmCIF')
for frame in trajectory:
    pass
end = time.time()

print(end - start)
