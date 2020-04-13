# Calculate the mean value from a file of numbers

import sys
import numpy as np

with open(sys.argv[1]) as in_file:
    vals = [float(line) for line in in_file]

print(np.mean(vals))
