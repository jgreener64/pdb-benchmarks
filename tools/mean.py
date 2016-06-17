# Calculate the mean value in a file of numbers

import sys
import numpy as np

in_file = open(sys.argv[1])
vals = []
for line in in_file:
    vals.append(float(line))
in_file.close()

print(np.mean(vals))

