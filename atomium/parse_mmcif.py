# Benchmark the parsing of a mmCIF file given as an argument

import sys
import time
import atomium

mmcif_filepath = sys.argv[1]

start = time.time()
atomium.open(mmcif_filepath)
end = time.time()

print(end - start)
