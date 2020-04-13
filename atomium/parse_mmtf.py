# Benchmark the parsing of a MMTF file given as an argument

import sys
import time
import atomium

mmtf_filepath = sys.argv[1]

start = time.time()
atomium.open(mmtf_filepath)
end = time.time()

print(end - start)
