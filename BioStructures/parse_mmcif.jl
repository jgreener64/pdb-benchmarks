# Benchmark the parsing of a mmCIF file given as an argument

using BioStructures

mmcif_filepath = ARGS[1]

# Run to JIT compile
read(mmcif_filepath, MMCIF)

elapsed = @elapsed struc = read(mmcif_filepath, MMCIF)

println(elapsed)
