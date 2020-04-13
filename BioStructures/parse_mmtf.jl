# Benchmark the parsing of a MMTF file given as an argument

using BioStructures

mmtf_filepath = ARGS[1]

# Run to JIT compile
read(mmtf_filepath, MMTF)

elapsed = @elapsed struc = read(mmtf_filepath, MMTF)

println(elapsed)
