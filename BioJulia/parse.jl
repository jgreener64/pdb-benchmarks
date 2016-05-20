# Benchmark the parsing of a PDB file given as an argument

using Bio.Structure

pdb_filepath = ARGS[1]

# Run to JIT compile
read(pdb_filepath, PDB)

elapsed = @elapsed struc = read(pdb_filepath, PDB)

println(elapsed)
