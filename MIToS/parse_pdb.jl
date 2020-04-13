# Benchmark the parsing of a PDB file given as an argument

using MIToS.PDB

pdb_filepath = ARGS[1]

# Run to JIT compile
read(pdb_filepath, PDBFile)

elapsed = @elapsed struc = read(pdb_filepath, PDBFile)

println(elapsed)
