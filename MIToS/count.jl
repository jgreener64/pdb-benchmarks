# Benchmark the counting of alanine residues in a PDB file

using MIToS.PDB

pdb_filepath = "pdbs/1AKE.pdb"
struc = read(pdb_filepath, PDBFile)

counter() = count(res -> res.id.name == "ALA", struc)

# Run to JIT compile
counter()

elapsed = @elapsed counter()

println(elapsed)
