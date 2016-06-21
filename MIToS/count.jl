# Benchmark the counting of alanine residues in a PDB file

using MIToS.PDB

pdb_filepath = "pdbs/1AKE.pdb"
struc = read(pdb_filepath, PDBFile)

count() = length(findobjects(struc, Is(:name, "ALA")))

# Run to JIT compile
count()

elapsed = @elapsed count()

println(elapsed)
