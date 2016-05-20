# Benchmark the calculation of Ramachandran phi/psi angles from a PDB file

using Bio.Structure

pdb_filepath = "pdbs/1AKE.pdb"
struc = read(pdb_filepath, PDB)

# Run to JIT compile
ramachandranangles(struc, stdresselector)

elapsed = @elapsed ramachandranangles(struc, stdresselector)

println(elapsed)
