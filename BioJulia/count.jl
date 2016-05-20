# Benchmark the counting of alanine residues in a PDB file

using Bio.Structure

pdb_filepath = "pdbs/1AKE.pdb"
struc = read(pdb_filepath, PDB)

# Run to JIT compile
countresidues(struc, res -> resname(res) == "ALA")

elapsed = @elapsed countresidues(struc, res -> resname(res) == "ALA")

println(elapsed)
