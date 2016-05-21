# Benchmark the counting of alanine residues in a PDB file

using Bio.Structure

pdb_filepath = "pdbs/1AKE.pdb"
struc = read(pdb_filepath, PDB)

function count()
    alanineselector(res::AbstractResidue) = resnameselector(res, ["ALA"])
    return countresidues(struc, alanineselector)
end

# Run to JIT compile
count()

elapsed = @elapsed count()

println(elapsed)
