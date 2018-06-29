# Benchmark the counting of alanine residues in a PDB file

using BioStructures

pdb_filepath = "pdbs/1AKE.pdb"
struc = read(pdb_filepath, PDB)

function counter()
    alanineselector(res::AbstractResidue) = resnameselector(res, ["ALA"])
    return countresidues(struc, alanineselector)
end

# Run to JIT compile
counter()

elapsed = @elapsed counter()

println(elapsed)
