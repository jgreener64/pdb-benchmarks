# Benchmark the counting of alanine residues in a PDB file

using Bio.Structure

pdb_filepath = "pdbs/1AKE.pdb"
runs = 1000

struc = read(pdb_filepath, PDB)

function multicount(n::Integer)
    for i in 1:n
        count = countresidues(struc, res -> resname(res) == "ALA")
    end
end

# Run to JIT compile
multicount(1)

tic()
multicount(runs)
elapsed = toq()

println("Average time per run: ", elapsed / runs)
