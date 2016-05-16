# Benchmark the counting of alanine residues in a PDB file

using Bio.Structure

pdb_filepath = "pdbs/1AKE.pdb"
runs = 100

struc = read(pdb_filepath, PDB)
times = Float64[]

# Run to JIT compile
c = countresidues(struc, res -> resname(res) == "ALA")

for i in 1:runs
    tic()
    c = countresidues(struc, res -> resname(res) == "ALA")
    elapsed = toq()
    push!(times, elapsed)
end

println("Average time per run: ", mean(times))
