# Benchmark the calculation of Ramachandran phi/psi angles from a PDB file

using Bio.Structure

pdb_filepath = "pdbs/1AKE.pdb"
runs = 100

struc = read(pdb_filepath, PDB)
times = Float64[]

# Run to JIT compile
phi, psi = ramachandranangles(struc, stdresselector)

for i in 1:runs
    tic()
    phi, psi = ramachandranangles(struc, stdresselector)
    elapsed = toq()
    push!(times, elapsed)
end

println("Average time per run: ", mean(times))
