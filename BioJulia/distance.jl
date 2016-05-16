# Benchmark the calculation of a distance in a PDB file
# The distance is the closest distance between any atoms of residues 50 and 60
#   in 1AKE

using Bio.Structure

pdb_filepath = "pdbs/1AKE.pdb"
runs = 1

struc = read(pdb_filepath, PDB)
times = Float64[]

# Run to JIT compile
d = distance(struc['A'][50], struc['A'][60])

for i in 1:runs
    tic()
    d = distance(struc['A'][50], struc['A'][60])
    elapsed = toq()
    push!(times, elapsed)
end

println("Average time per run: ", mean(times))
