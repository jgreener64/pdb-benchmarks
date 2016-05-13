# Benchmark the calculation of a distance in a PDB file
# The distance is the closest distance between any atoms of residues 50 and 60
#   in 1AKE

using Bio.Structure

pdb_filepath = "pdbs/1AKE.pdb"
runs = 10000

struc = read(pdb_filepath, PDB)

function multidistance(n::Integer)
    for i in 1:n
        dist = distance(struc['A'][50], struc['A'][60])
    end
end

# Run to JIT compile
multidistance(1)

tic()
multidistance(runs)
elapsed = toq()

println("Average time per run: ", elapsed / runs)
