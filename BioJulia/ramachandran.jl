# Benchmark the calculation of Ramachandran phi/psi angles from a PDB file

using Bio.Structure

pdb_filepath = "pdbs/1AKE.pdb"
runs = 100

struc = read(pdb_filepath, PDB)

function multiramachandran(n::Integer)
    for i in 1:n
        phi_angles, psi_angles = ramachandranangles(struc, stdresselector)
    end
end

# Run to JIT compile
multiramachandran(1)

tic()
multiramachandran(runs)
elapsed = toq()

println("Average time per run: ", elapsed / runs)
