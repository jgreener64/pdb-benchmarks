# Benchmark the parsing of a PDB file given as an argument

using Bio.Structure

const pdb_filepath = ARGS[1]
runs = 100

function multiparse(n::Integer)
    for i in 1:n
        struc = read(pdb_filepath, PDB)
    end
end

# Run to JIT compile
multiparse(1)

tic()
multiparse(runs)
elapsed = toq()

println("Average time per run: ", elapsed / runs)
