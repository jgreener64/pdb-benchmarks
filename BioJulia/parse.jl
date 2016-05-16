# Benchmark the parsing of a PDB file given as an argument

using Bio.Structure

const pdb_filepath = ARGS[1]
runs = 100

times = Float64[]

# Run to JIT compile
struc = read(pdb_filepath, PDB)

for i in 1:runs
    tic()
    struc = read(pdb_filepath, PDB)
    elapsed = toq()
    push!(times, elapsed)
end

println("Average time per run: ", mean(times))
