# Benchmark the parsing of a PDB file given as an argument

require 'bio'
require 'benchmark'

pdb_filepath = ARGV[0]
runs = 1

times = 0.0

for i in 1..runs
    elapsed = Benchmark.realtime {
        pdb = Bio::PDB.new(File.read(pdb_filepath))
    }
    times = times + elapsed
end

print "Average time per run: ", times / runs, "\n"
