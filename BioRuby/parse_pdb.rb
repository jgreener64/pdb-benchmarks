# Benchmark the parsing of a PDB file given as an argument

require "bio"
require "benchmark"

pdb_filepath = ARGV[0]

elapsed = Benchmark.realtime {
    Bio::PDB.new(File.read(pdb_filepath))
}

print elapsed, "\n"
