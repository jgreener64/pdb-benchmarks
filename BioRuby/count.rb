# Benchmark the counting of alanine residues in a PDB file

require "bio"
require "benchmark"

pdb_filepath = "pdbs/1AKE.pdb"
runs = 1

times = 0.0
pdb = Bio::PDB.new(File.read(pdb_filepath))

for i in 1..runs
    elapsed = Benchmark.realtime {
        pdb.find_residue { |res| res.resName == "ALA" }.length
    }
    times = times + elapsed
end

print "Average time per run: ", times / runs, "\n"
