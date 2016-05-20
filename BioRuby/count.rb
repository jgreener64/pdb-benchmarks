# Benchmark the counting of alanine residues in a PDB file

require "bio"
require "benchmark"

pdb_filepath = "pdbs/1AKE.pdb"
pdb = Bio::PDB.new(File.read(pdb_filepath))

elapsed = Benchmark.realtime {
    pdb.find_residue { |res| res.resName == "ALA" }.length
}

print elapsed, "\n"
