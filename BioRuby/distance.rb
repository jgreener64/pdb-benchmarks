# Benchmark the calculation of a distance in a PDB file
# The distance is the closest distance between any atoms of residues 50 and 60
#   of chain A in 1AKE

require "bio"
require "benchmark"
include Bio::PDB::Utils

pdb_filepath = "data/1AKE.pdb"
pdb = Bio::PDB.new(File.read(pdb_filepath))

elapsed = Benchmark.realtime {
    res_50 = pdb.find_residue { |res| res.resSeq == 50 and res.chain.id == "A" }[0]
    res_60 = pdb.find_residue { |res| res.resSeq == 60 and res.chain.id == "A" }[0]
    min_dist = Float::INFINITY
    res_50.each_atom do |atom_50|
        res_60.each_atom do |atom_60|
            if distance(atom_50, atom_60) < min_dist
                min_dist = distance(atom_50, atom_60)
            end
        end
    end
    min_dist
}

print elapsed, "\n"
