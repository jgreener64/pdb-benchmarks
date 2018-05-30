# Benchmark the calculation of a distance in a PDB file
# The distance is the closest distance between any atoms of residues 50 and 60
#   of chain A in 1AKE

using BioStructures

pdb_filepath = "pdbs/1AKE.pdb"
struc = read(pdb_filepath, PDB)

# Run to JIT compile
distance(struc['A'][50], struc['A'][60])

elapsed = @elapsed distance(struc['A'][50], struc['A'][60])

println(elapsed)
