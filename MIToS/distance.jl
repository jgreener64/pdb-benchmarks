# Benchmark the calculation of a distance in a PDB file
# The distance is the closest distance between any atoms of residues 50 and 60
#   of chain A in 1AKE

using MIToS.PDB

pdb_filepath = "pdbs/1AKE.pdb"
struc = read(pdb_filepath, PDBFile, model="1", chain="A", group="ATOM")

# Run to JIT compile
distance(struc[50], struc[60])

elapsed = @elapsed distance(struc[50], struc[60])

println(elapsed)
