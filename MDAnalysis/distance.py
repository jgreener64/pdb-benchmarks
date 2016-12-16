# Benchmark the calculation of a distance in a PDB file
# The distance is the closest distance between any atoms of residues 50 and 60
#   of chain A in 1AKE

import time
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

pdb_filepath = "pdbs/1AKE.pdb"
u = mda.Universe(pdb_filepath)

def distance():
    segA = u.segments[0]
    r50 = segA.atoms.select_atoms("resid 50")
    r60 = segA.atoms.select_atoms("resid 60")
    da = distance_array(r50.positions, r60.positions)
    return da.min()

start = time.time()
distance()
elapsed = time.time() - start

print elapsed
