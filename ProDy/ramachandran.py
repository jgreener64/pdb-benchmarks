# Benchmark the calculation of Ramachandran phi/psi angles from a PDB file

import time
from prody import *

pdb_filepath = "pdbs/1AKE.pdb"
runs = 10

struc = parsePDB(pdb_filepath)
times = []

def ramachandran():
    phi_angles = []
    psi_angles = []
    hv = struc.getHierView()
    for res in hv.iterResidues():
        try:
            phi_angle = calcPhi(res)
            psi_angle = calcPsi(res)
            phi_angles.append(phi_angle)
            psi_angles.append(psi_angle)
        except:
            pass
    return phi_angles, psi_angles

for i in range(runs):
    start = time.time()
    phi, psi = ramachandran()
    elapsed = time.time() - start
    times.append(elapsed)

print "Average time per run:", sum(times) / runs
