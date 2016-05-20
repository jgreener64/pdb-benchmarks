# Benchmark the calculation of Ramachandran phi/psi angles from a PDB file

import time
from prody import *

pdb_filepath = "pdbs/1AKE.pdb"
struc = parsePDB(pdb_filepath)

def ramachandran():
    phi_angles = []
    psi_angles = []
    for res in struc.getHierView().iterResidues():
        try:
            phi_angle = calcPhi(res)
            psi_angle = calcPsi(res)
            phi_angles.append(phi_angle)
            psi_angles.append(psi_angle)
        except:
            pass
    return phi_angles, psi_angles

start = time.time()
ramachandran()
elapsed = time.time() - start

print elapsed
