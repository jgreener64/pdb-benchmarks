# Benchmark the calculation of Ramachandran phi/psi angles from a PDB file

import time
import MDAnalysis as mda

pdb_filepath = "data/1AKE.pdb"
u = mda.Universe(pdb_filepath)

def ramachandran():
    phi_angles = []
    psi_angles = []
    for res in u.residues:
        try:
            phi = res.phi_selection()
        except:
            pass
        else:
            if not phi is None:
                phi_angles.append(phi.dihedral.value())
        try:
            psi = res.psi_selection()
        except:
            pass
        else:
            if not psi is None:
                psi_angles.append(psi.dihedral.value())
    return phi_angles, psi_angles

start = time.time()
ramachandran()
end = time.time()

print(end - start)
