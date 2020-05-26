# Benchmark the calculation of Ramachandran phi/psi angles from a PDB file
import time
from chemfiles import Trajectory, Selection


def ramachandran(frame):
    phi_selection = Selection("dihedrals: name(#1) C and name(#2) N and name(#3) CA and name(#4) C")
    phi_angles = []
    for (i, j, k, m) in phi_selection.evaluate(frame):
        # 57.29578 to convert from radians to degrees
        phi_angles.append(frame.dihedral(i, j, k, m) * 57.29578)

    psi_selection = Selection("dihedrals: name(#1) N and name(#2) CA and name(#3) C and name(#4) N")
    psi_angles = []
    for (i, j, k, m) in psi_selection.evaluate(frame):
        psi_angles.append(frame.dihedral(i, j, k, m) * 57.29578)

    # FIXME: the sign of the angles is inverted w.r.t. the MDAnalysis results
    return phi_angles, psi_angles


pdb_filepath = "data/1AKE.pdb"
frame = Trajectory(pdb_filepath).read()

start = time.time()
ramachandran(frame)
end = time.time()

print(end - start)
