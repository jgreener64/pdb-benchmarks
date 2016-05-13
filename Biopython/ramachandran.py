# Benchmark the calculation of Ramachandran phi/psi angles from a PDB file

import time
from Bio.PDB import PDBParser
from Bio.PDB.Vector import calc_dihedral

pdb_filepath = "pdbs/1AKE.pdb"
runs = 100

parser = PDBParser()
struc = parser.get_structure('test', pdb_filepath)

def multiramachandran(n):
    for i in range(n):
        phi_angles = []
        psi_angles = []
        residues = list(struc.get_residues())
        for i in range(1,len(residues)-1):
            res = residues[i]
            res_prev = residues[i-1]
            res_next = residues[i+1]
            # Check residues are not hetero residues and that they have sequential residue numbers
            if res.get_id()[0] == " " and res_prev.get_id()[0] == " " and res_next.get_id()[0] == " " and res.get_id()[1] == res_prev.get_id()[1]+1 and res_next.get_id()[1] == res.get_id()[1]+1:
                phi_angles.append(calc_dihedral(res_prev["C"].get_vector(), res["N"].get_vector(), res["CA"].get_vector(), res["C"].get_vector()))
                psi_angles.append(calc_dihedral(res["N"].get_vector(), res["CA"].get_vector(), res["C"].get_vector(), res_next["N"].get_vector()))

start = time.time()
multiramachandran(runs)
elapsed = time.time() - start

print "Average time per run:", elapsed / runs
