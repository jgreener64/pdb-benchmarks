import time
import MDAnalysis as mda

u = mda.Universe("pdbs/4AKE.pdb")

def count():
    """ah ah ah ah ah"""
    return (u.residues.resnames == "ALA").sum()

start = time.time()
count()
elapsed = time.time() - start

print elapsed
