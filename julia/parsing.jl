# Parse a PDB file given as an argument

using Bio.Structure

pdb_filepath = ARGS[1]

read(pdb_filepath, PDB)
@time struc = read(pdb_filepath, PDB)
