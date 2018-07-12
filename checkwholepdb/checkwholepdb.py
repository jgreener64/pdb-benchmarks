# Test which PDB entries error on PDB/mmCIF parsers
# Writes output to a file labelled with the week

import os
from datetime import datetime
from math import ceil
from Bio.PDB import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmtf import MMTFParser

start = datetime.now()
basedir = "."
pdbl = PDBList()
pdblist = pdbl.get_all_entries()

outstrs = ["Checking all PDB entries at {}".format(start.isoformat()),
        "Checking {} entries".format(len(pdblist))]

pdb_parser = PDBParser()
mmcif_parser = MMCIFParser()

for pu in sorted(pdblist):
    p = pu.lower()
    try:
        pdbl.retrieve_pdb_file(p, pdir=basedir, file_format="pdb")
    except:
        # Not having a PDB file is acceptable, though a failure to download an
        #   available file may hide an error in parsing
        try:
            os.remove("{}/pdb{}.ent".format(basedir, p))
        except:
            pass
    if os.path.isfile("{}/pdb{}.ent".format(basedir, p)):
        try:
            s = pdb_parser.get_structure("", "{}/pdb{}.ent".format(basedir, p))
        except:
            outstrs.append("{} - PDB parsing error".format(pu))
        os.remove("{}/pdb{}.ent".format(basedir, p))
    try:
        pdbl.retrieve_pdb_file(p, pdir=basedir, file_format="mmCif")
    except:
        try:
            os.remove("{}/{}.cif".format(basedir, p))
        except:
            pass
        outstrs.append("{} - no mmCIF download".format(pu))
    if os.path.isfile("{}/{}.cif".format(basedir, p)):
        try:
            s = mmcif_parser.get_structure("", "{}/{}.cif".format(basedir, p))
        except:
            outstrs.append("{} - mmCIF parsing error".format(pu))
        os.remove("{}/{}.cif".format(basedir, p))

if len(outstrs) == 2:
    outstrs.append("All entries read fine")

end = datetime.now()
outstrs.append("Time taken - {} minute(s)".format(int(ceil((end - start).seconds / 60))))

datestr = str(end.date()).replace("-", "")
# This overwrites any existing file
with open("{}/wholepdb_py_{}.txt".format(basedir, datestr), "w") as f:
    for l in outstrs:
        f.write(l + "\n")
