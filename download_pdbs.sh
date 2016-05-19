# Download the PDB files required for benchmarking into directory pdbs

[ ! -d pdbs ] && mkdir pdbs
for pdbid in 1CRN 3JYV 1HTQ 1AKE
do
    echo "Downloading $pdbid"
    wget "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=$pdbid" -O pdbs/$pdbid.pdb
done
