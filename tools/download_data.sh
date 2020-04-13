# Download the PDB files required for benchmarking into data directory

if [ ! -d data ]; then
    mkdir data
fi

for pdbid in 1CRN 3JYV 1HTQ 1AKE
do
    echo "Downloading $pdbid"
    wget "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=$pdbid" -O data/$pdbid.pdb
done
