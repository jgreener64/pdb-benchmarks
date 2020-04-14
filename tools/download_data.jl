# Download the data files required for benchmarking into data directory

using BioStructures

out_dir = "data"

if !isdir(out_dir)
    mkdir(out_dir)
end

for pdbid in ("1CRN", "1HTQ")
    for format in (PDB, MMCIF, MMTF)
        downloadpdb(pdbid, format=format, dir=out_dir)
    end
end

downloadpdb("1AKE", dir=out_dir)
