# Test which PDB entries error on PDB/mmCIF parsers
# Writes output to a file labelled with the week

using BioStructures

start = now()
basedir = "."
pdblist = pdbentrylist()

outstrs = ["Checking all PDB entries at $(now())",
        "Checking $(length(pdblist)) entries"]

for p in sort(pdblist)
    try
        downloadpdb(p, pdb_dir=basedir, file_format=PDB)
    catch
        # Not having a PDB file is acceptable, though a failure to download an
        #   available file may hide an error in parsing
        rm("$basedir/$p.pdb", force=true)
    end
    if isfile("$basedir/$p.pdb")
        try
            s = read("$basedir/$p.pdb", PDB)
        catch
            push!(outstrs, "$p - PDB parsing error")
        end
        rm("$basedir/$p.pdb")
    end
    try
        downloadpdb(p, pdb_dir=basedir, file_format=MMCIF)
    catch
        rm("$basedir/$p.cif", force=true)
        push!(outstrs, "$p - no mmCIF download")
    end
    if isfile("$basedir/$p.cif")
        try
            s = read("$basedir/$p.cif", MMCIF)
        catch
            push!(outstrs, "$p - mmCIF parsing error")
        end
        rm("$basedir/$p.cif")
    end
end

if length(outstrs) == 2
    push!(outstrs, "All entries read fine")
end

push!(outstrs, "Time taken - $(ceil(now() - start, Dates.Minute))")

datestr = replace(string(Date(now())), "-", "")
# This overwrites any existing file
open("$basedir/wholepdb_jl_$datestr.txt", "w") do f
    for l in outstrs
        println(f, l)
    end
end
