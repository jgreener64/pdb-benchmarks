# Print the benchmark results as a markdown table

times = Dict{String, Dict}()

open("benchmarks.csv") do f
    for line in eachline(f)
        if !startswith(line, "Package")
            software, benchmark, runtime = split(line, ",")
            if haskey(times, benchmark)
                times[benchmark][software] = parse(Float64, runtime)
            else
                times[benchmark] = Dict(software=> parse(Float64, runtime))
            end
        end
    end
end

for (benchmark, label, millisecond) in (
                            ("Parse PDB 1CRN"  , "Parse PDB 1CRN / ms"  , true ),
                            ("Parse PDB 1HTQ"  , "Parse PDB 1HTQ / s"   , false),
                            ("Parse mmCIF 1CRN", "Parse mmCIF 1CRN / ms", true ),
                            ("Parse mmCIF 1HTQ", "Parse mmCIF 1HTQ / s" , false),
                            ("Parse MMTF 1CRN" , "Parse MMTF 1CRN / ms" , true ),
                            ("Parse MMTF 1HTQ" , "Parse MMTF 1HTQ / s"  , false),
                            ("Count"           , "Count / ms"           , true ),
                            ("Distance"        , "Distance / ms"        , true ),
                            ("Ramachandran"    , "Ramachandran / ms"    , true ))
    print("| $(rpad(label, 21)) |")
    for software in ("BioStructures", "MIToS", "Biopython", "ProDy", "MDAnalysis", "biotite",
                        "atomium", "Bio3D", "Rpdb", "BioJava", "BioPerl", "BioRuby", "GEMMI",
                        "Victor", "ESBTL")
        if haskey(times[benchmark], software)
            if millisecond
                val = string(round(1000 * times[benchmark][software], sigdigits=2))
            else
                val = string(round(times[benchmark][software], sigdigits=2))
            end
        else
            val = "-"
        end
        print(" $(rpad(val, 13)) |")
    end
    println()
end
