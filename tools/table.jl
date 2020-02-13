# Print the benchmark results as a markdown table

times = Dict{String, Dict}()

open("benchmarks.csv") do f
    for line in eachline(f)
        if !startswith(line, "package")
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
                            ("parse 1CRN"  , "Parse 1CRN / ms"  , true ),
                            ("parse 3JYV"  , "Parse 3JYV / s"   , false),
                            ("parse 1HTQ"  , "Parse 1HTQ / s"   , false),
                            ("count"       , "Count / ms"       , true ),
                            ("distance"    , "Distance / ms"    , true ),
                            ("ramachandran", "Ramachandran / ms", true ))
    print("| $(rpad(label, 21)) |")
    for software in ("BioStructures", "MIToS", "Biopython", "ProDy", "MDAnalysis", "Bio3D",
                        "Rpdb", "BioPerl", "BioRuby", "Victor", "ESBTL")
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
