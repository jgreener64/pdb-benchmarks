# Run all benchmarks, save the results and form a csv file for plotting
# Requires all packages to be installed and compiled where applicable

echo "Running benchmarks"

# Number of runs for each benchmark apart from 1HTQ parsing
nb=10
# Number of runs for 1HTQ parsing
ns=3

# Remove current data files
rm */*.dat

# Remove current plot
rm plot/plot.png

# Reset benchmarking results
echo "Package,Benchmark,Runtime" > benchmarks.csv

# Run a benchmark
# Arguments are number of runs, benchmark command, output data file, csv file columns
function run_benchmark {
    for i in $(seq 1 $1)
    do
        eval $2 | tail -n1 >> $3
    done
    echo -n $4 >> benchmarks.csv
    python tools/mean.py $3 >> benchmarks.csv
}

# BioStructures
run_benchmark $nb "julia BioStructures/parse.jl pdbs/1CRN.pdb"    "BioStructures/parse_1CRN.dat"     "BioStructures,Parse 1CRN,"
run_benchmark $nb "julia BioStructures/parse.jl pdbs/3JYV.pdb"    "BioStructures/parse_3JYV.dat"     "BioStructures,Parse 3JYV,"
run_benchmark $ns "julia BioStructures/parse.jl pdbs/1HTQ.pdb"    "BioStructures/parse_1HTQ.dat"     "BioStructures,Parse 1HTQ,"
run_benchmark $nb "julia BioStructures/count.jl"                  "BioStructures/count.dat"          "BioStructures,Count,"
run_benchmark $nb "julia BioStructures/distance.jl"               "BioStructures/distance.dat"       "BioStructures,Distance,"
run_benchmark $nb "julia BioStructures/ramachandran.jl"           "BioStructures/ramachandran.dat"   "BioStructures,Ramachandran,"
echo "BioStructures benchmarks done"

# MIToS
run_benchmark $nb "julia MIToS/parse.jl pdbs/1CRN.pdb"      "MIToS/parse_1CRN.dat"        "MIToS,Parse 1CRN,"
run_benchmark $nb "julia MIToS/parse.jl pdbs/3JYV.pdb"      "MIToS/parse_3JYV.dat"        "MIToS,Parse 3JYV,"
run_benchmark $ns "julia MIToS/parse.jl pdbs/1HTQ.pdb"      "MIToS/parse_1HTQ.dat"        "MIToS,Parse 1HTQ,"
run_benchmark $nb "julia MIToS/count.jl"                    "MIToS/count.dat"             "MIToS,Count,"
run_benchmark $nb "julia MIToS/distance.jl"                 "MIToS/distance.dat"          "MIToS,Distance,"
echo "MIToS benchmarks done"

# Biopython
run_benchmark $nb "python Biopython/parse.py pdbs/1CRN.pdb"  "Biopython/parse_1CRN.dat"    "Biopython,Parse 1CRN,"
run_benchmark $nb "python Biopython/parse.py pdbs/3JYV.pdb"  "Biopython/parse_3JYV.dat"    "Biopython,Parse 3JYV,"
run_benchmark $ns "python Biopython/parse.py pdbs/1HTQ.pdb"  "Biopython/parse_1HTQ.dat"    "Biopython,Parse 1HTQ,"
run_benchmark $nb "python Biopython/count.py"                "Biopython/count.dat"         "Biopython,Count,"
run_benchmark $nb "python Biopython/distance.py"             "Biopython/distance.dat"      "Biopython,Distance,"
run_benchmark $nb "python Biopython/ramachandran.py"         "Biopython/ramachandran.dat"  "Biopython,Ramachandran,"
echo "Biopython benchmarks done"

# ProDy
run_benchmark $nb "python ProDy/parse.py pdbs/1CRN.pdb"      "ProDy/parse_1CRN.dat"        "ProDy,Parse 1CRN,"
run_benchmark $nb "python ProDy/parse.py pdbs/3JYV.pdb"      "ProDy/parse_3JYV.dat"        "ProDy,Parse 3JYV,"
run_benchmark $ns "python ProDy/parse.py pdbs/1HTQ.pdb"      "ProDy/parse_1HTQ.dat"        "ProDy,Parse 1HTQ,"
run_benchmark $nb "python ProDy/count.py"                    "ProDy/count.dat"             "ProDy,Count,"
run_benchmark $nb "python ProDy/distance.py"                 "ProDy/distance.dat"          "ProDy,Distance,"
run_benchmark $nb "python ProDy/ramachandran.py"             "ProDy/ramachandran.dat"      "ProDy,Ramachandran,"
echo "ProDy benchmarks done"

# MDAnalysis
run_benchmark $nb "python MDAnalysis/parse.py pdbs/1CRN.pdb" "MDAnalysis/parse_1CRN.dat"   "MDAnalysis,Parse 1CRN,"
run_benchmark $nb "python MDAnalysis/parse.py pdbs/3JYV.pdb" "MDAnalysis/parse_3JYV.dat"   "MDAnalysis,Parse 3JYV,"
run_benchmark $ns "python MDAnalysis/parse.py pdbs/1HTQ.pdb" "MDAnalysis/parse_1HTQ.dat"   "MDAnalysis,Parse 1HTQ,"
run_benchmark $nb "python MDAnalysis/count.py"               "MDAnalysis/count.dat"        "MDAnalysis,Count,"
run_benchmark $nb "python MDAnalysis/distance.py"            "MDAnalysis/distance.dat"     "MDAnalysis,Distance,"
run_benchmark $nb "python MDAnalysis/ramachandran.py"        "MDAnalysis/ramachandran.dat" "MDAnalysis,Ramachandran,"
echo "MDAnalysis benchmarks done"

# Bio3D
run_benchmark $nb "Rscript Bio3D/parse.R pdbs/1CRN.pdb"      "Bio3D/parse_1CRN.dat"        "Bio3D,Parse 1CRN,"
run_benchmark $nb "Rscript Bio3D/parse.R pdbs/3JYV.pdb"      "Bio3D/parse_3JYV.dat"        "Bio3D,Parse 3JYV,"
run_benchmark $ns "Rscript Bio3D/parse.R pdbs/1HTQ.pdb"      "Bio3D/parse_1HTQ.dat"        "Bio3D,Parse 1HTQ,"
run_benchmark $nb "Rscript Bio3D/count.R"                    "Bio3D/count.dat"             "Bio3D,Count,"
run_benchmark $nb "Rscript Bio3D/distance.R"                 "Bio3D/distance.dat"          "Bio3D,Distance,"
echo "Bio3D benchmarks done"

# Rpdb
run_benchmark $nb "Rscript Rpdb/parse.R pdbs/1CRN.pdb"       "Rpdb/parse_1CRN.dat"         "Rpdb,Parse 1CRN,"
run_benchmark $nb "Rscript Rpdb/parse.R pdbs/3JYV.pdb"       "Rpdb/parse_3JYV.dat"         "Rpdb,Parse 3JYV,"
run_benchmark $ns "Rscript Rpdb/parse.R pdbs/1HTQ.pdb"       "Rpdb/parse_1HTQ.dat"         "Rpdb,Parse 1HTQ,"
run_benchmark $nb "Rscript Rpdb/count.R"                     "Rpdb/count.dat"              "Rpdb,Count,"
run_benchmark $nb "Rscript Rpdb/distance.R"                  "Rpdb/distance.dat"           "Rpdb,Distance,"
echo "Rpdb benchmarks done"

# BioPerl
run_benchmark $nb "perl BioPerl/parse.pl pdbs/1CRN.pdb"      "BioPerl/parse_1CRN.dat"      "BioPerl,Parse 1CRN,"
run_benchmark $nb "perl BioPerl/parse.pl pdbs/3JYV.pdb"      "BioPerl/parse_3JYV.dat"      "BioPerl,Parse 3JYV,"
run_benchmark $ns "perl BioPerl/parse.pl pdbs/1HTQ.pdb"      "BioPerl/parse_1HTQ.dat"      "BioPerl,Parse 1HTQ,"
run_benchmark $nb "perl BioPerl/count.pl"                    "BioPerl/count.dat"           "BioPerl,Count,"
run_benchmark $nb "perl BioPerl/distance.pl"                 "BioPerl/distance.dat"        "BioPerl,Distance,"
echo "BioPerl benchmarks done"

# BioRuby
run_benchmark $nb "ruby BioRuby/parse.rb pdbs/1CRN.pdb"      "BioRuby/parse_1CRN.dat"      "BioRuby,Parse 1CRN,"
run_benchmark $nb "ruby BioRuby/parse.rb pdbs/3JYV.pdb"      "BioRuby/parse_3JYV.dat"      "BioRuby,Parse 3JYV,"
run_benchmark $ns "ruby BioRuby/parse.rb pdbs/1HTQ.pdb"      "BioRuby/parse_1HTQ.dat"      "BioRuby,Parse 1HTQ,"
run_benchmark $nb "ruby BioRuby/count.rb"                    "BioRuby/count.dat"           "BioRuby,Count,"
run_benchmark $nb "ruby BioRuby/distance.rb"                 "BioRuby/distance.dat"        "BioRuby,Distance,"
echo "BioRuby benchmarks done"

# Victor
run_benchmark $nb "Victor/parse pdbs/1CRN.pdb"               "Victor/parse_1CRN.dat"       "Victor,Parse 1CRN,"
run_benchmark $nb "Victor/parse pdbs/3JYV.pdb"               "Victor/parse_3JYV.dat"       "Victor,Parse 3JYV,"
run_benchmark $ns "Victor/parse pdbs/1HTQ.pdb"               "Victor/parse_1HTQ.dat"       "Victor,Parse 1HTQ,"
echo "Victor benchmarks done"

# ESBTL
run_benchmark $nb "ESBTL/parse pdbs/1CRN.pdb"                "ESBTL/parse_1CRN.dat"        "ESBTL,Parse 1CRN,"
run_benchmark $nb "ESBTL/parse pdbs/3JYV.pdb"                "ESBTL/parse_3JYV.dat"        "ESBTL,Parse 3JYV,"
echo "ESBTL benchmarks done"

# Plot results
julia plot/plot.jl
echo "Results plotted"
