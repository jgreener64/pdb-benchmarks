# Run all benchmarks, save the results and form a csv file for plotting
# Requires all packages to be installed and compiled where applicable

# Number of runs for each benchmark apart from 1HTQ parsing
nb=10
# Number of runs for 1HTQ parsing
ns=3

# Remove current data files
rm */*.dat

# Remove current plot
rm plot/plot.png

# Reset benchmarking results
echo "package,benchmark,time" > benchmarks.csv

# Run a benchmark
# Arguments are number of runs, benchmark command, output data file, csv file columns
function run_benchmark {
    for i in $(seq 1 $1)
    do
        eval $2 >> $3
    done
    echo -n $4 >> benchmarks.csv
    python tools/mean.py $3 >> benchmarks.csv
}

# BioJulia
run_benchmark $nb "julia BioJulia/parse.jl pdbs/1CRN.pdb"          "BioJulia/parse_1CRN.dat"    "BioJulia,parse 1CRN,"
run_benchmark $nb "julia BioJulia/parse.jl pdbs/3JYV.pdb"          "BioJulia/parse_3JYV.dat"    "BioJulia,parse 3JYV,"
run_benchmark $ns "julia BioJulia/parse.jl pdbs/1HTQ.pdb"          "BioJulia/parse_1HTQ.dat"    "BioJulia,parse 1HTQ,"
run_benchmark $nb "julia BioJulia/count.jl"                        "BioJulia/count.dat"         "BioJulia,count,"
run_benchmark $nb "julia BioJulia/distance.jl"                     "BioJulia/distance.dat"      "BioJulia,distance,"
run_benchmark $nb "julia BioJulia/ramachandran.jl"                 "BioJulia/ramachandran.dat"  "BioJulia,ramachandran,"
echo "BioJulia benchmarks done"

# MIToS
run_benchmark $nb "julia MIToS/parse.jl pdbs/1CRN.pdb"             "MIToS/parse_1CRN.dat"       "MIToS,parse 1CRN,"
run_benchmark $nb "julia MIToS/parse.jl pdbs/3JYV.pdb"             "MIToS/parse_3JYV.dat"       "MIToS,parse 3JYV,"
run_benchmark $ns "julia MIToS/parse.jl pdbs/1HTQ.pdb"             "MIToS/parse_1HTQ.dat"       "MIToS,parse 1HTQ,"
run_benchmark $nb "julia MIToS/count.jl"                           "MIToS/count.dat"            "MIToS,count,"
run_benchmark $nb "julia MIToS/distance.jl"                        "MIToS/distance.dat"         "MIToS,distance,"
echo "MIToS benchmarks done"

# Biopython
run_benchmark $nb "python Biopython/parse.py pdbs/1CRN.pdb"        "Biopython/parse_1CRN.dat"   "Biopython,parse 1CRN,"
run_benchmark $nb "python Biopython/parse.py pdbs/3JYV.pdb"        "Biopython/parse_3JYV.dat"   "Biopython,parse 3JYV,"
run_benchmark $ns "python Biopython/parse.py pdbs/1HTQ.pdb"        "Biopython/parse_1HTQ.dat"   "Biopython,parse 1HTQ,"
run_benchmark $nb "python Biopython/count.py"                      "Biopython/count.dat"        "Biopython,count,"
run_benchmark $nb "python Biopython/distance.py"                   "Biopython/distance.dat"     "Biopython,distance,"
run_benchmark $nb "python Biopython/ramachandran.py"               "Biopython/ramachandran.dat" "Biopython,ramachandran,"
echo "Biopython benchmarks done"

# ProDy
run_benchmark $nb "python ProDy/parse.py pdbs/1CRN.pdb"            "ProDy/parse_1CRN.dat"       "ProDy,parse 1CRN,"
run_benchmark $nb "python ProDy/parse.py pdbs/3JYV.pdb"            "ProDy/parse_3JYV.dat"       "ProDy,parse 3JYV,"
run_benchmark $ns "python ProDy/parse.py pdbs/1HTQ.pdb"            "ProDy/parse_1HTQ.dat"       "ProDy,parse 1HTQ,"
run_benchmark $nb "python ProDy/count.py"                          "ProDy/count.dat"            "ProDy,count,"
run_benchmark $nb "python ProDy/distance.py"                       "ProDy/distance.dat"         "ProDy,distance,"
run_benchmark $nb "python ProDy/ramachandran.py"                   "ProDy/ramachandran.dat"     "ProDy,ramachandran,"
echo "ProDy benchmarks done"

# Bio3D
run_benchmark $nb "Rscript Bio3D/parse.R pdbs/1CRN.pdb | tail -n1" "Bio3D/parse_1CRN.dat"       "Bio3D,parse 1CRN,"
run_benchmark $nb "Rscript Bio3D/parse.R pdbs/3JYV.pdb | tail -n1" "Bio3D/parse_3JYV.dat"       "Bio3D,parse 3JYV,"
run_benchmark $ns "Rscript Bio3D/parse.R pdbs/1HTQ.pdb | tail -n1" "Bio3D/parse_1HTQ.dat"       "Bio3D,parse 1HTQ,"
run_benchmark $nb "Rscript Bio3D/count.R | tail -n1"               "Bio3D/count.dat"            "Bio3D,count,"
run_benchmark $nb "Rscript Bio3D/distance.R | tail -n1"            "Bio3D/distance.dat"         "Bio3D,distance,"
echo "Bio3D benchmarks done"

# Rpdb
run_benchmark $nb "Rscript Rpdb/parse.R pdbs/1CRN.pdb"             "Rpdb/parse_1CRN.dat"        "Rpdb,parse 1CRN,"
run_benchmark $nb "Rscript Rpdb/parse.R pdbs/3JYV.pdb"             "Rpdb/parse_3JYV.dat"        "Rpdb,parse 3JYV,"
run_benchmark $ns "Rscript Rpdb/parse.R pdbs/1HTQ.pdb"             "Rpdb/parse_1HTQ.dat"        "Rpdb,parse 1HTQ,"
run_benchmark $nb "Rscript Rpdb/count.R"                           "Rpdb/count.dat"             "Rpdb,count,"
run_benchmark $nb "Rscript Rpdb/distance.R"                        "Rpdb/distance.dat"          "Rpdb,distance,"
echo "Rpdb benchmarks done"

# BioPerl
run_benchmark $nb "perl BioPerl/parse.pl pdbs/1CRN.pdb"            "BioPerl/parse_1CRN.dat"     "BioPerl,parse 1CRN,"
run_benchmark $nb "perl BioPerl/parse.pl pdbs/3JYV.pdb"            "BioPerl/parse_3JYV.dat"     "BioPerl,parse 3JYV,"
run_benchmark $ns "perl BioPerl/parse.pl pdbs/1HTQ.pdb"            "BioPerl/parse_1HTQ.dat"     "BioPerl,parse 1HTQ,"
run_benchmark $nb "perl BioPerl/count.pl"                          "BioPerl/count.dat"          "BioPerl,count,"
run_benchmark $nb "perl BioPerl/distance.pl"                       "BioPerl/distance.dat"       "BioPerl,distance,"
echo "BioPerl benchmarks done"

# BioRuby
run_benchmark $nb "ruby BioRuby/parse.rb pdbs/1CRN.pdb"            "BioRuby/parse_1CRN.dat"     "BioRuby,parse 1CRN,"
run_benchmark $nb "ruby BioRuby/parse.rb pdbs/3JYV.pdb"            "BioRuby/parse_3JYV.dat"     "BioRuby,parse 3JYV,"
run_benchmark $ns "ruby BioRuby/parse.rb pdbs/1HTQ.pdb"            "BioRuby/parse_1HTQ.dat"     "BioRuby,parse 1HTQ,"
run_benchmark $nb "ruby BioRuby/count.rb"                          "BioRuby/count.dat"          "BioRuby,count,"
run_benchmark $nb "ruby BioRuby/distance.rb"                       "BioRuby/distance.dat"       "BioRuby,distance,"
echo "BioRuby benchmarks done"

# Victor
run_benchmark $nb "Victor/parse pdbs/1CRN.pdb"                     "Victor/parse_1CRN.dat"      "Victor,parse 1CRN,"
run_benchmark $nb "Victor/parse pdbs/3JYV.pdb"                     "Victor/parse_3JYV.dat"      "Victor,parse 3JYV,"
run_benchmark $ns "Victor/parse pdbs/1HTQ.pdb"                     "Victor/parse_1HTQ.dat"      "Victor,parse 1HTQ,"
echo "Victor benchmarks done"

# ESBTL
run_benchmark $nb "ESBTL/parse pdbs/1CRN.pdb"                      "ESBTL/parse_1CRN.dat"       "ESBTL,parse 1CRN,"
run_benchmark $nb "ESBTL/parse pdbs/3JYV.pdb"                      "ESBTL/parse_3JYV.dat"       "ESBTL,parse 3JYV,"
echo "ESBTL benchmarks done"

# Plot results
julia plot/plot.jl
echo "Results plotted"
