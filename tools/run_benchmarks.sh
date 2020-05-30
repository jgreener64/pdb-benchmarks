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
run_benchmark $nb "julia BioStructures/parse_pdb.jl data/1CRN.pdb"   "BioStructures/parse_pdb_1CRN.dat"   "BioStructures,Parse PDB 1CRN,"
run_benchmark $ns "julia BioStructures/parse_pdb.jl data/1HTQ.pdb"   "BioStructures/parse_pdb_1HTQ.dat"   "BioStructures,Parse PDB 1HTQ,"
run_benchmark $nb "julia BioStructures/parse_mmcif.jl data/1CRN.cif" "BioStructures/parse_mmcif_1CRN.dat" "BioStructures,Parse mmCIF 1CRN,"
run_benchmark $ns "julia BioStructures/parse_mmcif.jl data/1HTQ.cif" "BioStructures/parse_mmcif_1HTQ.dat" "BioStructures,Parse mmCIF 1HTQ,"
run_benchmark $nb "julia BioStructures/parse_mmtf.jl data/1CRN.mmtf" "BioStructures/parse_mmtf_1CRN.dat"  "BioStructures,Parse MMTF 1CRN,"
run_benchmark $ns "julia BioStructures/parse_mmtf.jl data/1HTQ.mmtf" "BioStructures/parse_mmtf_1HTQ.dat"  "BioStructures,Parse MMTF 1HTQ,"
run_benchmark $nb "julia BioStructures/count.jl"                     "BioStructures/count.dat"            "BioStructures,Count,"
run_benchmark $nb "julia BioStructures/distance.jl"                  "BioStructures/distance.dat"         "BioStructures,Distance,"
run_benchmark $nb "julia BioStructures/ramachandran.jl"              "BioStructures/ramachandran.dat"     "BioStructures,Ramachandran,"
echo "BioStructures benchmarks done"

# MIToS
run_benchmark $nb "julia MIToS/parse_pdb.jl data/1CRN.pdb" "MIToS/parse_pdb_1CRN.dat" "MIToS,Parse PDB 1CRN,"
run_benchmark $ns "julia MIToS/parse_pdb.jl data/1HTQ.pdb" "MIToS/parse_pdb_1HTQ.dat" "MIToS,Parse PDB 1HTQ,"
run_benchmark $nb "julia MIToS/count.jl"                   "MIToS/count.dat"          "MIToS,Count,"
run_benchmark $nb "julia MIToS/distance.jl"                "MIToS/distance.dat"       "MIToS,Distance,"
echo "MIToS benchmarks done"

# Biopython
run_benchmark $nb "python Biopython/parse_pdb.py data/1CRN.pdb"   "Biopython/parse_pdb_1CRN.dat"   "Biopython,Parse PDB 1CRN,"
run_benchmark $ns "python Biopython/parse_pdb.py data/1HTQ.pdb"   "Biopython/parse_pdb_1HTQ.dat"   "Biopython,Parse PDB 1HTQ,"
run_benchmark $nb "python Biopython/parse_mmcif.py data/1CRN.cif" "Biopython/parse_mmcif_1CRN.dat" "Biopython,Parse mmCIF 1CRN,"
run_benchmark $ns "python Biopython/parse_mmcif.py data/1HTQ.cif" "Biopython/parse_mmcif_1HTQ.dat" "Biopython,Parse mmCIF 1HTQ,"
run_benchmark $nb "python Biopython/parse_mmtf.py data/1CRN.mmtf" "Biopython/parse_mmtf_1CRN.dat"  "Biopython,Parse MMTF 1CRN,"
run_benchmark $ns "python Biopython/parse_mmtf.py data/1HTQ.mmtf" "Biopython/parse_mmtf_1HTQ.dat"  "Biopython,Parse MMTF 1HTQ,"
run_benchmark $nb "python Biopython/count.py"                     "Biopython/count.dat"            "Biopython,Count,"
run_benchmark $nb "python Biopython/distance.py"                  "Biopython/distance.dat"         "Biopython,Distance,"
run_benchmark $nb "python Biopython/ramachandran.py"              "Biopython/ramachandran.dat"     "Biopython,Ramachandran,"
echo "Biopython benchmarks done"

# ProDy
run_benchmark $nb "python ProDy/parse_pdb.py data/1CRN.pdb" "ProDy/parse_pdb_1CRN.dat" "ProDy,Parse PDB 1CRN,"
run_benchmark $ns "python ProDy/parse_pdb.py data/1HTQ.pdb" "ProDy/parse_pdb_1HTQ.dat" "ProDy,Parse PDB 1HTQ,"
run_benchmark $nb "python ProDy/count.py"                   "ProDy/count.dat"          "ProDy,Count,"
run_benchmark $nb "python ProDy/distance.py"                "ProDy/distance.dat"       "ProDy,Distance,"
run_benchmark $nb "python ProDy/ramachandran.py"            "ProDy/ramachandran.dat"   "ProDy,Ramachandran,"
echo "ProDy benchmarks done"

# MDAnalysis
run_benchmark $nb "python MDAnalysis/parse_pdb.py data/1CRN.pdb" "MDAnalysis/parse_pdb_1CRN.dat" "MDAnalysis,Parse PDB 1CRN,"
run_benchmark $ns "python MDAnalysis/parse_pdb.py data/1HTQ.pdb" "MDAnalysis/parse_pdb_1HTQ.dat" "MDAnalysis,Parse PDB 1HTQ,"
run_benchmark $nb "python MDAnalysis/count.py"                   "MDAnalysis/count.dat"          "MDAnalysis,Count,"
run_benchmark $nb "python MDAnalysis/distance.py"                "MDAnalysis/distance.dat"       "MDAnalysis,Distance,"
run_benchmark $nb "python MDAnalysis/ramachandran.py"            "MDAnalysis/ramachandran.dat"   "MDAnalysis,Ramachandran,"
echo "MDAnalysis benchmarks done"

# biotite
run_benchmark $nb "python biotite/parse_pdb.py data/1CRN.pdb"   "biotite/parse_pdb_1CRN.dat"   "biotite,Parse PDB 1CRN,"
run_benchmark $ns "python biotite/parse_pdb.py data/1HTQ.pdb"   "biotite/parse_pdb_1HTQ.dat"   "biotite,Parse PDB 1HTQ,"
run_benchmark $nb "python biotite/parse_mmcif.py data/1CRN.cif" "biotite/parse_mmcif_1CRN.dat" "biotite,Parse mmCIF 1CRN,"
run_benchmark $ns "python biotite/parse_mmcif.py data/1HTQ.cif" "biotite/parse_mmcif_1HTQ.dat" "biotite,Parse mmCIF 1HTQ,"
run_benchmark $nb "python biotite/parse_mmtf.py data/1CRN.mmtf" "biotite/parse_mmtf_1CRN.dat"  "biotite,Parse MMTF 1CRN,"
run_benchmark $ns "python biotite/parse_mmtf.py data/1HTQ.mmtf" "biotite/parse_mmtf_1HTQ.dat"  "biotite,Parse MMTF 1HTQ,"
echo "biotite benchmarks done"

# atomium
run_benchmark $nb "python atomium/parse_pdb.py   data/1CRN.pdb" "atomium/parse_pdb_1CRN.dat"   "atomium,Parse PDB 1CRN,"
run_benchmark $ns "python atomium/parse_pdb.py   data/1HTQ.pdb" "atomium/parse_pdb_1HTQ.dat"   "atomium,Parse PDB 1HTQ,"
run_benchmark $nb "python atomium/parse_mmcif.py data/1CRN.cif" "atomium/parse_mmcif_1CRN.dat" "atomium,Parse mmCIF 1CRN,"
run_benchmark $ns "python atomium/parse_mmcif.py data/1HTQ.cif" "atomium/parse_mmcif_1HTQ.dat" "atomium,Parse mmCIF 1HTQ,"
run_benchmark $nb "python atomium/parse_mmtf.py data/1CRN.mmtf" "atomium/parse_mmtf_1CRN.dat"  "atomium,Parse MMTF 1CRN,"
run_benchmark $ns "python atomium/parse_mmtf.py data/1HTQ.mmtf" "atomium/parse_mmtf_1HTQ.dat"  "atomium,Parse MMTF 1HTQ,"
echo "atomium benchmarks done"

# Bio3D
run_benchmark $nb "Rscript Bio3D/parse_pdb.R data/1CRN.pdb" "Bio3D/parse_pdb_1CRN.dat" "Bio3D,Parse PDB 1CRN,"
run_benchmark $ns "Rscript Bio3D/parse_pdb.R data/1HTQ.pdb" "Bio3D/parse_pdb_1HTQ.dat" "Bio3D,Parse PDB 1HTQ,"
run_benchmark $nb "Rscript Bio3D/count.R"                   "Bio3D/count.dat"          "Bio3D,Count,"
run_benchmark $nb "Rscript Bio3D/distance.R"                "Bio3D/distance.dat"       "Bio3D,Distance,"
echo "Bio3D benchmarks done"

# Rpdb
run_benchmark $nb "Rscript Rpdb/parse_pdb.R data/1CRN.pdb" "Rpdb/parse_pdb_1CRN.dat" "Rpdb,Parse PDB 1CRN,"
run_benchmark $ns "Rscript Rpdb/parse_pdb.R data/1HTQ.pdb" "Rpdb/parse_pdb_1HTQ.dat" "Rpdb,Parse PDB 1HTQ,"
run_benchmark $nb "Rscript Rpdb/count.R"                   "Rpdb/count.dat"          "Rpdb,Count,"
run_benchmark $nb "Rscript Rpdb/distance.R"                "Rpdb/distance.dat"       "Rpdb,Distance,"
echo "Rpdb benchmarks done"

# BioJava
run_benchmark $nb "java -cp BioJava/target/pdb-benchmarks-1.0-SNAPSHOT.jar com.jgreener.pdb.parse_pdb data/1CRN.pdb"   "BioJava/parse_pdb_1CRN.dat"   "BioJava,Parse PDB 1CRN,"
run_benchmark $ns "java -cp BioJava/target/pdb-benchmarks-1.0-SNAPSHOT.jar com.jgreener.pdb.parse_pdb data/1HTQ.pdb"   "BioJava/parse_pdb_1HTQ.dat"   "BioJava,Parse PDB 1HTQ,"
run_benchmark $nb "java -cp BioJava/target/pdb-benchmarks-1.0-SNAPSHOT.jar com.jgreener.pdb.parse_mmcif data/1CRN.cif" "BioJava/parse_mmcif_1CRN.dat" "BioJava,Parse mmCIF 1CRN,"
run_benchmark $ns "java -cp BioJava/target/pdb-benchmarks-1.0-SNAPSHOT.jar com.jgreener.pdb.parse_mmcif data/1HTQ.cif" "BioJava/parse_mmcif_1HTQ.dat" "BioJava,Parse mmCIF 1HTQ,"
run_benchmark $nb "java -cp BioJava/target/pdb-benchmarks-1.0-SNAPSHOT.jar com.jgreener.pdb.parse_mmtf data/1CRN.mmtf" "BioJava/parse_mmtf_1CRN.dat"  "BioJava,Parse MMTF 1CRN,"
run_benchmark $ns "java -cp BioJava/target/pdb-benchmarks-1.0-SNAPSHOT.jar com.jgreener.pdb.parse_mmtf data/1HTQ.mmtf" "BioJava/parse_mmtf_1HTQ.dat"  "BioJava,Parse MMTF 1HTQ,"
echo "BioJava benchmarks done"

# BioPerl
run_benchmark $nb "perl BioPerl/parse_pdb.pl data/1CRN.pdb" "BioPerl/parse_pdb_1CRN.dat" "BioPerl,Parse PDB 1CRN,"
run_benchmark $ns "perl BioPerl/parse_pdb.pl data/1HTQ.pdb" "BioPerl/parse_pdb_1HTQ.dat" "BioPerl,Parse PDB 1HTQ,"
run_benchmark $nb "perl BioPerl/count.pl"                   "BioPerl/count.dat"          "BioPerl,Count,"
run_benchmark $nb "perl BioPerl/distance.pl"                "BioPerl/distance.dat"       "BioPerl,Distance,"
echo "BioPerl benchmarks done"

# BioRuby
run_benchmark $nb "ruby BioRuby/parse_pdb.rb data/1CRN.pdb" "BioRuby/parse_pdb_1CRN.dat" "BioRuby,Parse PDB 1CRN,"
run_benchmark $ns "ruby BioRuby/parse_pdb.rb data/1HTQ.pdb" "BioRuby/parse_pdb_1HTQ.dat" "BioRuby,Parse PDB 1HTQ,"
run_benchmark $nb "ruby BioRuby/count.rb"                   "BioRuby/count.dat"          "BioRuby,Count,"
run_benchmark $nb "ruby BioRuby/distance.rb"                "BioRuby/distance.dat"       "BioRuby,Distance,"
echo "BioRuby benchmarks done"

# GEMMI
run_benchmark $nb "GEMMI/parse_pdb data/1CRN.pdb"   "GEMMI/parse_pdb_1CRN.dat"   "GEMMI,Parse PDB 1CRN,"
run_benchmark $ns "GEMMI/parse_pdb data/1HTQ.pdb"   "GEMMI/parse_pdb_1HTQ.dat"   "GEMMI,Parse PDB 1HTQ,"
run_benchmark $nb "GEMMI/parse_mmcif data/1CRN.cif" "GEMMI/parse_mmcif_1CRN.dat" "GEMMI,Parse mmCIF 1CRN,"
run_benchmark $ns "GEMMI/parse_mmcif data/1HTQ.cif" "GEMMI/parse_mmcif_1HTQ.dat" "GEMMI,Parse mmCIF 1HTQ,"
echo "GEMMI benchmarks done"

# Victor
run_benchmark $nb "Victor/parse_pdb data/1CRN.pdb" "Victor/parse_pdb_1CRN.dat" "Victor,Parse PDB 1CRN,"
run_benchmark $ns "Victor/parse_pdb data/1HTQ.pdb" "Victor/parse_pdb_1HTQ.dat" "Victor,Parse PDB 1HTQ,"
echo "Victor benchmarks done"

# ESBTL
run_benchmark $nb "ESBTL/parse_pdb data/1CRN.pdb" "ESBTL/parse_pdb_1CRN.dat" "ESBTL,Parse PDB 1CRN,"
echo "ESBTL benchmarks done"

# chemfiles - Python
run_benchmark $nb "python chemfiles/parse_pdb.py data/1CRN.pdb"   "chemfiles/parse_pdb_1CRN_py.dat"   "chemfiles-python,Parse PDB 1CRN,"
# FIXME: this uncovered a bug in chemfiles, the bugfix will be avaible on
# chemfiles>=0.10 when released
#run_benchmark $ns "python chemfiles/parse_pdb.py data/1HTQ.pdb"   "chemfiles/parse_pdb_1HTQ_py.dat"   "chemfiles-python,Parse PDB 1HTQ,"
run_benchmark $nb "python chemfiles/parse_mmcif.py data/1CRN.cif" "chemfiles/parse_mmcif_1CRN_py.dat" "chemfiles-python,Parse mmCIF 1CRN,"
run_benchmark $ns "python chemfiles/parse_mmcif.py data/1HTQ.cif" "chemfiles/parse_mmcif_1HTQ_py.dat" "chemfiles-python,Parse mmCIF 1HTQ,"
run_benchmark $nb "python chemfiles/parse_mmtf.py data/1CRN.mmtf" "chemfiles/parse_mmtf_1CRN_py.dat"  "chemfiles-python,Parse MMTF 1CRN,"
#run_benchmark $ns "python chemfiles/parse_mmtf.py data/1HTQ.mmtf" "chemfiles/parse_mmtf_1HTQ_py.dat"  "chemfiles-python,Parse MMTF 1HTQ,"
run_benchmark $nb "python chemfiles/count.py"                     "chemfiles/count_py.dat"            "chemfiles-python,Count,"
run_benchmark $nb "python chemfiles/distance.py"                  "chemfiles/distance_py.dat"         "chemfiles-python,Distance,"
run_benchmark $nb "python chemfiles/ramachandran.py"              "chemfiles/ramachandran_py.dat"     "chemfiles-python,Ramachandran,"
echo "chemfiles-python benchmarks done"

# chemfiles - C++
run_benchmark $nb "chemfiles/parse_pdb data/1CRN.pdb"   "chemfiles/parse_pdb_1CRN_cxx.dat"   "chemfiles-cxx,Parse PDB 1CRN,"
#run_benchmark $ns "chemfiles/parse_pdb data/1HTQ.pdb"   "chemfiles/parse_pdb_1HTQ_cxx.dat"   "chemfiles-cxx,Parse PDB 1HTQ,"
run_benchmark $nb "chemfiles/parse_mmcif data/1CRN.cif" "chemfiles/parse_mmcif_1CRN_cxx.dat" "chemfiles-cxx,Parse mmCIF 1CRN,"
run_benchmark $ns "chemfiles/parse_mmcif data/1HTQ.cif" "chemfiles/parse_mmcif_1HTQ_cxx.dat" "chemfiles-cxx,Parse mmCIF 1HTQ,"
run_benchmark $nb "chemfiles/parse_mmtf data/1CRN.mmtf" "chemfiles/parse_mmtf_1CRN_cxx.dat"  "chemfiles-cxx,Parse MMTF 1CRN,"
#run_benchmark $ns "chemfiles/parse_mmtf data/1HTQ.mmtf" "chemfiles/parse_mmtf_1HTQ_cxx.dat"  "chemfiles-cxx,Parse MMTF 1HTQ,"
run_benchmark $nb "chemfiles/count"                     "chemfiles/count_cxx.dat"            "chemfiles-cxx,Count,"
run_benchmark $nb "chemfiles/distance"                  "chemfiles/distance_cxx.dat"         "chemfiles-cxx,Distance,"
run_benchmark $nb "chemfiles/ramachandran"              "chemfiles/ramachandran_cxx.dat"     "chemfiles-cxx,Ramachandran,"
echo "chemfiles-cxx benchmarks done"

# Plot results
julia plot/plot.jl
echo "Results plotted"
