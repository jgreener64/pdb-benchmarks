# Run all benchmarks and save results as .dat files

# Number of runs
n=10
# Number of runs for 1HTQ parsing
ns=2

# Remove current data files
rm */*.dat

# BioJulia
for i in $(seq 1 $n); do julia BioJulia/parse.jl pdbs/1CRN.pdb >> BioJulia/parse_1CRN.dat; done
for i in $(seq 1 $n); do julia BioJulia/parse.jl pdbs/3JYV.pdb >> BioJulia/parse_3JYV.dat; done
for i in $(seq 1 $ns); do julia BioJulia/parse.jl pdbs/1HTQ.pdb >> BioJulia/parse_1HTQ.dat; done
for i in $(seq 1 $n); do julia BioJulia/count.jl >> BioJulia/count.dat; done
for i in $(seq 1 $n); do julia BioJulia/distance.jl >> BioJulia/distance.dat; done
for i in $(seq 1 $n); do julia BioJulia/ramachandran.jl >> BioJulia/ramachandran.dat; done

# Biopython
for i in $(seq 1 $n); do python Biopython/parse.py pdbs/1CRN.pdb >> Biopython/parse_1CRN.dat; done
for i in $(seq 1 $n); do python Biopython/parse.py pdbs/3JYV.pdb >> Biopython/parse_3JYV.dat; done
for i in $(seq 1 $ns); do python Biopython/parse.py pdbs/1HTQ.pdb >> Biopython/parse_1HTQ.dat; done
for i in $(seq 1 $n); do python Biopython/count.py >> Biopython/count.dat; done
for i in $(seq 1 $n); do python Biopython/distance.py >> Biopython/distance.dat; done
for i in $(seq 1 $n); do python Biopython/ramachandran.py >> Biopython/ramachandran.dat; done

# ProDy
for i in $(seq 1 $n); do python ProDy/parse.py pdbs/1CRN.pdb >> ProDy/parse_1CRN.dat; done
for i in $(seq 1 $n); do python ProDy/parse.py pdbs/3JYV.pdb >> ProDy/parse_3JYV.dat; done
for i in $(seq 1 $ns); do python ProDy/parse.py pdbs/1HTQ.pdb >> ProDy/parse_1HTQ.dat; done
for i in $(seq 1 $n); do python ProDy/count.py >> ProDy/count.dat; done
for i in $(seq 1 $n); do python ProDy/distance.py >> ProDy/distance.dat; done
for i in $(seq 1 $n); do python ProDy/ramachandran.py >> ProDy/ramachandran.dat; done

# Bio3D
for i in $(seq 1 $n); do Rscript Bio3D/parse.R pdbs/1CRN.pdb | tail -n1 >> Bio3D/parse_1CRN.dat; done
for i in $(seq 1 $n); do Rscript Bio3D/parse.R pdbs/3JYV.pdb | tail -n1 >> Bio3D/parse_3JYV.dat; done
for i in $(seq 1 $ns); do Rscript Bio3D/parse.R pdbs/1HTQ.pdb | tail -n1 >> Bio3D/parse_1HTQ.dat; done
for i in $(seq 1 $n); do Rscript Bio3D/count.R | tail -n1 >> Bio3D/count.dat; done
for i in $(seq 1 $n); do Rscript Bio3D/distance.R | tail -n1 >> Bio3D/distance.dat; done

# Rpdb
for i in $(seq 1 $n); do Rscript Rpdb/parse.R pdbs/1CRN.pdb >> Rpdb/parse_1CRN.dat; done
for i in $(seq 1 $n); do Rscript Rpdb/parse.R pdbs/3JYV.pdb >> Rpdb/parse_3JYV.dat; done
for i in $(seq 1 $ns); do Rscript Rpdb/parse.R pdbs/1HTQ.pdb >> Rpdb/parse_1HTQ.dat; done
for i in $(seq 1 $n); do Rscript Rpdb/count.R >> Rpdb/count.dat; done
for i in $(seq 1 $n); do Rscript Rpdb/distance.R >> Rpdb/distance.dat; done

# BioPerl
for i in $(seq 1 $n); do perl BioPerl/parse.pl pdbs/1CRN.pdb >> BioPerl/parse_1CRN.dat; done
for i in $(seq 1 $n); do perl BioPerl/parse.pl pdbs/3JYV.pdb >> BioPerl/parse_3JYV.dat; done
for i in $(seq 1 $ns); do perl BioPerl/parse.pl pdbs/1HTQ.pdb >> BioPerl/parse_1HTQ.dat; done
for i in $(seq 1 $n); do perl BioPerl/count.pl >> BioPerl/count.dat; done
for i in $(seq 1 $n); do perl BioPerl/distance.pl >> BioPerl/distance.dat; done

# BioRuby
for i in $(seq 1 $n); do ruby BioRuby/parse.rb pdbs/1CRN.pdb >> BioRuby/parse_1CRN.dat; done
for i in $(seq 1 $n); do ruby BioRuby/parse.rb pdbs/3JYV.pdb >> BioRuby/parse_3JYV.dat; done
for i in $(seq 1 $ns); do ruby BioRuby/parse.rb pdbs/1HTQ.pdb >> BioRuby/parse_1HTQ.dat; done
for i in $(seq 1 $n); do ruby BioRuby/count.rb >> BioRuby/count.dat; done
for i in $(seq 1 $n); do ruby BioRuby/distance.rb >> BioRuby/distance.dat; done
