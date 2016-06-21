# PDB benchmarks

Open source software packages to parse [Protein Data Bank](http://www.rcsb.org/pdb/home/home.do) (PDB) files and manipulate protein structures exist in many languages, often as part of Bio* projects.

This repository aims to collate benchmarks for common tasks across various languages and packages. The collection of scripts may also be useful to get an idea how each package works.

Please feel free to contribute scripts from other packages, or submit improvements the scripts already present - I'm looking for the fastest implementation for each software that makes use of the provided API.

Disclosure: I contributed the `Bio.Structure` module to BioJulia.


## Tests

* Parsing 3 PDB files, taken from the benchmarking in [1]:
  * [1CRN](http://www.rcsb.org/pdb/explore/explore.do?structureId=1crn) - hydrophobic protein (327 atoms)
  * [3JYV](http://www.rcsb.org/pdb/explore/explore.do?structureId=3jyv) - 80S rRNA (57,327 atoms)
  * [1HTQ](http://www.rcsb.org/pdb/explore/explore.do?structureId=1htq) - multicopy glutamine synthetase (10 models of 97,872 atoms)
* Counting the number of alanine residues in adenylate kinase ([1AKE](http://www.rcsb.org/pdb/explore/explore.do?structureId=1ake))
* Calculating the distance between residues 50 and 60 of chain A in adenylate kinase ([1AKE](http://www.rcsb.org/pdb/explore/explore.do?structureId=1ake))
* Calculating the Ramachandran phi/psi angles in adenylate kinase ([1AKE](http://www.rcsb.org/pdb/explore/explore.do?structureId=1ake))

[1] Gajda MJ, hPDB - Haskell library for processing atomic biomolecular structures in protein data bank format, *BMC Research Notes* 2013, **6**:483 | [link](http://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-6-483)

The PDB files can be downloaded to directory `pdbs` by running `source tools/download_pdbs.sh` from this directory. If you have all the software installed, and compiled where applicable, you can use the script `tools/run_benchmarks.sh` from this directory to run the benchmarks. The mean over a number of runs is taken for each benchmark to obtain the values below.

Benchmarks were carried out on a 3.1 GHz Intel Core i7 processor with 16 GB 1867 MHz DDR3 RAM. The operating system was Mac OS X Yosemite 10.10.5. Time is the elapsed time.


## Software

* [BioJulia](https://biojulia.github.io/Bio.jl/) development version running on Julia v0.4.6 (times measured after JIT compilation)
* [Biopython](http://biopython.org/wiki/Biopython) v1.66 running on Python v2.7.10
* [ProDy](http://prody.csb.pitt.edu/) v1.7 running on Python v2.7.10
* [Bio3D](http://thegrantlab.org/bio3d/index.php) v2.2-2 running on R v3.2.2
* [Rpdb](https://cran.r-project.org/web/packages/Rpdb/index.html) v2.2 running on R v3.2.2
* [BioPerl](http://bioperl.org/index.html) v1.6.924 running on Perl v5.18.2
* [BioRuby](http://bioruby.org/) v1.5.0 running on Ruby v2.0.0
* [Victor](http://protein.bio.unipd.it/victor/index.php/Main_Page) v1.0 compiled with g++ v6.1.0
* [ESBTL](http://esbtl.sourceforge.net/index.html) v1.0-beta01 compiled with g++ v6.1.0


## Comparison

Note that direct comparison between these times should be treated with caution, as each package does something slightly different. For example, things that increase parsing time include:

* Parsing the PDB header
* Accounting for disorder at both the atom and residue (point mutation) level
* Forming a heirarchical model of the protein that makes access to specific residues, atoms etc. easier and faster after parsing
* Checking that the PDB format is adhered to at various levels of strictness

Each package supports these to varying degrees.

|                       | BioJulia     | Biopython    | ProDy        | Bio3D        | Rpdb         | BioPerl      | BioRuby      | Victor       | ESBTL        |
| :-------------------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- |
| Parse 1CRN / ms       | 3.0          | 10           | 2.4          | 33           | 18           | 62           | 27           | 11           | 6.2          |
| Parse 3JYV / s        | 0.63         | 1.0          | 0.30         | 14           | 2.1          | 3.8          | 1.0          | 8.1          | 0.96         |
| Parse 1HTQ / s        | 18           | 23           | 1.8          | 56           | 32           | 66           | 18           | 17           | -            |
| Count / ms            | 0.45         | 0.45         | 9.4          | 0.51         | 0.41         | 0.86         | 0.25         | -            | -            |
| Distance / ms         | 0.028        | 0.35         | 5.6          | 1.1          | 1.6          | 0.93         | 0.56         | -            | -            |
| Ramachandran / ms     | 8.8          | 150          | 210          | -            | -            | -            | -            | -            | -            |
| Language              | Julia        | Python       | Python       | R            | R            | Perl         | Ruby         | C++          | C++          |
| Parses header         | ✗            | ✓            | ✓            | ✓            | ✓            | ✗            | ✓            | ✓            | ✗            |
| Heirarchichal parsing | ✓            | ✓            | ✓            | ✗            | ✗            | ✓            | ✓            | ✓            | ✓            |
| Supports disorder     | ✓            | ✓            | ✗            | ✗            | ✗            | ✗            | ✗            | ✗            | ✓            |
| Writes PDBs           | ✓            | ✓            | ✓            | ✓            | ✓            | ✓            | ✗            | ✓            | ✓            |
| Superimposition       | ✗            | ✓            | ✓            | ✓            | ✗            | ✗            | ✗            | ✗            | ✗            |
| PCA                   | ✗            | ✗            | ✓            | ✓            | ✗            | ✗            | ✗            | ✗            | ✗            |
| License               | MIT          | Biopython    | MIT          | GPLv2        | GPL          | GPL/Artistic | Ruby         | GPLv3        | GPLv3        |

Benchmarks as a plot:

![benchmarks](plot/plot.png "benchmarks")


## Opinions

* For most purposes, particularly work on small numbers of files, the speed of the programs will not hold you back. In this case use the language/package you are most familiar with.
* If you are analysing ensembles of proteins use packages with that functionality, such as ProDy or Bio3D, rather than writing the code yourself.


## Contributing

If you want to contribute benchmarks for a package, please make a pull request with the script(s) in a directory like the other packages. I will run the benchmarks again and change the README, thanks.


## Plans

* Finish Ramachandran scripts for remaining languages.
* Test BioJava, hPDB, possibly others.
* Run methods on the whole of the PDB to look at how they deal with errors.


## Resources

* A list of PDB parsing packages, particularly in C/C++, can be found [here](http://bioinf.org.uk/software/bioplib/libraries/).
* The Biopython [documentation](http://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ) has a useful discussion on disorder at the atom and residue level.
* Sets of utility scripts exist including [pdbtools](https://github.com/harmslab/pdbtools), [pdb-tools](https://github.com/JoaoRodrigues/pdb-tools) and [PDBFixer](https://github.com/pandegroup/pdbfixer).
