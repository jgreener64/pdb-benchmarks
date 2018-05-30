# PDB benchmarks

Open source software packages to parse [Protein Data Bank](http://www.rcsb.org/pdb/home/home.do) (PDB) files and manipulate protein structures exist in many languages, often as part of Bio* projects.

This repository aims to collate benchmarks for common tasks across various languages and packages. The collection of scripts may also be useful to get an idea how each package works.

**NB - Some of the packages tested use out of date versions. I will look to fix this soon.**

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

* [BioJulia](https://biojulia.github.io/Bio.jl/) master branch running on Julia v0.5.0 (times measured after JIT compilation)
* [MIToS](https://github.com/diegozea/MIToS.jl) v1.2.3 running on Julia v0.4.6 (times measured after JIT compilation)
* [Biopython](http://biopython.org/wiki/Biopython) v1.66 running on Python v2.7.11
* [ProDy](http://prody.csb.pitt.edu/) v1.7 running on Python v2.7.11
* [MDAnalysis](http://www.mdanalysis.org/) v0.15.0 running on Python v2.7.11
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

|                       | BioJulia     | MIToS        | Biopython    | ProDy        | MDAnalysis   | Bio3D        | Rpdb         | BioPerl       | BioRuby      | Victor        | ESBTL        |
| :-------------------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :------------ | :----------- | :------------ | :----------- |
| Parse 1CRN / ms       | 1.4          | 2.4          | 9.1          | 2.2          | 6.4          | 31           | 19           | 63            | 25           | 10            | 6.8          |
| Parse 3JYV / s        | 0.49         | 0.74         | 1.0          | 0.28         | 0.80         | 14           | 2.2          | 3.8           | 0.98         | 7.7           | 0.95         |
| Parse 1HTQ / s        | 27           | 28           | 25           | 1.7          | 3.0          | 60           | 34           | 71            | 18           | 17            | -            |
| Count / ms            | 0.91         | 0.16         | 0.48         | 8.9          | 5.7          | 0.53         | 0.46         | 0.79          | 0.19         | -             | -            |
| Distance / ms         | 0.11         | 0.011        | 0.39         | 5.6          | 3.3          | 1.4          | 1.9          | 0.85          | 0.51         | -             | -            |
| Ramachandran / ms     | 13           | -            | 130          | 180          | 3500         | -            | -            | -             | -            | -             | -            |
| Language              | Julia        | Julia        | Python       | Python       | Python       | R            | R            | Perl          | Ruby         | C++           | C++          |
| Parses header         | ✗            | ✗            | ✓            | ✓            | ✗            | ✓            | ✓            | ✗             | ✓            | ✓             | ✗            |
| Hierarchichal parsing | ✓            | ✗            | ✓            | ✓            | ✓            | ✗            | ✗            | ✓             | ✓            | ✓             | ✓            |
| Supports disorder     | ✓            | ✗            | ✓            | ✗            | ✗            | ✗            | ✗            | ✗             | ✗            | ✗             | ✓            |
| Writes PDBs           | ✓            | ✓            | ✓            | ✓            | ✓            | ✓            | ✓            | ✓             | ✗            | ✓             | ✓            |
| Superimposition       | ✗            | ✓            | ✓            | ✓            | ✓            | ✓            | ✗            | ✗             | ✗            | ✗             | ✗            |
| PCA                   | ✗            | ✗            | ✗            | ✓            | ✓            | ✓            | ✗            | ✗             | ✗            | ✗             | ✗            |
| License               | MIT          | MIT          | Biopython    | MIT          | GPLv2        | GPLv2        | GPL          | GPL/Artistic  | Ruby         | GPLv3         | GPLv3        |

Benchmarks as a plot:

![benchmarks](plot/plot.png "benchmarks")


## Opinions

* For most purposes, particularly work on small numbers of files, the speed of the programs will not hold you back. In this case use the language/package you are most familiar with.
* If you are analysing ensembles of proteins use packages with that functionality, such as ProDy or Bio3D, rather than writing the code yourself.
* For fast parsing, consider using a binary format such as [MMTF](http://mmtf.rcsb.org/) or [binaryCIF](https://github.com/dsehnal/BinaryCIF).


## Contributing

If you want to contribute benchmarks for a package, please make a pull request with the script(s) in a directory like the other packages. I will run the benchmarks again and change the README, thanks.


## Plans

* Test BioJava, hPDB, possibly others.
* Finish Ramachandran scripts for remaining languages.
* Run methods on the whole of the PDB to look at how they deal with errors.
* Add benchmarks for parsing mmCIF, the standard PDB archive format.
* Add benchmarks for parsing binary formats, e.g. MMTF.


## Resources

* Benchmarks for mmCIF parsing can be found [here](https://github.com/project-gemmi/mmcif-benchmark).
* The PDB file format documentation can be found [here](http://www.wwpdb.org/documentation/file-format).
* A list of PDB parsing packages, particularly in C/C++, can be found [here](http://bioinf.org.uk/software/bioplib/libraries/).
* The Biopython [documentation](http://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ) has a useful discussion on disorder at the atom and residue level.
* Sets of utility scripts exist including [pdbtools](https://github.com/harmslab/pdbtools), [pdb-tools](https://github.com/JoaoRodrigues/pdb-tools) and [PDBFixer](https://github.com/pandegroup/pdbfixer).
