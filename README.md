# PDB benchmarks

**Warning: under development. Don't trust these benchmarks yet.**

Open source software to parse [Protein Data Bank](http://www.rcsb.org/pdb/home/home.do) (PDB) files and manipulate protein structures exist in many languages, often as part of Bio* projects.

This repository aims to collate benchmarks for common tasks across various languages. The collection of scripts may also be useful to get an idea how each package functions.

Please feel free to contribute scripts from other languages, or submit improvements the scripts already present - I'm looking for the fastest implementation for each software that makes use of the provided API.

Disclosure: I contributed the `Bio.Structure` module to BioJulia.


## Tests

* Parsing 3 PDB files, taken from the benchmarking in [1]:
  * [1CRN](http://www.rcsb.org/pdb/explore/explore.do?structureId=1crn) - hydrophobic protein (327 atoms)
  * [3JYV](http://www.rcsb.org/pdb/explore/explore.do?structureId=3jyv) - 80S rRNA (57,327 atoms)
  * [1HTQ](http://www.rcsb.org/pdb/explore/explore.do?structureId=1htq) - multicopy glutamine synthetase (10 models of 97,872 atoms)
* Counting the number of alanine residues in adenylate kinase ([1AKE](http://www.rcsb.org/pdb/explore/explore.do?structureId=1ake))
* Calculating the distance between residues 50 and 60 in adenylate kinase ([1AKE](http://www.rcsb.org/pdb/explore/explore.do?structureId=1ake))
* Calculating the Ramachandran phi/psi angles in adenylate kinase ([1AKE](http://www.rcsb.org/pdb/explore/explore.do?structureId=1ake))

[1] Gajda MJ, hPDB - Haskell library for processing atomic biomolecular structures in protein data bank format, *BMC Research Notes* 2013, **6**:483
Link ^^

The PDB files can be downloaded to directory `pdbs` by running `source download_pdbs.sh` from this directory.


## Software

* [BioJulia](https://biojulia.github.io/Bio.jl/) development version running on Julia v0.4.6 (times measured after JIT compilation)
* [Biopython](http://biopython.org/wiki/Biopython) v1.66 running on Python v2.7.10
* [ProDy](http://prody.csb.pitt.edu/) v1.7 running on Python v2.7.10
* [Bio3D](http://thegrantlab.org/bio3d/index.php) v2.2-2 running on R v3.2.2
* [Rpdb](https://cran.r-project.org/web/packages/Rpdb/index.html) v2.2 running on R v3.2.2
* [BioPerl](http://bioperl.org/index.html) v1.6.924 running on Perl v5.18.2
* [BioJava](http://biojava.org/) v... running on Java v1.8.0_91


## Comparison

Benchmarks were carried out on a 3.1 GHz Intel Core i7 processor with 16 GB 1867 MHz DDR3 RAM. The operating system was Mac OS X Yosemite 10.10.5.

Time measured is time to completion, so CPU...

Note you can't just compare as some parsing does different things, e.g. headers... for example ProDy doesn't do heirarchical below.

All times in seconds.

|                       | BioJulia     | Biopython    | ProDy        | Bio3D        | Rpdb         | BioPerl      | BioJava      |
| :-------------------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- |
| Parse 1CRN            |              |              |              |              |              |              |              |
| Parse 3JYV            |              |              |              |              |              |              |              |
| Parse 1HTQ            |              |              |              |              |              |              |              |
| Count                 |              |              |              |              |              |              |              |
| Distance              |              |              |              |              |              |              |              |
| Ramachandran          |              |              |              |              |              |              |              |
| ----                  |              |              |              |              |              |              |              |
| Language              | Julia        | Python       | Python       | R            | R            | Perl         | Java         |
| Parses header         | ✗            | ✓            | ✓            |              |              |              |              |
| Heirarchichal parsing | ✓            | ✓            |              |              |              |              |              |
| Writes PDBs           | ✓            | ✓            |              |              |              | ✓            |              |
| Superimposition       | ✗            | ✓            | ✓            |              |              |              |              |
| Supports disorder     | ✓            | ✓            |              |              |              |              |              |
| License               | MIT          | Biopython    | MIT          | GPLv2        | GPL          | GPL/Artistic | LGPL         |

![benchmarks](plot/plot.png "benchmarks")
