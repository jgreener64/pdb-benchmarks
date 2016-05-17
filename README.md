# PDB benchmarks

**Warning: under development. Don't trust these benchmarks yet.**

Open source software to parse [Protein Data Bank](http://www.rcsb.org/pdb/home/home.do) (PDB) files and manipulate protein structures exist in many languages, often as part of Bio* projects.

This repository aims to collate benchmarks for common tasks across various languages.

Please feel free to contribute scripts from other languages, or improve the scripts already present - I'm looking for the fastest implementation for each software that makes use of the provided API.

Disclosure: I contributed the `Bio.Structure` module to BioJulia.


## Tests

* Parsing 3 PDB files, taken from the benchmarking in [1]:
  - 1CRN, hydrophobic protein (327 atoms)
  - 3JYV, 80S rRNA (57,327 atoms)
  - 1HTQ, multicopy glutamine synthetase (97,872 atoms)
* Counting the number of alanine residues in adenylate kinase (1AKE)
* Calculating the distance between residues 50 and 60 in adenylate kinase (1AKE)
* Calculating the Ramachandran phi/psi angles in adenylate kinase (1AKE)

[1] Gajda MJ, hPDB – Haskell library for processing atomic biomolecular structures in protein data bank format, *BMC Research Notes* 2013, **6**:483


## Software

* BioJulia development version running on Julia v0.4.6 (times measured after JIT compilation)
* Biopython v1.66 running on Python v2.7.10
* ProDy v1.7 running on Python v2.7.10
* Bio3D v2.2-2 running on R v3.2.2
* Rpdb v2.2 running on R v3.2.2


## Comparison

Benchmarks were carried out on a 3.1 GHz Intel Core i7 processor with 16 GB 1867 MHz DDR3 RAM. The operating system was Mac OS X Yosemite 10.10.5.

Time measured is time to completion, so CPU...

Note you can't just compare as some parsing does different things, e.g. headers... for example ProDy doesn't do heirarchical below.

All times in seconds.

|                       | BioJulia    | Biopython   | ProDy       | Bio3D       | Rpdb        |
| :-------------------- | :---------- | :---------- | :---------- | :---------- | :---------- |
| Parse 1CRN            |             |             |             |             |             |
| Parse 3JYV            |             |             |             |             |             |
| Parse 1HTQ            |             |             |             |             |             |
| Count                 |             |             |             |             |             |
| Distance              |             |             |             |             |             |
| Ramachandran          |             |             |             |             |             |
| ----                  |             |             |             |             |             |
| Language              | Julia       | Python      | Python      | R           | R           |
| Parses header         | ✗           | ✓           | ✓           |             |             |
| Heirarchichal parsing | ✓           | ✓           |             |             |             |
| Writes PDBs           | ✓           | ✓           |             |             |             |
| Superimposition       | ✗           | ✓           | ✓           |             |             |
| Supports disorder     | ✓           | ✓           |             |             |             |
| License               | MIT         | Biopython   | MIT         | GPLv2       | GPL         |

![benchmarks](plot/plot.png "benchmarks")
