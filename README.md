# PDB benchmarks

Open source software packages to parse files in various formats from the [Protein Data Bank](http://www.rcsb.org/pdb/home/home.do) (PDB) and manipulate protein structures exist in many languages, often as part of Bio* projects.

This repository aims to collate benchmarks for common tasks across various languages and packages. The collection of scripts may also be useful to get an idea how each package works.

Please feel free to contribute scripts from other packages, or submit improvements to the scripts already present - I'm looking for the fastest implementation for each software that makes use of the provided API.

Disclosure: I contributed the BioStructures.jl package to BioJulia and have made contributions to Biopython.

## Tests

* Parsing 2 PDB entries, taken from the benchmarking in [1], in the PDB, mmCIF and MMTF formats:
  * [1CRN](http://www.rcsb.org/pdb/explore/explore.do?structureId=1crn) - hydrophobic protein (327 atoms).
  * [1HTQ](http://www.rcsb.org/pdb/explore/explore.do?structureId=1htq) - multicopy glutamine synthetase (10 models of 97,872 atoms).
* Counting the number of alanine residues in adenylate kinase ([1AKE](http://www.rcsb.org/pdb/explore/explore.do?structureId=1ake)).
* Calculating the distance between residues 50 and 60 of chain A in adenylate kinase ([1AKE](http://www.rcsb.org/pdb/explore/explore.do?structureId=1ake)).
* Calculating the Ramachandran phi/psi angles in adenylate kinase ([1AKE](http://www.rcsb.org/pdb/explore/explore.do?structureId=1ake)).

[1] Gajda MJ, hPDB - Haskell library for processing atomic biomolecular structures in protein data bank format, *BMC Research Notes* 2013, **6**:483 - [link](http://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-6-483)

The PDB files can be downloaded to directory `data` by running `julia tools/download_data.jl` from this directory. If you have all the software installed, and compiled where applicable, you can run `sh tools/run_benchmarks.sh` from this directory to run the benchmarks and store the results in `benchmarks.csv`. The mean over a number of runs is taken for each benchmark to obtain the values below.

Benchmarks were carried out on an Intel Xeon CPU E5-1620 v3 3.50GHz x 8 processor with 32 GB 2400 MHz DDR4 RAM. The operating system was CentOS v8.1. Time is the elapsed time.

## Software

Currently, 15 packages across 7 programming languages are included in the benchmarks:
* [BioStructures](https://github.com/BioJulia/BioStructures.jl) v0.10.0 running on Julia v1.3.1; times measured after JIT compilation.
* [MIToS](https://github.com/diegozea/MIToS.jl) v2.4.0 running on Julia v1.3.1; times measured after JIT compilation.
* [Biopython](http://biopython.org/wiki/Biopython) v1.76 running on Python v3.7.4.
* [ProDy](http://prody.csb.pitt.edu) v1.10.11 running on Python v3.7.4.
* [MDAnalysis](http://www.mdanalysis.org) v0.20.1 running on Python v3.7.4.
* [biotite](https://www.biotite-python.org) v0.20.1 running on Python v3.7.4.
* [atomium](https://github.com/samirelanduk/atomium) v1.0.2 running on Python v3.7.4.
* [Bio3D](http://thegrantlab.org/bio3d/index.php) v2.4.1 running on R v3.6.2.
* [Rpdb](https://cran.r-project.org/web/packages/Rpdb/index.html) v2.3 running on R v3.6.2.
* [BioJava](https://biojava.org) v5.3.0 running on Java v1.8.0.
* [BioPerl](http://bioperl.org/index.html) v1.007002 running on Perl v5.26.3.
* [BioRuby](http://bioruby.org) v2.0.1 running on Ruby v2.5.5.
* [GEMMI](https://gemmi.readthedocs.io/en/latest/index.html) v0.3.6 compiled with gcc v8.3.1; there is also a Python interface but benchmarking was done in C++.
* [Victor](http://protein.bio.unipd.it/victor/index.php/Main_Page) v1.0 compiled with gcc v7.3.1.
* [ESBTL](http://esbtl.sourceforge.net/index.html) v1.0-beta01 compiled with gcc v7.3.1.

## Results

Note that direct comparison between these times should be treated with caution, as each package does something slightly different. For example, things that increase parsing time include:

* Parsing the header information.
* Accounting for disorder at both the atom and residue (point mutation) level.
* Forming a heirarchical model of the protein that makes access to specific residues, atoms etc. easier and faster after parsing.
* Allowing models in a file to have different atoms present.
* Checking that the file format is adhered to at various levels of strictness.

Each package supports these to varying degrees.

|                       | BioStructures | MIToS         | Biopython     | ProDy         | MDAnalysis    | biotite       | atomium       | Bio3D         | Rpdb          | BioJava       | BioPerl       | BioRuby       | GEMMI         | Victor        | ESBTL         |
| :-------------------- | :------------ | :------------ | :------------ | :------------ | :------------ | :------------ | :------------ | :------------ | :------------ | :------------ | :------------ | :------------ | :------------ | :------------ | :------------ |
| Parse PDB 1CRN / ms   | 0.76          | 0.61          | 7.8           | 3.1           | 5.0           | 4.6           | 6.9           | 9.8           | 9.7           | 7.5           | 42.0          | 21.0          | 0.23          | 7.6           | 2.8           |
| Parse PDB 1HTQ / s    | 3.5           | 3.0           | 17.0          | 2.2           | 1.5           | 4.8           | 20.0          | 2.9           | 14.0          | 1.2           | 49.0          | 13.0          | 0.34          | 11.0          | -             |
| Parse mmCIF 1CRN / ms | 1.9           | -             | 18.0          | -             | -             | 4.8           | 14.0          | -             | -             | 38.0          | -             | -             | 1.0           | -             | -             |
| Parse mmCIF 1HTQ / s  | 8.3           | -             | 46.0          | -             | -             | 8.9           | 36.0          | -             | -             | 17.0          | -             | -             | 1.7           | -             | -             |
| Parse MMTF 1CRN / ms  | 1.1           | -             | 5.2           | -             | -             | 1.3           | 4.7           | -             | -             | 4.5           | -             | -             | -             | -             | -             |
| Parse MMTF 1HTQ / s   | 3.4           | -             | 16.0          | -             | -             | 0.17          | 44.0          | -             | -             | 0.74          | -             | -             | -             | -             | -             |
| Count / ms            | 0.17          | 0.017         | 0.24          | 9.8           | 0.074         | -             | -             | 0.17          | 0.18          | -             | 0.47          | 0.076         | -             | -             | -             |
| Distance / ms         | 0.012         | 0.0042        | 0.28          | 51.0          | 0.65          | -             | -             | 19.0          | 1.3           | -             | 0.55          | 0.33          | -             | -             | -             |
| Ramachandran / ms     | 1.4           | -             | 120.0         | 220.0         | 1300.0        | -             | -             | -             | -             | -             | -             | -             | -             | -             | -             |
| Language              | Julia         | Julia         | Python        | Python        | Python        | Python        | Python        | R             | R             | Java          | Perl          | Ruby          | C++/Python    | C++           | C++           |
| License               | MIT           | MIT           | Biopython     | MIT           | GPLv2         | BSD 3-Clause  | MIT           | GPLv2         | GPLv2/GPLv3   | LGPLv2.1      | GPL/Artistic  | Ruby          | MPLv2/LGPLv3  | GPLv3         | GPLv3         |
| Hierarchichal parsing | ✓             | ✗             | ✓             | ✓             | ✓             | ✗             | ✓             | ✗             | ✗             | ✓             | ✓             | ✓             | ✓             | ✓             | ✓             |
| Supports disorder     | ✓             | ✗             | ✓             | ✗             | ✗             | ✗             | ✗             | ✗             | ✗             | ✗             | ✗             | ✗             | ✓             | ✗             | ✓             |
| Writes PDBs           | ✓             | ✓             | ✓             | ✓             | ✓             | ✓             | ✓             | ✓             | ✓             | ✓             | ✓             | ✗             | ✓             | ✓             | ✓             |
| Parses PDB header     | ✗             | ✗             | ✓             | ✓             | ✗             | ✗             | ✓             | ✓             | ✓             | ✓             | ✗             | ✓             | ✓             | ✓             | ✗             |
| Superimposition       | ✓             | ✓             | ✓             | ✓             | ✓             | ✓             | ✗             | ✓             | ✗             | ✓             | ✗             | ✗             | ✗             | ✗             | ✗             |
| PCA                   | ✗             | ✗             | ✗             | ✓             | ✓             | ✗             | ✗             | ✓             | ✗             | ✗             | ✗             | ✗             | ✗             | ✗             | ✗             |

Benchmarks as a plot, sorted by increasing time to parse PDB 1CRN:

![benchmarks](plot/plot.png "benchmarks")

## Parsing the whole PDB

It is instructive to run parsers over the whole PDB to see where errors arise. This approach has led to me submitting corrections for small mistakes (e.g. duplicate atoms, residue number errors) in a few PDB structures. As of July 2018, the PDB entries that error with the Biopython (permissive mode) and BioJulia parsers are:
* 4UDF - mmCIF file errors in Biopython and BioJulia due to duplicate C and O atoms in Lys91 of chains B, F etc.
* 1EJG - mmCIF file errors in Biopython due to blank and non-blank alt loc IDs at residue Pro22/Ser22.
* 5O61 - mmCIF file errors in Biopython due to an incorrect residue number at line 165,223.

Running Biopython in non-permissive mode picks up more potential problems such as broken chains and mixed blank/non-blank alt loc IDs. For further discussion on errors in PDB files see the Biopython [documentation](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf). The scripts to reproduce the whole PDB checking can be found in `checkwholepdb`. There is also a script to check recent PDB changes that can be run as a CRON job.

## Opinions

* For most purposes, particularly work on small numbers of files, the speed of the programs will not hold you back. In this case use the language/package you are most familiar with.
* For fast parsing, use a binary format such as [MMTF](http://mmtf.rcsb.org) or [binaryCIF](https://github.com/dsehnal/BinaryCIF).
* Whilst mmCIF became the standard PDB archive format in 2014, and is a very flexible archive format, that does not mean that it is the best choice for all of bioinformatics. mmCIF files take up a lot of space on disk, are slowest to read and do not yet work with many bioinformatics tools.
* If you are analysing ensembles of proteins then use packages with that functionality, such as ProDy or Bio3D, rather than writing the code yourself.

## Contributing

If you want to contribute benchmarks for a package, please make a pull request with the script(s) in a directory like the other packages. I will run the benchmarks again and change the README, thanks.

## Resources

* Information on file formats for [PDB](http://www.wwpdb.org/documentation/file-format), [mmCIF](http://mmcif.wwpdb.org) and [MMTF](https://github.com/rcsb/mmtf).
* Benchmarks for mmCIF parsing can be found [here](https://github.com/project-gemmi/mmcif-benchmark).
* A list of PDB parsing packages, particularly in C/C++, can be found [here](http://bioinf.org.uk/software/bioplib/libraries).
* The Biopython [documentation](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf) has a useful discussion on disorder at the atom and residue level.
* Sets of utility scripts exist including [pdbtools](https://github.com/harmslab/pdbtools), [pdb-tools](https://github.com/JoaoRodrigues/pdb-tools) and [PDBFixer](https://github.com/pandegroup/pdbfixer).
