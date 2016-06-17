// Benchmark the parsing of a PDB file given as an argument

#include <ESBTL/default.h>
#include <iostream>
#include <time.h>

int main( int argc, char* argv[] ) {
    std::string pdb_filepath = argv[1];
    std::cout.setstate(std::ios_base::failbit);
    clock_t tStart = clock();
    ESBTL::PDB_line_selector sel;
    std::vector<ESBTL::Default_system> systems;
    ESBTL::All_atom_system_builder<ESBTL::Default_system> builder(systems, sel.max_nb_systems());
    ESBTL::read_a_pdb_file(pdb_filepath, sel, builder, ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> >());
    clock_t tFinish = clock();
    std::cout.clear();
    printf("%.6f\n", (double)(tFinish - tStart)/CLOCKS_PER_SEC);
}
