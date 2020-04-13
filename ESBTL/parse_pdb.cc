// Benchmark the parsing of a PDB file given as an argument

#include <ESBTL/default.h>
#include <iostream>
#include <time.h>

int main( int argc, char* argv[] ) {
    std::string pdb_filepath = argv[1];
    std::cout.setstate(std::ios_base::failbit);
    struct timespec tstart, tend;
    clock_gettime(CLOCK_REALTIME, &tstart);
    ESBTL::PDB_line_selector sel;
    std::vector<ESBTL::Default_system> systems;
    ESBTL::All_atom_system_builder<ESBTL::Default_system> builder(systems, sel.max_nb_systems());
    ESBTL::read_a_pdb_file(pdb_filepath, sel, builder, ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> >());
    clock_gettime(CLOCK_REALTIME, &tend);
    std::cout.clear();
    printf("%.6f\n", (tend.tv_sec - tstart.tv_sec) + (tend.tv_nsec - tstart.tv_nsec)/1E9);
}
