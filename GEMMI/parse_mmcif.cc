// Benchmark the parsing of a mmCIF file given as an argument

#include <gemmi/mmread.hpp>
#include <iostream>
#include <time.h>

int main( int argc, char* argv[] ) {
    std::string mmcif_filepath = argv[1];
    struct timespec tstart, tend;
    clock_gettime(CLOCK_REALTIME, &tstart);
    gemmi::Structure st = gemmi::read_structure_file(mmcif_filepath);
    clock_gettime(CLOCK_REALTIME, &tend);
    printf("%.6f\n", (tend.tv_sec - tstart.tv_sec) + (tend.tv_nsec - tstart.tv_nsec)/1E9);
}
