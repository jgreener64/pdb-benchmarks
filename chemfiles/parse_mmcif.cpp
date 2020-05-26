// Benchmark the parsing of a mmCIF file given as an argument
#include <cstdio>
#include <string>
#include <time.h>

#include <chemfiles.hpp>

int main(int argc, char* argv[]) {
    std::string mmcif_filepath = argv[1];

    struct timespec tstart, tend;
    clock_gettime(CLOCK_REALTIME, &tstart);
    auto frame = chemfiles::Trajectory(mmcif_filepath, 'r', "mmCIF").read();
    clock_gettime(CLOCK_REALTIME, &tend);

    printf("%.6f\n", (tend.tv_sec - tstart.tv_sec) + (tend.tv_nsec - tstart.tv_nsec) / 1e9);
    return 0;
}
