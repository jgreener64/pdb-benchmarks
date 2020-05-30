// Benchmark the counting of alanine residues in a PDB file
#include <cstdio>
#include <string>

#include <time.h>

#include <chemfiles.hpp>

using namespace chemfiles;

static size_t count(Frame& frame) {
    auto selection = Selection("resname ALA");
    return selection.list(frame).size();
}

int main() {
    auto pdb_filepath = "data/1AKE.pdb";
    auto frame = Trajectory(pdb_filepath).read();

    timespec tstart, tend;
    clock_gettime(CLOCK_REALTIME, &tstart);
    count(frame);
    clock_gettime(CLOCK_REALTIME, &tend);

    printf("%.6f\n", (tend.tv_sec - tstart.tv_sec) + (tend.tv_nsec - tstart.tv_nsec) / 1e9);
    return 0;
}
