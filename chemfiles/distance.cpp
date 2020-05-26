// Benchmark the calculation of a distance in a PDB file
// The distance is the closest distance between any atoms of residues 50 and 60
//   of chain A in 1AKE
#include <cstdio>
#include <cmath>
#include <string>

#include <time.h>

#include <chemfiles.hpp>

using namespace chemfiles;

static double distance(Frame& frame) {
    // FIXME: this should use Selection("resid 50 and [chainname] A") which will
    // be available in chemfiles 0.10 (the next release)
    auto r50 = Selection("resid 50 and index < 1000").list(frame);
    auto r60 = Selection("resid 60 and index < 1000").list(frame);

    double min = INFINITY;
    for (auto i: r50) {
        for (auto j: r60) {
            auto r = frame.distance(i, j);
            if (r < min) {
                min = r;
            }
        }
    }

    return min;
}

int main() {
    auto pdb_filepath = "data/1AKE.pdb";
    auto frame = Trajectory(pdb_filepath).read();

    timespec tstart, tend;
    clock_gettime(CLOCK_REALTIME, &tstart);
    distance(frame);
    clock_gettime(CLOCK_REALTIME, &tend);

    printf("%.6f\n", (tend.tv_sec - tstart.tv_sec) + (tend.tv_nsec - tstart.tv_nsec) / 1e9);
    return 0;
}
