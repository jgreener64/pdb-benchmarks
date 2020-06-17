// Benchmark the calculation of a distance in a PDB file
// The distance is the closest distance between any atoms of residues 50 and 60
//   of chain A in 1AKE
#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <cmath>
#include <string>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

static double distance(gemmi::Structure& st) {
    gemmi::Chain* a = st.first_model().find_chain("A");
    gemmi::Residue& r50 = a->find_residue_group(gemmi::SeqId(50,' '))[0];
    gemmi::Residue& r60 = a->find_residue_group(gemmi::SeqId(60,' '))[0];
    double min_dist_sq = INFINITY;
    for (const gemmi::Atom& a: r50.atoms)
        for (const gemmi::Atom& b: r60.atoms) {
            double d2 = a.pos.dist_sq(b.pos);
            if (d2 < min_dist_sq)
                min_dist_sq = d2;
        }
    return std::sqrt(min_dist_sq);
}

int main() {
    std::string pdb_filepath = "data/1AKE.pdb";
    gemmi::Structure st = gemmi::read_pdb_file(pdb_filepath);
    timespec tstart, tend;
    clock_gettime(CLOCK_REALTIME, &tstart);
    double d = distance(st);
    clock_gettime(CLOCK_REALTIME, &tend);
    assert(std::fabs(d - 9.57605) < 1e-5);
    printf("%.6f\n", (tend.tv_sec - tstart.tv_sec) + (tend.tv_nsec - tstart.tv_nsec) / 1e9);
    return 0;
}
