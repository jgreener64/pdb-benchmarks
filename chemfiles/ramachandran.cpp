// Benchmark the calculation of Ramachandran phi/psi angles from a PDB file
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>

#include <time.h>

#include <chemfiles.hpp>

using namespace chemfiles;

using ramachandran_t = std::pair<std::vector<double>, std::vector<double>>;
static ramachandran_t ramachandran(Frame& frame) {
    auto phi_selection = Selection("dihedrals: name(#1) C and name(#2) N and name(#3) CA and name(#4) C");
    auto phi_list = phi_selection.evaluate(frame);

    auto phi_angles = std::vector<double>();
    phi_angles.reserve(phi_list.size());

    for (const auto& phi: phi_list) {
        // 57.29578 to convert from radians to degrees
        phi_angles.push_back(frame.dihedral(phi[0], phi[1], phi[2], phi[3]) * 57.29578);
    }


    auto psi_selection = Selection("dihedrals: name(#1) N and name(#2) CA and name(#3) C and name(#4) N");
    auto psi_list = psi_selection.evaluate(frame);

    auto psi_angles = std::vector<double>();
    psi_angles.reserve(psi_list.size());

    for (const auto& phi: psi_list) {
        // 57.29578 to convert from radians to degrees
        phi_angles.push_back(frame.dihedral(phi[0], phi[1], phi[2], phi[3]) * 57.29578);
    }

    // FIXME: the sign of the angles is inverted w.r.t. the MDAnalysis results
    return {phi_angles, psi_angles};
}

int main() {
    auto pdb_filepath = "data/1AKE.pdb";
    auto frame = Trajectory(pdb_filepath).read();

    timespec tstart, tend;
    clock_gettime(CLOCK_REALTIME, &tstart);
    ramachandran(frame);
    clock_gettime(CLOCK_REALTIME, &tend);

    printf("%.6f\n", (tend.tv_sec - tstart.tv_sec) + (tend.tv_nsec - tstart.tv_nsec) / 1e9);
    return 0;
}
