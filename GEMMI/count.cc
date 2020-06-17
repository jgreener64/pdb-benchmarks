// Benchmark the counting of alanine residues in a PDB file
#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <string>
#include <gemmi/pdb.hpp>

static int count(const gemmi::Structure& st) {
    int counter = 0;
    const std::string resname = "ALA";
    for (const gemmi::Chain& chain : st.first_model().chains)
        for (const gemmi::Residue& residue : chain.residues)
            if (residue.name == resname)
                ++counter;
    return counter;
}

int main() {
    std::string pdb_filepath = "data/1AKE.pdb";
    gemmi::Structure st = gemmi::read_pdb_file(pdb_filepath);
    timespec tstart, tend;
    clock_gettime(CLOCK_REALTIME, &tstart);
    int n = count(st);
    clock_gettime(CLOCK_REALTIME, &tend);
    assert(n == 38);
    printf("%.6f\n", (tend.tv_sec - tstart.tv_sec) + (tend.tv_nsec - tstart.tv_nsec) / 1e9);
    return 0;
}
