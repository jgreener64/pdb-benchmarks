// Benchmark the parsing of a PDB file given as an argument

#include <PdbLoader.h>
#include <Protein.h>
#include <iostream>
#include <time.h>

using namespace Victor::Biopool;
using namespace Victor;

int main( int argc, char* argv[] ) {
    string pdb_filepath = argv[1];
    clock_t tStart = clock();
    ifstream inFile( pdb_filepath.c_str() );
    // See options at
    // http://protein.bio.unipd.it/victor_doxygen/classVictor_1_1Biopool_1_1PdbLoader.html
    PdbLoader pl(inFile, true, true, false, true, false, false, false, true);
    Protein prot;
    prot.load( pl );
    clock_t tFinish = clock();
    printf("%.6f\n", (double)(tFinish - tStart)/CLOCKS_PER_SEC);
}
