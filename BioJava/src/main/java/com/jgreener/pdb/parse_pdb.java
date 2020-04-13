// Benchmark the parsing of a PDB file given as an argument

package com.jgreener.pdb;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.PDBFileReader;

public class parse_pdb
{
    public static void main( String[] args )
    {
        String pdb_filepath = args[0];
        // Run once to trigger illegal reflective access warning
        PDBFileReader pdbreader1 = new PDBFileReader();
        try {
            Structure structure = pdbreader1.getStructure(pdb_filepath);
        } catch (Exception e) {
            e.printStackTrace();
        }
        long startTime = System.nanoTime();
        PDBFileReader pdbreader2 = new PDBFileReader();
        try {
            Structure structure = pdbreader2.getStructure(pdb_filepath);
        } catch (Exception e) {
            e.printStackTrace();
        }
        long endTime = System.nanoTime();
        System.out.println((endTime - startTime) / 1000000000.0);
    }
}
