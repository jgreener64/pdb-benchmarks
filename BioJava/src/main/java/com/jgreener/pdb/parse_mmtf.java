// Benchmark the parsing of a MMTF file given as an argument

package com.jgreener.pdb;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.mmtf.MmtfActions;

import java.nio.file.Paths;

public class parse_mmtf
{
    public static void main( String[] args )
    {
        String mmtf_filepath = args[0];
        // Run once to trigger illegal reflective access warning
        try {
            Structure structure = MmtfActions.readFromFile(Paths.get(mmtf_filepath));
        } catch (Exception e) {
            e.printStackTrace();
        }
        long startTime = System.nanoTime();
        try {
            Structure structure = MmtfActions.readFromFile(Paths.get(mmtf_filepath));
        } catch (Exception e) {
            e.printStackTrace();
        }
        long endTime = System.nanoTime();
        System.out.println((endTime - startTime) / 1000000000.0);
    }
}
