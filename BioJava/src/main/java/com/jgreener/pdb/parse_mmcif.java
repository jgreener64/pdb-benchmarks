// Benchmark the parsing of a mmCIF file given as an argument

package com.jgreener.pdb;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.StructureIOFile;
import org.biojava.nbio.structure.io.MMCIFFileReader;

public class parse_mmcif
{
    public static void main( String[] args )
    {
        String mmcif_filepath = args[0];
        // Run once to trigger illegal reflective access warning
        StructureIOFile reader1 = new MMCIFFileReader();
        try {
            Structure structure = reader1.getStructure(mmcif_filepath);
        } catch (Exception e) {
            e.printStackTrace();
        }
        long startTime = System.nanoTime();
        StructureIOFile reader2 = new MMCIFFileReader();
        try {
            Structure structure = reader2.getStructure(mmcif_filepath);
        } catch (Exception e) {
            e.printStackTrace();
        }
        long endTime = System.nanoTime();
        System.out.println((endTime - startTime) / 1000000000.0);
    }
}
