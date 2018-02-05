/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xuebo.analysis.annotation;

//import static net.maizegenetics.dna.read.ReadUtils.baseCompleMap;

/**
 * Holding DNA sequence and providing operation functions for DNA
 * @author Fei Lu 
 */
public class Sequence {
    String seq = null;
    
    public Sequence (String seq) {
        this.seq = seq;
    }
    
    /**
     * Return a random DNA sequence
     * @param length
     * @return 
     */
    public static String getRandomSequence (int length) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < length; i++) {
            double r = Math.random();
            if (r < 0.25) sb.append("A");
            else if (r < 0.5) sb.append("T");
            else if (r < 0.75) sb.append("G");
            else sb.append("C");
        }
        return sb.toString();
    }
    
    /**
     * Return DNA sequence
     * @return 
     */
    public String getSeq () {
        return this.seq;
    }
    
    /**
     * Return an array of fragmented sequence
     * @param fragmentLength
     * @return 
     */
    public String[] getFragments (int fragmentLength) {
        int len = seq.length();
        int left = len%fragmentLength;
        int fragmentNumber;
        if (left == 0) {
            fragmentNumber = len/fragmentLength;
        }
        else {
            fragmentNumber = len/fragmentLength+1;
        }
        String[] fragments = new String[fragmentNumber];
        if (left == 0) {
            for (int i = 0; i < fragmentNumber; i++) {
                fragments[i] = seq.substring(i*fragmentLength, i*fragmentLength+fragmentLength);
            }
        }
        else {
            for (int i = 0; i < fragmentNumber-1; i++) {
                fragments[i] = seq.substring(i*fragmentLength, i*fragmentLength+fragmentLength);
            }
            fragments[fragmentNumber-1] = seq.substring((fragmentNumber-1)*fragmentLength, seq.length());
        }
        return fragments;
    }
    
    /**
     * Return fragmented sequences in Fasta format, easier to output
     * @param fragmentLength
     * @param prefix
     * @return 
     */
    public Fasta getFragmentsFasta (int fragmentLength, String prefix) {
        String[] fragments = this.getFragments(fragmentLength);
        String[] names = new String[fragments.length];
        int[] ids = new int[fragments.length];
        for (int i = 0; i < fragments.length; i++) {
            names[i] = prefix+String.valueOf(i+1);
            ids[i] = i+1;
        }
        return new Fasta(names, fragments, ids);
    }
    
    /**
    * Return reverse complementary sequence
    * @return 
    */
    public String getReverseComplementarySeq () {
        return this.getReverseComplementarySeq(0, this.seq.length());
    }
    
    /**
     * Return reverse complementary sequence
     * @param startIndex inclusive
     * @param endIndex exclusive
     * @return 
     */
    public String getReverseComplementarySeq (int startIndex, int endIndex) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < endIndex - startIndex; i++) {
            sb.append(ReadUtils.baseCompleMap.get(String.valueOf(this.seq.charAt(i+startIndex))));
        }
        return sb.reverse().toString();
    }
}
