/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.MaizeGeneticLoad;

import pgl.infra.range.Range;
import gnu.trove.set.hash.TIntHashSet;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;


/**
 * Hold categorical range attributes, non-overlap ranges
 * @author Fei Lu
 */
public class Rangescp {
    public String annotation = "unknown";
    Range[] ranges = null;
    ArrayList<String> note = new ArrayList();
    
    public Rangescp () {}
   
    public Rangescp (Range[] ranges, String annotation) {
        this.ranges = ranges;
        this.annotation = annotation;
        this.sortByStartPosition();
    }
    
    public Rangescp (List<Range> rList, String annotation) {
        this(rList.toArray(new Range[rList.size()]), annotation);
    }
    
    public Rangescp (String infileS, IOFileFormat format) {
        this.readFile(infileS, format);
    }
    
    public static String getRangePositionString (Range[] rs) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < rs.length; i++) {
            sb.append(rs[i].start).append(":").append(rs[i].end).append(";");
        }
        if (rs.length>0) {
            sb.deleteCharAt(sb.length()-1);
        }
        return sb.toString();
    }
    
    public static String getRangePositionString (List<Range> rList) {
        return getRangePositionString(rList.toArray(new Range[rList.size()]));
    }

    public void clearNote () {
        this.note.clear();
    }
    
    public void addNote (String s) {
        this.note.add(s);
    }
    
    /**
     * Return new Ranges object from merging two objects
     * @param another
     * @return 
     */
    public Rangescp merge (Rangescp another) {  
        if (!this.annotation.equals(another.annotation)) {
            return null;
        }
        Range[] rs = new Range[this.getRangeNumber()+another.getRangeNumber()];
        System.arraycopy(this.ranges, 0, rs, 0, this.getRangeNumber());
        System.arraycopy(another.ranges, 0, rs, this.getRangeNumber(), another.getRangeNumber());
        Rangescp r = new Rangescp(rs, this.getAnnotation());
        for (int i = 0; i < this.note.size(); i++) r.addNote(this.note.get(i));
        for (int i = 0; i < another.note.size(); i++) r.addNote(another.note.get(i));
        r.sortByStartPosition();
        return r;
    }
    
    /**
     * Return new Ranges object of a chromosome
     * @param chromosome
     * @return 
     */
    public Rangescp getRangesByChromosome (int chromosome) {
        int startIndex = this.getStartIndexOfChromosome(chromosome);
        int endIndex = this.getEndIndexOfChromosome(chromosome);
        int size = endIndex-startIndex;
        Range[] rs = new Range[size];
        for (int i = 0; i < rs.length; i++) {
            rs[i] = this.ranges[i+startIndex];
        }
        return new Rangescp(rs, this.getAnnotation());
    }
    
    /**
     * Collapse the ranges to non-overlap ranges
     */
    public void collapse () {
        System.out.println("Starting collasping ranges");
        //this.sortByStartPosition();
        int[] chromosomes = this.getChromosomes();
        int[] mins = new int[chromosomes.length];
        int[] maxs = new int[chromosomes.length];
        for (int i = 0; i < mins.length; i++) {
            mins[i] = Integer.MAX_VALUE;
            maxs[i] = Integer.MIN_VALUE;
        }
        for (int i = 0; i < this.getRangeNumber(); i++) {
            int index = Arrays.binarySearch(chromosomes, this.getRangeChromosome(i));
            int v = this.getRangeStart(i);
            if (v < mins[index]) mins[index] = v;
            v = this.getRangeEnd(i);
            if (v > maxs[index]) maxs[index] = v;
        }
        ArrayList<Range> rList = new ArrayList();
        for (int i = 0; i < chromosomes.length; i++) {
            System.out.println("Start collasping chromosome " + String.valueOf(chromosomes[i]));
            int base = mins[i];
            int length = maxs[i] - mins[i];
            byte[] status = new byte[length];
            int startIndex = this.getStartIndexOfChromosome(chromosomes[i]);
            int endIndex = this.getEndIndexOfChromosome(chromosomes[i]);
            for (int j = startIndex; j < endIndex; j++) {
                for (int k = this.getRangeStart(j); k < this.getRangeEnd(j); k++) {
                    status[k-base] = 1;
                }
            }
            int current = 0;
            while (current < length) {
                if (status[current] != 0) {
                    int start = current+base;
                    while (current < length && status[current] == 1) {
                        current++;
                    }
                    int end = current+base;
                    rList.add(new Range(chromosomes[i], start, end));
                }
                current++;
            }
        }
        ranges = rList.toArray(new Range[rList.size()]);
        this.sortByStartPosition();
    }
    
    /**
     * Return the starting index of ranges on a certain chromosome
     * @param chromosome
     * @return 
     */
    public int getStartIndexOfChromosome (int chromosome) {
        Range query = new Range(chromosome, Integer.MIN_VALUE, Integer.MIN_VALUE);
        int hit  = Arrays.binarySearch(ranges, query);
        int index = -hit-1;
        if (this.getRangeChromosome(index) == chromosome) return index;
        return hit;
    }
    
    /**
     * Return the starting index of ranges on a certain chromosome, exclusive
     * @param chromosome
     * @return 
     */
    public int getEndIndexOfChromosome (int chromosome) {
        Range query = new Range(chromosome+1, Integer.MIN_VALUE, Integer.MIN_VALUE);
        int hit  = Arrays.binarySearch(ranges, query);
        int index = -hit-1;
        if (this.getRangeChromosome(index-1) == chromosome) return index;
        return hit;
    }
    
    public long getTotalRangeSize () {
        long sum = 0;
        for (int i = 0; i < this.getRangeNumber(); i++) {
            sum+=this.getRangeSize(i);
        }
        return sum;
    }
    
    public int getRangeSize (int index) {
        return ranges[index].getRangeSize();
    }
    
    public int getChromosomeNumber () {
        return this.getChromosomes().length;
    }
    
    public int[] getChromosomes () {
        TIntHashSet chrSet = new TIntHashSet();
        for (int i = 0; i < this.getRangeNumber(); i++) {
            chrSet.add(this.getRangeChromosome(i));
        }
        int[] chromosomes = chrSet.toArray();
        Arrays.sort(chromosomes);
        return chromosomes;
    }
    
    public void sortByStartPosition () {
        System.out.println("Start sorting ranges");
        Arrays.sort(ranges);
        System.out.println("Finished sorting ranges");
    }
    
    public String getAnnotation () {
        return this.annotation;
    }
    
    public int getRangeNumber () {
        return ranges.length;
    }
    
    public Range getRange (int index) {
        return ranges[index];
    }
    
    public int getRangeChromosome (int index) {
        return ranges[index].chr;
    }
    
    public int getRangeStart (int index) {
        return ranges[index].start;
    }
    
    public int getRangeEnd (int index) {
        return ranges[index].end;
    }
    
    /**
     * Return the index of a range which a position falls in, return negative value if the position is not in any range
     * @param chromosome
     * @param position
     * @return 
     */
    public int getRangeIndex (int chromosome, int position) {
        Range query = new Range(chromosome, position, position);
        int index = Arrays.binarySearch(ranges, query, new RangeSearchComparator());
        if (index < 0) {
            index = -index-2;
        }
        if (index > -1 && position < this.getRangeEnd(index) && this.getRangeChromosome(index) == chromosome) return index;
        return -2-index;
    }
    
    class RangeSearchComparator implements Comparator<Range>{
        @Override
        public int compare(Range o1, Range o2) {
            if (o1.chr == o2.chr) {
                if (o1.start == o2.start) return 0;
                else if (o1.start < o2.start) return -1;
                return 1;
            }
            return o1.chr-o2.chr;
        }
    }
    
    public boolean isInRanges (int chromosome, int position) {
        int index = this.getRangeIndex(chromosome, position);
        if (index < 0) return false;
        return true;
    }
    
    public void readFile (String infileS, IOFileFormat format) {
        if (format == IOFileFormat.Binary) {
            this.readBinaryFile(infileS);
        }
        else if (format == IOFileFormat.Text) {
            this.readTextFile(infileS);
        }
        else {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
        this.sortByStartPosition();
        System.out.println("Read in " + infileS);
    }
    
    private void readBinaryFile (String infileS) {
        try {
            DataInputStream dis = IOUtils.getBinaryReader(infileS);
            this.annotation = dis.readUTF();
            this.ranges = new Range[dis.readInt()];
            int noteNumber = dis.readInt();
            int n100 = ranges.length/100;
            for (int i = 0; i < noteNumber; i++) this.addNote(dis.readUTF());
            int nPro = 0;
            for (int i = 0; i < this.getRangeNumber(); i++) {
                ranges[i] = new Range(dis.readInt(), dis.readInt(), dis.readInt());
                if (i%n100 == 0) {
                    System.out.println("Read in %" + String.valueOf(nPro++));
                }
            }
            dis.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void readTextFile (String infileS) {
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            this.annotation = br.readLine();
            int n = Integer.valueOf(br.readLine());
            String temp;
            while ((temp = br.readLine()).startsWith("#")) {
                this.addNote(temp.replaceFirst("#", ""));
            }
            this.ranges = new Range[n];
            for (int i = 0; i < n; i++) {
                String[] tem = br.readLine().split("\t");
                ranges[i] = new Range(Integer.valueOf(tem[0]), Integer.valueOf(tem[1]), Integer.valueOf(tem[2]));
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
//    public void writeFile (String outfileS, IOFileFormat format) {
//        if (format == IOFileFormat.Binary) {
//            this.writeBinaryFile(outfileS);
//        }
//        else if (format == IOFileFormat.Text) {
//            this.writeTextFile(outfileS);
//        }
//        else {
//            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//        }
//    }
//    
    private void writeBinaryFile (String outfileS) {
        try {
            DataOutputStream dos = IOUtils.getBinaryWriter(outfileS);
            dos.writeUTF(this.annotation);
            dos.writeInt(this.getRangeNumber());
            dos.writeInt(this.note.size());
            for (int i = 0; i < this.note.size(); i++) dos.writeUTF(note.get(i));
            for (int i = 0; i < this.getRangeNumber(); i++) {
                dos.writeInt(this.getRangeChromosome(i));
                dos.writeInt(this.getRangeStart(i));
                dos.writeInt(this.getRangeEnd(i));
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
//    private void writeTextFile (String outfileS) {
//        try {
//            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
//            bw.write(String.valueOf(this.annotation));
//            bw.newLine();
//            bw.write(String.valueOf(this.getRangeNumber()));
//            bw.newLine();
//            for (int i = 0; i < this.note.size(); i++) {
//                bw.write("#"+this.note.get(i));
//                bw.newLine();
//            }
//            bw.write("Chromosome\tStart\tEnd");
//            bw.newLine();
//            for (int i = 0; i < this.getRangeNumber(); i++) {
//                bw.write(ranges[i].getOutputString());
//                bw.newLine();
//            }
//            bw.flush();
//            bw.close();
//        }
//        catch (Exception e) {
//            e.printStackTrace();
//        }
//    }
}
