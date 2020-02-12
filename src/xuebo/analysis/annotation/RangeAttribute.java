/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.annotation;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.util.Arrays;
import java.util.List;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;

/**
 * Hold range attributes with strand and value in each range, non-overlap ranges
 * @author Fei Lu
 */
public class RangeAttribute extends Ranges {
    /**MinValue means not defined. 1 means plus, 0 mean minus*/
    byte[] strands = null;
    /**NaN means not defined*/
    float[] values = null;
    
    public RangeAttribute (Ranges r, byte[] strands, float[] values) {
        super.ranges = r.ranges;
        super.annotation = r.annotation;
        this.strands = strands;
        this.values = values;
        if (values.length != ranges.length || strands.length != ranges.length) {
            System.out.println("The sizes of ranges and attribute do not match, program quit");
            System.exit(1);
        }
    }
    
    public RangeAttribute (Range[] ranges, String annotation, byte[] strands, float[] values) {
        super.ranges = ranges;
        super.annotation = annotation;
        this.strands = strands;
        this.values = values;
        if (values.length != ranges.length || strands.length != ranges.length) {
            System.out.println("The sizes of ranges and attribute do not match, program quit");
            System.exit(1);
        }
    }
    
    public RangeAttribute(String infileS, IOFileFormat format) {
        this.readFile(infileS, format);
    }
    
    public RangeAttribute merge (RangeAttribute another) {
        if (!this.annotation.equals(another.annotation)) {
            return null;
        }
        Range[] newr = new Range[this.getRangeNumber()+another.getRangeNumber()];
        byte[] news = new byte[this.getRangeNumber()+another.getRangeNumber()];
        float[] newv = new float[this.getRangeNumber()+another.getRangeNumber()];
        System.arraycopy(this.ranges, 0, newr, 0, this.getRangeNumber());
        System.arraycopy(another.ranges, 0, newr, this.getRangeNumber(), another.getRangeNumber());
        System.arraycopy(this.strands, 0, news, 0, this.getRangeNumber());
        System.arraycopy(another.strands, 0, news, this.getRangeNumber(), another.getRangeNumber());
        System.arraycopy(this.values, 0, newv, 0, this.getRangeNumber());
        System.arraycopy(another.values, 0, newv, this.getRangeNumber(), another.getRangeNumber());
        RangeAttribute ra = new RangeAttribute(newr, this.getAnnotation(), news, newv);
        for (int i = 0; i < this.note.size(); i++) ra.addNote(this.note.get(i));
        for (int i = 0; i < another.note.size(); i++) ra.addNote(another.note.get(i));
        ra.sortByStartPosition();
        return ra;
    }
    
    public RangeAttribute getRangeAttributeByChromosome (int chromosome) {
        int startIndex = this.getStartIndexOfChromosome(chromosome);
        int endIndex = this.getEndIndexOfChromosome(chromosome);
        int size = endIndex-startIndex;
        Range[] rs = new Range[size];
        byte[] ns = new byte[size];
        float[] nv = new float[size];
        for (int i = 0; i < rs.length; i++) {
            rs[i] = this.ranges[i+startIndex];
            ns[i] = this.strands[i+startIndex];
            nv[i] = this.values[i+startIndex];
        }
        return new RangeAttribute(new Ranges(rs, this.getAnnotation()), ns, nv);
    }
    
    public float getValue (int index) {
        return this.values[index];
    }
    
    public byte getStrand (int index) {
        return this.strands[index];
    }
    
    @Override
    public void collapse () {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    @Override
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
            for (int i = 0; i < noteNumber; i++) this.addNote(dis.readUTF());
            this.strands = new byte[this.getRangeNumber()];
            this.values = new float[this.getRangeNumber()];
            for (int i = 0; i < this.getRangeNumber(); i++) {
                ranges[i] = new Range(dis.readInt(), dis.readInt(), dis.readInt());
                strands[i] = dis.readByte();
                values[i] = dis.readFloat();
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
            this.strands = new byte[this.getRangeNumber()];
            this.values = new float[this.getRangeNumber()];
            for (int i = 0; i < n; i++) {
                String[] tem = br.readLine().split("\t");
                ranges[i] = new Range(Integer.valueOf(tem[0]), Integer.valueOf(tem[1]), Integer.valueOf(tem[2]));
                strands[i] = Byte.valueOf(tem[3]);
                values[i] = Float.valueOf(tem[4]);
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    @Override
    public void writeFile (String outfileS, IOFileFormat format) {
        if (format == IOFileFormat.Binary) {
            this.writeBinaryFile(outfileS);
        }
        else if (format == IOFileFormat.Text) {
            this.writeTextFile(outfileS);
        }
        else {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
    }
    
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
                dos.writeByte(this.getStrand(i));
                dos.writeFloat(this.getValue(i));
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void writeTextFile (String outfileS) {
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(String.valueOf(this.annotation));
            bw.newLine();
            bw.write(String.valueOf(this.getRangeNumber()));
            bw.newLine();
            for (int i = 0; i < this.note.size(); i++) {
                bw.write("#"+this.note.get(i));
                bw.newLine();
            }
            bw.write("Chromosome\tStart\tEnd\tStrand\tValue");
            bw.newLine();
            for (int i = 0; i < this.getRangeNumber(); i++) {
                bw.write(ranges[i].getOutputString()+"\t"+String.valueOf(strands[i])+"\t"+String.valueOf(values[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    @Override
    public void sortByStartPosition () {
        GenericSorting.quickSort(0, ranges.length, compByStartPosition, swapper); 
    }
    
    Swapper swapper = new Swapper() { 
        @Override
        public void swap(int a, int b) {
            Range r = ranges[a];
            ranges[a] = ranges[b];
            ranges[b] = r;
            byte tempB = strands[a];
            strands[a] = strands[b];
            strands[b] = tempB;
            float f = values[a];
            values[a] = values[b];
            values[b] = f;
        } 
    }; 

    IntComparator compByStartPosition = new IntComparator() { 
        @Override
        public int compare(int a, int b) { 
            return ranges[a].compareTo(ranges[b]);
        }    
    }; 
}
