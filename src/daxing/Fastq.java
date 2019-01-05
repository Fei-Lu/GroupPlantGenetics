package daxing;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import com.koloboke.collect.map.hash.HashByteByteMap;
import format.dna.BaseEncoder;
import utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/**
 *
 * @author xudaxing
 */
public class Fastq implements Swapper, IntComparator {
    private List<String> flowcellLaneIndex;
    private List<String> tileCoordinates;
    private String taxonName;
    private List<long[]> reads;
    private boolean isThisFastqOnlyContainOneFlowcellLaneIndex=false;

    //两种fastq构造器，一种是根据正常的fastq文件，另一种是根据简化板的fastq文件
    Fastq(String fastq1, String fastq2, String sampleID){
        File file1=new File(fastq1);
        File file2=new File(fastq2);
        if(file1.isFile() && file2.isFile()){
            if(fastq1.endsWith("sfq")){
                this.readSimpleFastq(fastq1, fastq2, sampleID);
            }
            else{
                this.readFastq(fastq1, fastq2, sampleID);
                this.isFromSameLibrary();
            }
        }
        else {
            System.out.println("please input two absolute paths of fastq files and a sampleID");
        }
    }

    private void readFastq(String fastq1, String fastq2, String sampleID){
        this.flowcellLaneIndex= new ArrayList<>();
        this.tileCoordinates=new ArrayList<>();
        this.reads=new ArrayList<>();
        this.taxonName=sampleID;
        BufferedReader br1;
        BufferedReader br2;
        String temp1;
        String temp2;
        int count=0;
        try{
            if(fastq1.endsWith("gz")){
                br1 = IOUtils.getTextGzipReader(fastq1);
                br2 = IOUtils.getTextGzipReader(fastq2);
            }
            else {
                br1 = IOUtils.getTextReader(fastq1);
                br2 = IOUtils.getTextReader(fastq2);
            }
            while((temp1=br1.readLine())!=null){
                temp2=br2.readLine();
                StringBuilder sb1=new StringBuilder();
                StringBuilder sb2=new StringBuilder();
                String[] line=temp1.split(":");
                sb1.append(line[2]).append("_").append(line[3]).append("_").append(line[9]);
                sb2.append(line[4]).append(":").append(line[5]).append(":").append(line[6]);
                this.flowcellLaneIndex.add(sb1.toString());
                this.tileCoordinates.add(sb2.toString().substring(0, sb2.length()-2));
                temp1=br1.readLine(); temp2=br2.readLine();
                HashByteByteMap ascIIByteMap = BaseEncoder.getAscIIByteMap();
                //reads原为146bp,现在改为128bp
                long[] tag = TagReads.getTagFromReads(temp1, temp2, ascIIByteMap, 5);
                this.reads.add(tag);
                count++;
                temp1=br1.readLine(); temp2=br2.readLine();
                temp1=br1.readLine(); temp2=br2.readLine();
            }
            System.out.println("A total of "+count+" tags were in "+ this.taxonName);
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }

    //这个构造函数创建以sfq结尾的fastq对象
    private void readSimpleFastq(String fastq1, String fastq2, String sampleID){
        this.isThisFastqOnlyContainOneFlowcellLaneIndex=true;
        this.flowcellLaneIndex= new ArrayList<>();
        this.tileCoordinates=new ArrayList<>();
        this.reads=new ArrayList<>();
        this.taxonName=sampleID;
        BufferedReader br1;
        BufferedReader br2;
        String temp1;
        String temp2;
        int count=0;
        try{
            br1 = IOUtils.getTextReader(fastq1);
            br2 = IOUtils.getTextReader(fastq2);
            while((temp1=br1.readLine())!=null){
                temp2=br2.readLine();
                String[] line=temp1.split("\t");
                flowcellLaneIndex.add(line[0]);
                tileCoordinates.add(line[1]);
                temp1=br1.readLine(); temp2=br2.readLine();
                HashByteByteMap ascIIByteMap = BaseEncoder.getAscIIByteMap();
                //reads原为146bp,现在改为150bp
                long[] tag = TagReads.getTagFromReads(temp1, temp2, ascIIByteMap, 5);
                this.reads.add(tag);
                count++;
            }
            br1.close();
            br2.close();
            System.out.println("A total of "+count+" tags were in "+this.getFlowcellLaneIndex()+"_"+this.taxonName+"_"+"R1.sfq");
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }

    private void isFromSameLibrary(){
        Set<String> s=new TreeSet<>(flowcellLaneIndex);
        if(s.size()==1) {
            isThisFastqOnlyContainOneFlowcellLaneIndex=true;
            System.out.println(this.getFlowcellLaneIndex()+"_"+this.taxonName+" fastq file was from the same library: "+flowcellLaneIndex.get(0));
        }
        else {
            isThisFastqOnlyContainOneFlowcellLaneIndex=false;
            System.out.println(this.taxonName+" fastq file was not from the same library"+", it has "+s.size()+" flowcellLaneIndex");
        }
    }

    public void addBarcode(String wellID, String barcodeForR1, String barcodeForR2){
        String read1;
        String read2;
        List<long[]> tags= new ArrayList<>();
        for(int i=0;i<reads.size();i++){
            String[] read1read2=TagReads.getReadsFromTag(reads.get(i),(short)160,(short)160);
            read1=read1read2[0];
            read2=read1read2[1];
            StringBuilder sb1=new StringBuilder();
            StringBuilder sb2=new StringBuilder();
            sb1.append(barcodeForR1).append(read1);
            sb2.append(barcodeForR2).append(read2);
            HashByteByteMap ascIIByteMap = BaseEncoder.getAscIIByteMap();
            tags.add(TagReads.getTagFromReads(sb1.toString(), sb2.toString(), ascIIByteMap, 5));
        }
        this.reads=tags;
    }

    @Override
    public void swap (int index1, int index2) {
        String s=flowcellLaneIndex.get(index1);
        flowcellLaneIndex.set(index1, flowcellLaneIndex.get(index2));
        flowcellLaneIndex.set(index2, s);
        s=tileCoordinates.get(index1);
        tileCoordinates.set(index1, tileCoordinates.get(index2));
        tileCoordinates.set(index2, s);
        long[] l=reads.get(index1);
        reads.set(index1, reads.get(index2));
        reads.set(index2, l);
    }

    @Override
    public int compare (int index1, int index2) {
        int value=flowcellLaneIndex.get(index1).compareTo(flowcellLaneIndex.get(index2));
        return value;
    }

    public void sort(){
        GenericSorting.quickSort(0, this.getNumberOfReads(), this, this);
    }

    public int getNumberOfReads(){
        return reads.size();
    }

    //将fastq文件分别按照flowcellLaneIndex进行输出
    public void writeFlowcellLaneIndexFastq(String outFileS){
        this.sort();
        Set<String> s=this.getAllDifferentFlowcellLaneIndex();
        if(s.size()<=1){
            System.out.println("This writeFlowcellLaneIndexFastq method must be used for fastq that have at least two flowcellLaneIndex");
            System.exit(1);
        }
        List<String> l = new ArrayList<>(s);
        int index1=0;
        int index2=0;
        try{
            //在相邻的flowcellLaneIndex处，遍历flowcellLaneIndex列表
            for(int i=0; i<(l.size()-1); i++){
                BufferedWriter bw1=IOUtils.getTextWriter(outFileS+"/"+l.get(i)+"_"+taxonName+"_"+"R1.sfq");
                BufferedWriter bw2=IOUtils.getTextWriter(outFileS+"/"+l.get(i)+"_"+taxonName+"_"+"R2.sfq");
                index1=flowcellLaneIndex.indexOf(l.get(i));
                index2=flowcellLaneIndex.indexOf(l.get(i+1));
                for(int j=index1; j<index2; j++){
                    String[] read1read2=TagReads.getReadsFromTag(reads.get(j),(short)160,(short)160);
                    bw1.write(flowcellLaneIndex.get(j)+"\t");
                    bw1.write(tileCoordinates.get(j));bw1.newLine();
                    bw1.write(read1read2[0]);bw1.newLine();
                    bw2.write(flowcellLaneIndex.get(j)+"\t");
                    bw2.write(tileCoordinates.get(j));bw2.newLine();
                    bw2.write(read1read2[1]);bw2.newLine();
                }
                bw2.flush();bw2.close();
                bw1.flush();bw1.close();
            }
            //将最后一个flowcellLaneIndex对应的reads写入文件
            BufferedWriter bw1=IOUtils.getTextWriter(outFileS+"/"+l.get(l.size()-1)+"_"+taxonName+"_"+"R1.sfq");
            BufferedWriter bw2=IOUtils.getTextWriter(outFileS+"/"+l.get(l.size()-1)+"_"+taxonName+"_"+"R2.sfq");
            for(int k=index2;k<flowcellLaneIndex.size();k++){
                String[] read1read2=TagReads.getReadsFromTag(reads.get(k),(short)160,(short)160);
                bw1.write(flowcellLaneIndex.get(k)+"\t");
                bw1.write(tileCoordinates.get(k));bw1.newLine();
                bw1.write(read1read2[0]);bw1.newLine();
                bw2.write(flowcellLaneIndex.get(k)+"\t");
                bw2.write(tileCoordinates.get(k));bw2.newLine();
                bw2.write(read1read2[1]);bw2.newLine();
            }
            bw2.flush();bw2.close();
            bw1.flush();bw1.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }

    public Set<String> getAllDifferentFlowcellLaneIndex(){
        Set<String> t=new TreeSet<>(flowcellLaneIndex);
        //System.out.println(taxonName+" have "+t.size()+" flowcellLaneIndex");
        return t;
    }

    public String getFlowcellLaneIndex(){
        if(this.isThisFastqOnlyContainOneFlowcellLaneIndex){
            return this.flowcellLaneIndex.get(0);
        }
        else{
            System.out.println("Current fastq does not have only one flowCellLaneIndex");
            System.exit(1);
        }
        return null;
    }

    public void mergeFastq(Fastq fq){
        if((this.isThisFastqOnlyContainOneFlowcellLaneIndex) &&(fq.isThisFastqOnlyContainOneFlowcellLaneIndex)){
            if(this.getFlowcellLaneIndex().equals(fq.getFlowcellLaneIndex())){
                this.flowcellLaneIndex.addAll(fq.flowcellLaneIndex);
                this.reads.addAll(fq.reads);
                this.tileCoordinates.addAll(fq.tileCoordinates);
                this.taxonName=this.getFlowcellLaneIndex();
                this.isThisFastqOnlyContainOneFlowcellLaneIndex=true;
            }
            else{
                System.out.println("These two fastq files did not hava same flowcellLaneIndex");
                System.exit(1);
            }
        }
        else{
            System.out.println("These two fastq files have at least two flowcellLaneIndex");
            System.exit(1);
        }
    }

    public void writeFastq(String r1OutFileS, String r2OutFileS){
        String[] read1Read2=null;
        try(BufferedWriter bw1=IOUtils.getTextWriter(r1OutFileS);
            BufferedWriter bw2=IOUtils.getTextWriter(r2OutFileS)
        ){
            for(int i=0;i<this.reads.size();i++){
                read1Read2=TagReads.getReadsFromTag(this.reads.get(i), (short)160, (short)160);
                bw1.write(this.getFlowcellLaneIndex()+"\t"+this.tileCoordinates.get(i));bw1.newLine();
                bw2.write(this.getFlowcellLaneIndex()+"\t"+this.tileCoordinates.get(i));bw2.newLine();
                bw1.write(read1Read2[0]); bw1.newLine();
                bw2.write(read1Read2[1]); bw2.newLine();
                bw1.write("+");bw1.newLine();
                bw2.write("+");bw2.newLine();
                bw1.write("B"); bw1.newLine();
                bw2.write("B"); bw2.newLine();
            }
            bw1.flush();bw2.flush();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }

    public String getTaxonName(){
        return this.taxonName;
    }

}

