package jijin;


import htsjdk.samtools.*;
import pgl.infra.utils.IOUtils;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

public class insVerification {
    public String chrom = null;
    public int breakpoint = 0;
    public String read=null;
    public String isIns = null ;
    public List<SAMRecord> reads = new ArrayList<SAMRecord>();
    public List<Integer>  isD = new ArrayList<>();
    public insVerification() {


    }

    public void getBreakpoint(String inputfile, int i) throws IOException {
        BufferedReader br = IOUtils.getTextReader(inputfile);
        String line = " ";
        int lines = 0;
        while ((line = br.readLine()) != null) {
            lines++;

            if (lines == i) {
                this.chrom = line.split("\t")[0];
                this.breakpoint = Integer.valueOf(line.split("\t")[1]);
                //System.out.println(this.breakpoint);
                break;
            }


        }
        //int[] a =new int[]{this.chrom,this.breakpoint};
        //return a;
    }

    public String getReads(String inputfile) {
        File bamFile = new File(inputfile);
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bamFile);
        //SAMRecordIterator reads = sr.iterator();
        //System.out.println(this.chrom);
        //System.out.println(this.breakpoint);
        SAMRecordIterator samIterator1 = sr.queryOverlapping(this.chrom, this.breakpoint, this.breakpoint + 1);
        while (samIterator1.hasNext()) {
            SAMRecord curSAM = samIterator1.next();
            //this.read=curSAM.getReadName();
            this.reads.add(curSAM);


            //System.out.println(curSAM.getReferenceName());

        }


        if (this.reads.isEmpty()) {  // 判断是否有reads
            //System.out.println(this.reads);
            isIns="0";
        } else {
            for (SAMRecord i : this.reads) {
                if (i.getCigar().isClipped()) {//i.getCigarString().contains("D") || i.getCigarString().contains("H")|| i.getCigarString().contains("S")
//                    if(i.getInferredInsertSize()<0){
//                        isD.add(0);
//                    }
//                    else{
//                        isD.add(1);
//                    }
                    //isD.add(1);
                    /*
                    System.out.println(i.getCigar().iterator());
                    System.out.println(this.breakpoint);
                    System.out.println(i.getStart()+i.getLengthOnReference());
                    final Iterator<CigarElement> cigar = i.getCigar().iterator();
                    */
                    List<Integer> breakpoint2 = new LinkedList<>();
                    List<CigarElement> cigar=i.getCigar().getCigarElements();
                    int bp=i.getAlignmentStart();
                    breakpoint2.add(bp);
                    for (CigarElement j : cigar) {
                        //System.out.println("len"+j.getLength());
                        //Integer value = Integer.valueOf(j.toString().substring(0, j.toString().length() - 1));
                        int value = j.getLength();
                        //System.out.println(j.getOperator());
                        String state = j.getOperator().toString();//j.toString().substring(j.toString().length() - 1, j.toString().length());
                        //System.out.println(state);
                        //breakpoint2.add(i.getAlignmentStart());

                        if (state == "M") {
                            bp = bp + value;
                            //System.out.println(bp);
                            breakpoint2.add(bp);


                        } else if (state == "I" || state == "P" || state == "N") {
                            continue;

                        } else if (state == "D") {
                            bp = bp + value;
                            breakpoint2.add(bp);

                        }else{
                            continue;
                        }

                    }
                    List<Integer> posacc = new ArrayList<>();
                    System.out.println(breakpoint2);
                    for (int m : breakpoint2){
                       posacc.add(m-this.breakpoint);
                    }
                    //System.out.println(posacc);
                    //未完
                    if(posacc.contains(3) || posacc.contains(2) || posacc.contains(1) || posacc.contains(0) || posacc.contains(-1) || posacc.contains(-2) || posacc.contains(-3) ){
                        isD.add(1);
                    }else{
                        isD.add(0);
                    }
                } else {
                    isD.add(0);
                }
//                    if(i.getInferredInsertSize()>0 && i.getInferredInsertSize()<350){
//                        isD.add(0);
//                    }
//                    else{

                   // }
            }

            int occurrences = Collections.frequency(isD,1);

            if(occurrences/isD.size()>=0.5){
                isIns="1";
            }
            else{
                isIns="0";
            }
        }
        //System.out.println(isD);
        return isIns;
    }




    public static void main (String[]args) throws IOException {

        FileWriter fw = new FileWriter("E:\\Desktop\\三代数据测试\\benchmark\\insertion\\ins.txt");
        BufferedWriter bw = new BufferedWriter(fw);
        int line = (int) Files.lines(Paths.get("E:\\Desktop\\三代数据测试\\benchmark\\insertion\\pbmm2-svim_pass_50bp_INS.vcf.bed")).count();
        System.out.println(line);
        for (int x = 1; x <= line ; x = x + 1) {
            insVerification iv = new insVerification();
            System.out.println(x);
            iv.getBreakpoint("E:\\Desktop\\三代数据测试\\benchmark\\insertion\\pbmm2-svim_pass_50bp_INS.vcf.bed", x);
            String a=iv.getReads("E:\\Desktop\\三代数据测试\\benchmark\\insertion\\75.bam");
            //System.out.println(a);
            bw.write(a);
            bw.newLine();
        }
        bw.flush();
        bw.close();
    }
}

