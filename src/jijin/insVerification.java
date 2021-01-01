package jijin;



import htsjdk.samtools.*;
import pgl.infra.utils.IOUtils;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
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
                if (i.getCigarString().contains("D")) {
                    isD.add(1);
                } else {
                    isD.add(0);
                }

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
        for (int x = 1; x < 1000; x = x + 1) {
            insVerification iv = new insVerification();
            System.out.println(x);
            iv.getBreakpoint("E:\\Desktop\\三代数据测试\\benchmark\\insertion\\pbmm2-svim_pass_50bp_INS.vcf.bed", x);
            String a=iv.getReads("E:\\Desktop\\三代数据测试\\benchmark\\insertion\\75.1.bam");
            System.out.println(a);
            bw.write(a);
            bw.newLine();
        }
        bw.flush();
        bw.close();
    }
}

