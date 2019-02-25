/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;

import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import utils.IOUtils;
import  lipengKang.analysis.KStringUtils;
import zhouyao.analysis.wheatHapMap.YaoIOUtils;

/**
 *
 * @author kanglipeng
 */
public class gerpAnalysis1 {

 
    HashMap<Integer, TIntArrayList> alignedBlockMap = new HashMap<>();
    HashMap<Integer, TIntArrayList> gerpBlockMap = new HashMap<>();
    TIntArrayList alignedPos = new TIntArrayList();
    ArrayList<Integer> alignedSegment = new ArrayList();
    ArrayList<String> gerpLine = new ArrayList();
    TIntArrayList gerpRange = new TIntArrayList();
    HashMap<Integer, ArrayList<String>> gerpLineMap = new HashMap<>();
    ArrayList<String> gerpPos = new ArrayList();
    int genomeSize = 0;
     String mafDir=null;
      String gerpDir=null;
      String outfile=null;
      String seqFile=null;
      String chromosomeNum=null;
       String subgenome=null;
      public gerpAnalysis1(String infiles){
        this.parseParameters(infiles);
      }

    public gerpAnalysis1(String mafDir, String gerpDir, String outfile, String seqFile,String chromosomeNum,String subgenome ) {
      
        this.countGenomeSize(seqFile,chromosomeNum,subgenome);
        this.getAlignedRange(mafDir,chromosomeNum,subgenome);
        this.parseGerpLine(gerpDir,subgenome);
        this.normalizeGerp(outfile,chromosomeNum,subgenome);

    }

    public void countGenomeSize(String seqFile,String chromosomeNum,String subgenome) {
        try {
            BufferedReader br = IOUtils.getTextReader(seqFile);
            String temp = null;

            while ((temp = br.readLine()) != null) {
                if (temp.startsWith(">")) {
                    continue;
                }
                genomeSize = genomeSize + temp.length();

            }

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
         System.out.println("...............normalizing wheat"+ subgenome+".chr"+chromosomeNum+" GERP++.rates....................");
        System.out.println("wheat"+ subgenome+".chr"+chromosomeNum+" genome size is " + String.valueOf(genomeSize) + "bp");

    }

    public HashMap<Integer, TIntArrayList> getAlignedRange(String mafDir,String chromosomeNum,String subgenome) {
     
       
        int alignedSegmentNum = 0;
        int gerpStartPos = 1;
        int gerpEndPos = 0;
        long alignedbps=0;

        try {
            BufferedReader br = IOUtils.getTextReader(mafDir);
            String temp = null;
            String[] tem = null;

            while ((temp = br.readLine()) != null) {

                List<String> tList = KStringUtils.fastSplit(temp);
                List<String> tListNew = new ArrayList<>();
                if (!temp.startsWith("s wheat"+subgenome+".chr")) {
                    continue;
                }
                if (!temp.startsWith("s wheat"+subgenome+".chr"+chromosomeNum)) {
                    for (int j = 0; j < tList.size(); j++) {

                        if (tList.get(j) != null && !tList.get(j).equals("")) {

                            tListNew.add(tList.get(j));
                        }
                    }
                    tem = tListNew.toArray(new String[tListNew.size()]);
                    if (tem == null || tem.length == 0) {
                        continue;
                    }
                     int alignedSeqSize = tem[6].length();
                    alignedbps=alignedbps+alignedSeqSize;
               
                    gerpStartPos = gerpEndPos + alignedSeqSize + 1;
                    gerpEndPos = gerpStartPos - 1;
                    continue;
                }

                for (int j = 0; j < tList.size(); j++) {

                    if (tList.get(j) != null && !tList.get(j).equals("")) {

                        tListNew.add(tList.get(j));
                    }
                }

                tem = tListNew.toArray(new String[tListNew.size()]);
                if (tem == null || tem.length == 0) {
                    continue;
                }
                int alignedSeqSize = tem[6].length();
                  alignedbps=alignedbps+alignedSeqSize;
                alignedPos = new TIntArrayList();
                gerpRange = new TIntArrayList();
                int startPos = Integer.parseInt(tem[2]) + 1;//maf produced by last aligner is 0-based, here we change to 1-based as CDS.fasta and gff3
                int endPos = startPos + alignedSeqSize - 1;
                alignedPos.add(startPos);
                alignedPos.add(endPos);
                gerpStartPos = gerpEndPos + 1;
                gerpEndPos = gerpStartPos + alignedSeqSize - 1;
                gerpRange.add(gerpStartPos);
                gerpRange.add(gerpEndPos);
                //StringBuilder blockNum = new StringBuilder();
                // blockNum.append(gerpStartPos).append("-").append(gerpEndPos);
                alignedSegment.add(alignedSegmentNum);
                alignedBlockMap.put(alignedSegmentNum, alignedPos);
                gerpBlockMap.put(alignedSegmentNum, gerpRange);
                //  gerpStartPos = gerpStartPos + alignedSeqSize;
                alignedSegmentNum++;
            }

           
        
                
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("alignedBlock established and .rates file should have at least "+alignedbps+" lines");
       
        return alignedBlockMap;//postion map from maf   1------list(start,end)  key is block number in gerpcol result;value is 1-based position in chromosome

    }

    public HashMap<Integer, ArrayList<String>> parseGerpLine(String gerpDir,String subgenome) {

        HashMap<Integer, TIntArrayList> gerpBlockMap = this.gerpBlockMap;
        ArrayList<Integer> alignedSegment = this.alignedSegment;
        System.out.println(gerpBlockMap.size());
        long lineFlag = 0;
        try {
            String temp = null;
            BufferedReader br = IOUtils.getTextReader(gerpDir);

            int gerpRangeIndex = 0;
            int GerpStartline = gerpBlockMap.get(alignedSegment.get(gerpRangeIndex)).get(0);
            int GerpEndline = gerpBlockMap.get(alignedSegment.get(gerpRangeIndex)).get(1);
            /*int GerpStartline = alignedBlockMap.get(alignedSegment.get(gerpRangeIndex)).get(0);
            int GerpEndline = alignedBlockMap.get(alignedSegment.get(gerpRangeIndex)).get(1);*/
            System.out.println("gerpstartline is "+GerpStartline);
             System.out.println(GerpEndline);
            while ((temp = br.readLine()) != null) {
                lineFlag++;
                if (lineFlag < GerpStartline) {
                    continue;
                }
                if (lineFlag < GerpEndline && lineFlag >= GerpStartline) {
                    gerpLine.add(temp);
                }
                if (lineFlag == GerpEndline) {
                    gerpLine.add(temp);
                    gerpLineMap.put(alignedSegment.get(gerpRangeIndex), gerpLine);
                    gerpLine = new ArrayList<>();
                    gerpRangeIndex++;
                    if (gerpRangeIndex < alignedSegment.size()) {
                        GerpStartline = gerpBlockMap.get(alignedSegment.get(gerpRangeIndex)).get(0);
                        GerpEndline = gerpBlockMap.get(alignedSegment.get(gerpRangeIndex)).get(1);
                    }

                }

       
            }
            System.out.println("gerpRange numbers is " + gerpRangeIndex);

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("wheat"+subgenome+" gerp score covers " + String.valueOf(lineFlag)+" bp");// the gerp.rates file covers
        return gerpLineMap;
    }

    public void normalizeGerp(String outfile,String chromosomeNum,String subgenome) {
     
        HashMap<Integer, ArrayList<String>> gerpLineMap = this.gerpLineMap;
        HashMap<Integer, TIntArrayList> alignedBlockMap = this.alignedBlockMap;
        ArrayList<Integer> alignedSegment = this.alignedSegment;
        int writedLineNum = 1;
        int genomeSize = this.genomeSize;
        System.out.println("alignedBlocksMap size is " + alignedBlockMap.size());
        System.out.println("gerpLineMap size is " + gerpLineMap.size());

        try {
            BufferedWriter bw = YaoIOUtils.getTextWriter(outfile);
            for (int i = 0; i < alignedSegment.size(); i++) {

                int startFlag = alignedBlockMap.get(alignedSegment.get(i)).get(0);
                int endFlag = alignedBlockMap.get(alignedSegment.get(i)).get(1);
                while (writedLineNum < startFlag) {

                    bw.write("0.00");
                    bw.write("\t");
                    bw.write("0.00");
                    bw.newLine();
                    writedLineNum++;

                }

                while (writedLineNum >= startFlag && writedLineNum <= endFlag) {
                    for (int j = 0; j < gerpLineMap.get(alignedSegment.get(i)).size(); j++) {
                        bw.write(gerpLineMap.get(alignedSegment.get(i)).get(j));
                        bw.newLine();
                        writedLineNum++;
                    }
                }
            }
            while (writedLineNum <= genomeSize) {
                writedLineNum++;
                bw.write("0.00");
                bw.write("\t");
                bw.write("0.00");
                bw.newLine();
            }

            System.out.println(writedLineNum - 1 + " normalized gerps in wheat"+subgenome+".chr"+chromosomeNum);
            bw.flush();
            bw.close();

        } catch (Exception e) {
            System.out.println("Error in Reading and Write MSA");
            e.printStackTrace();
        }
    }
    
    private void parseParameters(String infileS) {
        
        List<String> pLineList = new ArrayList<>();
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = null;
            boolean ifOut = false;
            if (!(temp = br.readLine()).equals("calling maf cover gene")) {
                ifOut = true;
            }
            if (!(temp = br.readLine()).equals("Author: clip")) {
                ifOut = true;
            }
            if (!(temp = br.readLine()).equals("Email: lipenkang@163.com")) {
                ifOut = true;
            }
            if (!(temp = br.readLine()).equals("Homepage: http://plantgeneticslab.weebly.com/")) {
                ifOut = true;
            }
            if (ifOut) {
                System.out.println("Thanks for using TEP.");
                System.out.println("Please keep the authorship in the parameter file. Program stops.");
                System.exit(0);
            }
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    continue;
                }
                if (temp.isEmpty()) {
                    continue;
                }
                pLineList.add(temp);
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        this.mafDir = pLineList.get(0);
        this.gerpDir = pLineList.get(1);
        this.outfile=pLineList.get(2);
        this.seqFile= pLineList.get(3);
        this.chromosomeNum=pLineList.get(4);
        this.subgenome=pLineList.get(5);
        
    }
    
    

    public static void main(String[] args) {
      /*  String mafDir="/data1/home/lipeng/result/CNEr/axt/A/formal/merge/gerp/wheatA1.maf";
        String gerpDir="/data1/home/lipeng/result/CNEr/axt/A/formal/merge/gerp/Triticeacea.cpy.mfa.rates";
        String outfile="/data1/home/lipeng/result/CNEr/axt/A/formal/merge/gerp/wheatA.chr1.gerp";
        String seqFile="/data1/home/lipeng/database/wheat/smRef/A/chr01A.fa";*/
       /* String mafDir = "/Users/kanglipeng/Desktop/test.maf";
        String gerpDir = "/Users/kanglipeng/Desktop/testout.txt";
        String outfile = "/Users/kanglipeng/Desktop/test.out";
        String seqFile = "/Users/kanglipeng/Desktop/test.fa";*/
       gerpAnalysis1 a= new gerpAnalysis1(args[0]);
       new gerpAnalysis1(a.mafDir, a.gerpDir, a.outfile, a.seqFile,a.chromosomeNum,a.subgenome );
       
    }

}
