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
import pgl.infra.utils.IOUtils;
import lipengKang.analysis.KStringUtils;

import zhouyao.analysis.wheatHapMap.YaoIOUtils;

/**
 *
 * @author kanglipeng
 */
public class gerpAnalysis {

    HashMap<Integer, TIntArrayList> alignedBlockMap = new HashMap<>();
    HashMap<Integer, ArrayList<Long>> gerpBlockMap = new HashMap<>();
    TIntArrayList alignedPos = new TIntArrayList();
    ArrayList<Integer> alignedSegment = new ArrayList();
    ArrayList<String> gerpLine = new ArrayList();
    ArrayList<Long> gerpRange = new ArrayList();
    HashMap<Integer, ArrayList<String>> gerpLineMap = new HashMap<>();
    ArrayList<String> gerpPos = new ArrayList();
    long genomeSize = 0;
    String mafDir = null;
    String gerpDir = null;
    String outfile = null;
    String seqFile = null;
    String chromosomeNum = null;
    String subgenome = null;

    public gerpAnalysis(String infiles) {
        this.parseParameters(infiles);
    }

    public gerpAnalysis(String mafDir, String gerpDir, String outfile, String seqFile, String chromosomeNum, String subgenome) {

        this.countGenomeSize(seqFile, chromosomeNum, subgenome);
       //this.getAlignedRange(mafDir, chromosomeNum, subgenome);
        //this.parseGerpLine(gerpDir, subgenome);
        //this.normalizeGerp(outfile, chromosomeNum, subgenome);

    }

    public void countGenomeSize(String seqFile, String chromosomeNum, String subgenome) {
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
        System.out.println("...............normalizing wheat" + subgenome + ".chr" + chromosomeNum + " GERP++.rates....................");
        System.out.println("wheat" + subgenome + ".chr" + chromosomeNum + " genome size is " + String.valueOf(genomeSize) + "bp");

    }

    public HashMap<Integer, TIntArrayList> getAlignedRange(String mafDir, String chromosomeNum, String subgenome) {
        /*
a score=3930.0
s wheatD.chr4        83 10 + 509857067 GAG-TGCATA
s wheatA.chr4 141979723 10 - 744588157 GAGGTGCATA
s wheatB.chr4   2952499 10 + 673617499 GAGGTGCATA
s barley.chr4 645142670 10 - 647060158 GAGCCACACA
 gerp.rates file is follow the sequece in  multiple maf file. So the first 
  alignedRange is 83-92 and the first gerpBlock is 1-10
        
         */

        int alignedSegmentNum = 0;
        long gerpStartPos = 0;
        long gerpEndPos = 0;

        int endPos = 0;


        try {
            BufferedReader br = IOUtils.getTextReader(mafDir);
            String temp = null;
            String[] tem = null;
            String actualSeq=null;
            while ((temp = br.readLine()) != null) {

                List<String> tList = KStringUtils.fastSplit(temp);
                List<String> tListNew = new ArrayList<>();
                //find reference lines in mafs
                if (!temp.startsWith("s wheat" + subgenome + ".chr")) {
                    continue;
                }
                //find specifc chromosomes of reference  in mafs(not normalizing chromosome)
                if (!temp.startsWith("s wheat" + subgenome + ".chr" + chromosomeNum)) {
                    // start with "s wheatA.chr" but not start with "s wheatA.chr1"; add all non-empty string length to gerpBlock
                    for (int j = 0; j < tList.size(); j++) {

                        if (tList.get(j) != null && !tList.get(j).equals("")) {

                            tListNew.add(tList.get(j));
                        }
                    }
                    tem = tListNew.toArray(new String[tListNew.size()]);
                    if (tem == null || tem.length == 0) {
                        continue;
                    }

                    //filter "N" and "n" in mafs
                    //filter "-" in mafs
                    int gerpClaculatedSize = tem[6].replaceAll("N", "").replaceAll("n", "").replaceAll("-", "").length();
                    gerpStartPos = gerpEndPos + gerpClaculatedSize + 1;
                    gerpEndPos = gerpStartPos - 1;
                    continue;
                }
                //find specifc chromosomes of reference  in mafs(normalizing chromosome)
                for (int j = 0; j < tList.size(); j++) {

                    if (tList.get(j) != null && !tList.get(j).equals("")) {

                        tListNew.add(tList.get(j));
                    }
                }

                tem = tListNew.toArray(new String[tListNew.size()]);
                if (tem == null || tem.length == 0) {
                    continue;
                }
                int startPos = Integer.parseInt(tem[2]) + 1;
               actualSeq= tem[6].replaceAll("-", "");
                char[] seqArr = actualSeq.toCharArray();
                if(seqArr.length==0) continue;
                if (seqArr.length == 1) {
                    startPos = startPos;
                    //single base sequence parsing
                    if ( seqArr[0] != 'N' && seqArr[0] != 'n') {
                        endPos = startPos;
                        alignedPos.add(startPos);
                        alignedPos.add(endPos);
                        gerpStartPos = gerpEndPos + 1;
                        alignedSegment.add(alignedSegmentNum);
                        alignedBlockMap.put(alignedSegmentNum, alignedPos);
                        gerpEndPos = gerpStartPos + endPos - alignedPos.get(0);
                        gerpRange.add(gerpStartPos);
                        gerpRange.add(gerpEndPos);
                        gerpBlockMap.put(alignedSegmentNum, gerpRange);
                        alignedSegmentNum++;
                        gerpRange = new ArrayList();
                        alignedPos = new TIntArrayList();
                    }
                } else {

                    //>=2 bases sequence parsing; using 2-mer consist of rediculous bases(n,N,-)and normal bses to calculate chromosome position and gerp position of an sub-aligned block
                    //"-A" type 2-mer is used to start a block and "A-" type 2-mer is used to end a block and "AA" type 2-mer is used to extending
                    //sequences ends with -A or AA and sequences start with A- or AA are error to finish start-extend-end flow, so these special situation are judged additionally
                    
                    for (int l = 0; l < seqArr.length - 1; l++) {
                        //sequences consist of 2-mers full of N and n ,such as "nn"
                        if ((seqArr[l] == 'n' || seqArr[l] == 'N') && (seqArr[l + 1] == 'n' || seqArr[l + 1] == 'N')) {
                            startPos = startPos + 1;
                        }
                        //sequences consist of 2-mers full of N ,n combined with A,T,C and G; such as "nA","NC","nT"
                        if ((seqArr[l] == 'n' || seqArr[l] == 'N' )&& ( seqArr[l + 1] != 'n' && seqArr[l + 1] != 'N')) {//nnnA
                            // sequences end with above-mentioned "nA" type 2-mer
                            if (l == seqArr.length - 2) {
                                startPos = startPos + 1;
                                endPos = startPos;
                                alignedPos.add(startPos);
                                alignedPos.add(endPos);
                                gerpStartPos = gerpEndPos + 1;
                                gerpEndPos = gerpStartPos + endPos - alignedPos.get(0);
                                alignedSegment.add(alignedSegmentNum);
                                alignedBlockMap.put(alignedSegmentNum, alignedPos);
                                gerpRange.add(gerpStartPos);
                                gerpRange.add(gerpEndPos);
                                gerpBlockMap.put(alignedSegmentNum, gerpRange);
                                gerpRange = new ArrayList();
                                alignedPos = new TIntArrayList();
                                alignedSegmentNum++;
                            } else {
                                //above-mentioned 2-mer exist in middle of sequences
                                startPos = startPos + 1;
                                gerpStartPos = gerpEndPos + 1;
                                alignedPos.add(startPos);
                                gerpRange.add(gerpStartPos);
                            }
                        }
                        //sequences consist of 2-mers full of A,T,C and G combined with N, n and -
                        if (( seqArr[l] != 'n' && seqArr[l] != 'N') && ( seqArr[l + 1] == 'n' || seqArr[l + 1] == 'N')) {//AAAnnn //Annn
                            // sequences start with above-mentioned "An" type 2-mer
                            if (l == 0) {
                                endPos = startPos;
                                gerpStartPos = gerpEndPos + 1;
                                alignedPos.add(startPos);
                                alignedPos.add(endPos);
                                gerpEndPos = gerpStartPos + endPos - alignedPos.get(0);
                                gerpRange.add(gerpStartPos);
                                gerpRange.add(gerpEndPos);
                                alignedSegment.add(alignedSegmentNum);
                                gerpBlockMap.put(alignedSegmentNum, gerpRange);
                                alignedBlockMap.put(alignedSegmentNum, alignedPos);
                                gerpRange = new ArrayList();
                                alignedPos = new TIntArrayList();
                                alignedSegmentNum++;
                                startPos = startPos + 1;
                            } else {
                                endPos = startPos;
                                alignedPos.add(endPos);
                                gerpStartPos = gerpEndPos + 1;
                                gerpEndPos = gerpStartPos + endPos - alignedPos.get(0);
                                startPos = startPos + 1;
                                gerpRange.add(gerpEndPos);
                                alignedSegment.add(alignedSegmentNum);
                                gerpBlockMap.put(alignedSegmentNum, gerpRange);
                                alignedBlockMap.put(alignedSegmentNum, alignedPos);
                                gerpRange = new ArrayList();
                                alignedPos = new TIntArrayList();
                                alignedSegmentNum++;
                            }
                        }
                        //sequences consist of 2-mers full of normal bases---A,T,C,G
                        if ((seqArr[l] != 'n' && seqArr[l] != 'N' )&& (seqArr[l + 1] != 'n' && seqArr[l + 1] != 'N')) {//---AAA //AA---
                            startPos = startPos + 1;
                            // sequence ends with "AA" type but suceed to start and fail to end
                            if (l == seqArr.length - 2 && l != 0) {
                                endPos = startPos;
                                alignedPos.add(endPos);
                                gerpEndPos = gerpStartPos + endPos - alignedPos.get(0);
                                gerpRange.add(gerpEndPos);
                                alignedSegment.add(alignedSegmentNum);
                                alignedBlockMap.put(alignedSegmentNum, alignedPos);
                                gerpBlockMap.put(alignedSegmentNum, gerpRange);
                                gerpRange = new ArrayList();
                                alignedPos = new TIntArrayList();
                                alignedSegmentNum++;
                            }
                            // sequence start with "AA" type but fail to start and suceed to end
                            if (l == 0 && l != seqArr.length - 2) {
                                gerpStartPos = gerpEndPos + 1;
                                alignedPos.add(startPos - 1);
                                gerpRange.add(gerpStartPos);
                            }
                            // sequence ends with "AA" type but fail to start and fail to end
                            if (l == 0 && l == seqArr.length - 2) {
                                endPos = startPos;
                                alignedPos.add(startPos - 1);
                                alignedPos.add(endPos);
                                gerpStartPos = gerpEndPos + 1;
                                gerpEndPos = gerpStartPos + endPos - alignedPos.get(0);
                                gerpRange.add(gerpStartPos);
                                gerpRange.add(gerpEndPos);
                                alignedSegment.add(alignedSegmentNum);
                                alignedBlockMap.put(alignedSegmentNum, alignedPos);
                                gerpBlockMap.put(alignedSegmentNum, gerpRange);
                                gerpRange = new ArrayList();
                                alignedPos = new TIntArrayList();
                                alignedSegmentNum++;
                            }
                        }

                    }
                }

            }

            br.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("alignedBlock established and .rates file should have at least " + gerpBlockMap.get(alignedSegment.size() - 1).get(1) + " lines");

        return alignedBlockMap;//postion map from maf   1------list(start,end)  key is block number in gerpcol result;value is 1-based position in chromosome

    }

    public HashMap<Integer, ArrayList<String>> parseGerpLine(String gerpDir, String subgenome) {

        HashMap<Integer, ArrayList<Long>> gerpBlockMap = this.gerpBlockMap;
        ArrayList<Integer> alignedSegment = this.alignedSegment;
        System.out.println("gerpBlock size is " + gerpBlockMap.size());
        System.out.println("alignedSegment size is " + alignedSegment.size());
        long lineFlag = 0;
        try {
            String temp = null;
            BufferedReader br = IOUtils.getTextReader(gerpDir);

            int gerpRangeIndex = 0;
            long GerpStartline = gerpBlockMap.get(alignedSegment.get(gerpRangeIndex)).get(0);
            long GerpEndline = gerpBlockMap.get(alignedSegment.get(gerpRangeIndex)).get(1);
            /*int GerpStartline = alignedBlockMap.get(alignedSegment.get(gerpRangeIndex)).get(0);
            int GerpEndline = alignedBlockMap.get(alignedSegment.get(gerpRangeIndex)).get(1);*/
            System.out.println("first gerpstartline is " + GerpStartline);

            System.out.println("first gerpEndline is " + GerpEndline);
            System.out.println("last gerpstartline is " + gerpBlockMap.get(alignedSegment.get(alignedSegment.size() - 1)).get(0));
            System.out.println("last gerpEndline is " + gerpBlockMap.get(alignedSegment.get(alignedSegment.size() - 1)).get(1));
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
                    while (gerpRangeIndex < alignedSegment.size()) {
                        if (gerpBlockMap.get(alignedSegment.get(gerpRangeIndex)).size() == 2) {
                            GerpStartline = gerpBlockMap.get(alignedSegment.get(gerpRangeIndex)).get(0);
                            GerpEndline = gerpBlockMap.get(alignedSegment.get(gerpRangeIndex)).get(1);
                            break;
                        } else {
                            System.out.println("curious gerpblock is " + gerpRangeIndex);
                            gerpRangeIndex++;

                        }
                    }

                }

            }
            br.close();
            System.out.println("gerpRange numbers is " + gerpRangeIndex);

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("wheat" + subgenome + " gerp score covers " + String.valueOf(lineFlag) + " bp");// the gerp.rates file covers
        return gerpLineMap;
    }

    public void normalizeGerp(String outfile, String chromosomeNum, String subgenome) {

        HashMap<Integer, ArrayList<String>> gerpLineMap = this.gerpLineMap;
        HashMap<Integer, TIntArrayList> alignedBlockMap = this.alignedBlockMap;
        ArrayList<Integer> alignedSegment = this.alignedSegment;
        int writedLineNum = 1;
        long genomeSize = this.genomeSize;
        System.out.println("alignedBlocksMap size is " + alignedBlockMap.size());
        System.out.println("gerpLineMap size is " + gerpLineMap.size());

        try {
            BufferedWriter bw = YaoIOUtils.getTextWriter(outfile);
            for (int i = 0; i < alignedSegment.size(); i++) {

                int startFlag = alignedBlockMap.get(alignedSegment.get(i)).get(0);
                int endFlag = alignedBlockMap.get(alignedSegment.get(i)).get(1);
                while (writedLineNum < startFlag) {

                    bw.write("0");
                    bw.write("\t");
                    bw.write("0");
                    bw.newLine();
                    writedLineNum++;

                }

                while (writedLineNum >= startFlag && writedLineNum <= endFlag) {
                    for (int j = writedLineNum-startFlag; j < gerpLineMap.get(alignedSegment.get(i)).size(); j++) {
                        bw.write(gerpLineMap.get(alignedSegment.get(i)).get(j));
                        bw.newLine();
                        writedLineNum++;
                    }
                }
            }
            while (writedLineNum <= genomeSize) {
                writedLineNum++;
                bw.write("0");
                bw.write("\t");
                bw.write("0");
                bw.newLine();
            }

            System.out.println(writedLineNum - 1 + " normalized gerps in wheat" + subgenome + ".chr" + chromosomeNum);
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
        this.outfile = pLineList.get(2);
        this.seqFile = pLineList.get(3);
        this.chromosomeNum = pLineList.get(4);
        this.subgenome = pLineList.get(5);

    }

    public static void main(String[] args) {
        /*  String mafDir="/data1/home/lipeng/result/CNEr/axt/A/formal/merge/gerp/wheatA1.maf";
        String gerpDir="/data1/home/lipeng/result/CNEr/axt/A/formal/merge/gerp/Triticeacea.cpy.mfa.rates";
        String outfile="/data1/home/lipeng/result/CNEr/axt/A/formal/merge/gerp/wheatA.chr1.gerp";
        String seqFile="/data1/home/lipeng/database/wheat/smRef/A/chr01A.fa";
        String mafDir = "/Users/kanglipeng/Desktop/test.maf";
        String gerpDir = "/Users/kanglipeng/Desktop/testout.txt";
        String outfile = "/Users/kanglipeng/Desktop/test.out";
        String seqFile = "/Users/kanglipeng/Desktop/test.fa";*/
        gerpAnalysis a = new gerpAnalysis(args[0]);
        new gerpAnalysis(a.mafDir, a.gerpDir, a.outfile, a.seqFile, a.chromosomeNum, a.subgenome);

    }

}
