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
import lipengKang.analysis.KStringUtils;
import zhouyao.analysis.wheatHapMap.YaoIOUtils;

/**
 *
 * @author kanglipeng
 */
public class gerpAnalysis {

    /* TIntArrayList chr1AlignedRange = new TIntArrayList();
    TIntArrayList chr2AlignedRange = new TIntArrayList();
    TIntArrayList chr3AlignedRange = new TIntArrayList();
    TIntArrayList chr4AlignedRange = new TIntArrayList();
    TIntArrayList chr5AlignedRange = new TIntArrayList();
    TIntArrayList chr6AlignedRange = new TIntArrayList();
    TIntArrayList chr7AlignedRange = new TIntArrayList();*/
    HashMap<Integer, TIntArrayList> alignedBlockMap = new HashMap<>();
    HashMap<Integer, TIntArrayList> gerpBlockMap = new HashMap<>();
    TIntArrayList alignedPos = new TIntArrayList();
    ArrayList<Integer> alignedSegment = new ArrayList();
    ArrayList<String> gerpLine = new ArrayList();
    TIntArrayList gerpRange = new TIntArrayList();
    HashMap<Integer, ArrayList<String>> gerpLineMap = new HashMap<>();
    ArrayList<String> gerpPos = new ArrayList();
    int genomeSize = 0;

    public gerpAnalysis(String mafDir, String gerpDir, String outfile, String seqFile) {
        this.countGenomeSize(seqFile);
        this.getAlignedRange(mafDir);
        this.parseGerpLine(gerpDir);
        this.normalizeGerp(outfile);

    }

    public void countGenomeSize(String seqFile) {
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
        System.out.println("genome size is " + String.valueOf(genomeSize) + "bp");

    }

    public HashMap<Integer, TIntArrayList> getAlignedRange(String mafDir) {
        int alignedSegmentNum = 0;
        int gerpStartPos = 1;
        int gerpEndPos = 0;

        try {
            BufferedReader br = IOUtils.getTextReader(mafDir);
            String temp = null;
            String[] tem = null;

            while ((temp = br.readLine()) != null) {

                List<String> tList = KStringUtils.fastSplit(temp);
                List<String> tListNew = new ArrayList<>();
                if (!temp.startsWith("s wheatA.chr")) {
                    continue;
                }
                if (!temp.startsWith("s wheatA.chr1")) {
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

            /*     if (tem[1].substring(10).equals("1")) {

                    chr1AlignedRange.add(startPos);
                    chr1AlignedRange.add(endPos);
                }
                if (tem[1].substring(10).equals("2")) {

                    chr2AlignedRange.add(startPos);
                    chr2AlignedRange.add(endPos);
                }
                if (tem[1].substring(10).equals("3")) {
                    chr3AlignedRange.add(startPos);
                    chr3AlignedRange.add(endPos);
                }
                if (tem[1].substring(10).equals("4")) {
                    chr4AlignedRange.add(startPos);
                    chr4AlignedRange.add(endPos);
                }
                if (tem[1].substring(10).equals("5")) {
                    chr5AlignedRange.add(startPos);
                    chr5AlignedRange.add(endPos);
                }
                if (tem[1].substring(10).equals("6")) {
                    chr6AlignedRange.add(startPos);
                    chr6AlignedRange.add(endPos);
                }
                if (tem[1].substring(10).equals("7")) {
                    chr7AlignedRange.add(startPos);
                    chr7AlignedRange.add(endPos);
                }

            }
            chr7AlignedRange.sort();
            chr6AlignedRange.sort();
            chr5AlignedRange.sort();
            chr4AlignedRange.sort();
            chr3AlignedRange.sort();
            chr2AlignedRange.sort();
            chr1AlignedRange.sort();

            alignedBlockMap.put(1, chr1AlignedRange);
            alignedBlockMap.put(2, chr2AlignedRange);
            alignedBlockMap.put(3, chr3AlignedRange);
            alignedBlockMap.put(4, chr4AlignedRange);
            alignedBlockMap.put(5, chr5AlignedRange);
            alignedBlockMap.put(6, chr6AlignedRange);
            alignedBlockMap.put(7, chr7AlignedRange);*/
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("alignedBlock established");
        return alignedBlockMap;//postion map from maf   1------list(start,end)  key is block number in gerpcol result;value is 1-based position in chromosome

    }

    public HashMap<Integer, ArrayList<String>> parseGerpLine(String gerpDir) {

        HashMap<Integer, TIntArrayList> gerpBlockMap = this.gerpBlockMap;
        ArrayList<Integer> alignedSegment = this.alignedSegment;
        long lineFlag = 0;
        try {
            String temp = null;
            BufferedReader br = IOUtils.getTextReader(gerpDir);

            int gerpRangeIndex = 0;
            int GerpStartline = gerpBlockMap.get(alignedSegment.get(gerpRangeIndex)).get(0);
            int GerpEndline = gerpBlockMap.get(alignedSegment.get(gerpRangeIndex)).get(1);
            /*int GerpStartline = alignedBlockMap.get(alignedSegment.get(gerpRangeIndex)).get(0);
            int GerpEndline = alignedBlockMap.get(alignedSegment.get(gerpRangeIndex)).get(1);*/
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

                /*  StringBuilder gerpChrPos = new StringBuilder();
                   gerpChrPos.append(String.valueOf(gerpChrStartPos)).append("-").append(String.valueOf(gerpChrEndPos));
                   gerpPos.add(gerpChrPos.toString());
 
                    GerpStartline = gerpRange.get(gerpRangeIndex);
                   GerpEndline = gerpRange.get(gerpRangeIndex + 1);*/
            }
            System.out.println("gerpRange numbers is " + gerpRangeIndex);

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("All gerp score covers " + String.valueOf(lineFlag));
        return gerpLineMap;
    }

    public void normalizeGerp(String outfile) {
        HashMap<Integer, ArrayList<String>> gerpLineMap = this.gerpLineMap;
        HashMap<Integer, TIntArrayList> alignedBlockMap = this.alignedBlockMap;
        ArrayList<Integer> alignedSegment = this.alignedSegment;
        int writedLineNum = 1;
        int genomeSize = this.genomeSize;
        System.out.println("aligned Blocks number is " + alignedSegment.size());
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

            System.out.println(writedLineNum - 1 + " gerps in the chromosome ");
            bw.flush();
            bw.close();

        } catch (Exception e) {
            System.out.println("Error in Reading and Write MSA");
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        String mafDir="/data1/home/lipeng/result/CNEr/axt/A/formal/merge/gerp/wheatA1.maf";
        String gerpDir="/data1/home/lipeng/result/CNEr/axt/A/formal/merge/gerp/Triticeacea.cpy.mfa.rates";
        String outfile="/data1/home/lipeng/result/CNEr/axt/A/formal/merge/gerp/wheatA.chr1.gerp";
        String seqFile="/data1/home/lipeng/database/wheat/smRef/A/chr01A.fa";
       /* String mafDir = "/Users/kanglipeng/Desktop/test.maf";
        String gerpDir = "/Users/kanglipeng/Desktop/testout.txt";
        String outfile = "/Users/kanglipeng/Desktop/test.out";
        String seqFile = "/Users/kanglipeng/Desktop/test.fa";*/
        new gerpAnalysis(mafDir, gerpDir, outfile, seqFile);
    }

}
