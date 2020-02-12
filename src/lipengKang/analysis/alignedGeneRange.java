/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;

import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import lipengKang.analysis.KStringUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author kanglipeng
 */
public class alignedGeneRange {

    String mafDir = null;
    String gff3Dir = null;

    HashMap<String, TIntArrayList> genePosMap = new HashMap<>();
    HashMap<String, TIntArrayList> alignedPosMap = new HashMap<>();
    TIntArrayList genePos = new TIntArrayList();

    List gene = new ArrayList();
    List alignedSegment = new ArrayList();
    TIntArrayList alignedPos = new TIntArrayList();
    TIntArrayList coverNum = new TIntArrayList();

    public alignedGeneRange(String parameterFileS) {
        this.parseParameters(parameterFileS);
    }

    public alignedGeneRange(String gff3Dir, String mafDir) {

        getGenePosMap(gff3Dir);
        getAlignedPosMap(mafDir);

    }

    public alignedGeneRange(alignedGeneRange li, String gff3Dir, String mafDir) {
        geneMatchPercent(li, gff3Dir, mafDir);
    }

    public HashMap<String, TIntArrayList> getGenePosMap(String gff3Dir) {

        try {
            BufferedReader br = IOUtils.getTextReader(gff3Dir);
            String temp = null;
            String[] tem = null;
            String te[] = null;
            String geneID = null;

            while ((temp = br.readLine()) != null) {
                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);

                if (!tem[0].startsWith("chr")) {
                    continue;
                }
                if (!tem[0].contains("A")) {
                    continue;
                }
                if (!tem[2].startsWith("gene")) {
                    continue;
                }
//split geneID
                te = tem[8].split(";");
                geneID = te[0];
                int min = Integer.parseInt(tem[3]);
                int max = Integer.parseInt(tem[4]);
                genePos = new TIntArrayList();
                genePos.add(min);
                genePos.add(max);

                gene.add(tem[0] + geneID);
                genePosMap.put(tem[0] + geneID, genePos);

            }
            

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        return genePosMap;

    }


    public HashMap<String, TIntArrayList> getAlignedPosMap(String mafDir) {

        try {
            BufferedReader br = IOUtils.getTextReader(mafDir);
            String temp = null;
            String[] tem = null;
            int alignedSegmentNum = 0;

            while ((temp = br.readLine()) != null) {

                List<String> tList = KStringUtils.fastSplit(temp);
                List<String> tListNew = new ArrayList<>();
                if (!temp.startsWith("s wheatA-chr")) {
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
                /* if (!tem[0].startsWith("s")) {
                        continue;
                    }
                    if (!tem[1].startsWith("wheatA.chr")) {
                        continue;
                    }*/
                alignedSegmentNum++;
                int alignedSeqSize = tem[6].length();

                if (tem[4].contains("+")) {
                    alignedPos = new TIntArrayList();

                    int startPos = Integer.parseInt(tem[2]);
                    int endPos = startPos + alignedSeqSize;
                    alignedPos.add(startPos);
                    alignedPos.add(endPos);

                } else {
                    alignedPos = new TIntArrayList();
                    int endPos = Integer.parseInt(tem[2]);
                    int startPos = endPos - alignedSeqSize;
                    alignedPos.add(startPos);
                    alignedPos.add(endPos);
                }

                alignedSegment.add(tem[1] + alignedSegmentNum);
                alignedPosMap.put(tem[1] + alignedSegmentNum, alignedPos);

            }

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        return alignedPosMap;
    }

    public void geneMatchPercent(alignedGeneRange li, String gff3Dir, String mafDir) {
        File maf = new File(mafDir);
        String mafname = maf.getName();

        HashMap<String, TIntArrayList> genePosMap = li.getGenePosMap(gff3Dir);
        gene = li.gene;
Set s=new HashSet(gene);
gene.clear();
gene.addAll(s);
        HashMap<String, TIntArrayList> alignedPosMap = li.getAlignedPosMap(mafDir);
        alignedSegment = li.alignedSegment;
Set x=new HashSet(alignedSegment);
alignedSegment.clear();
alignedSegment.addAll(x);
        int alignedStart = 0;
        int alignedEnd = 0;
        int geneStart = 0;
        int geneEnd = 0;
        int coverBaseNum = 0;
        int allCover = 0;
        int geneBaseNum = 0;
        for (int i = 0; i < alignedSegment.size(); i++) {

            for (int j = 0; j < gene.size(); j++) {
                for (int m = 1; m <= 7; m++) {
                    if (alignedSegment.get(i).toString().startsWith("wheatA-chr" + m) && gene.get(j).toString().startsWith("chr" + m + "A")) {
                        alignedStart = alignedPosMap.get(alignedSegment.get(i)).get(0);
                        alignedEnd = alignedPosMap.get(alignedSegment.get(i)).get(1);

                        // if (gene.get(j).toString().startsWith("chr"+m+"A")){
                        geneStart = genePosMap.get(gene.get(j)).get(0);
                        geneEnd = genePosMap.get(gene.get(j)).get(1);

                        if (alignedEnd < geneStart || alignedStart > geneEnd) {
                            continue;
                        }
                        if (alignedEnd <= geneEnd && alignedStart >= geneStart) {
                            coverBaseNum = alignedEnd - alignedStart;
                        }
                        if (alignedStart <= geneStart && alignedEnd <= geneEnd) {
                            coverBaseNum = alignedEnd - geneStart;
                        }
                        if (alignedStart >= geneStart && alignedEnd >= geneEnd) {
                            coverBaseNum = geneEnd - alignedStart;
                        }
                        coverNum.add(coverBaseNum);
                        allCover = allCover + coverBaseNum;
                    }
                }

            }
        }
        for (int j = 0; j < gene.size(); j++) {

            geneStart = genePosMap.get(gene.get(j)).get(0);
            geneEnd = genePosMap.get(gene.get(j)).get(1);
            geneBaseNum = geneBaseNum + geneEnd - geneStart;

        }
 DecimalFormat df=new DecimalFormat("0.00");
        // int geneMatchBaseNum = geneMatchPosSet.size();
        System.out.println(mafname + " maf covers " + 100 * (float)allCover/geneBaseNum + "% bases of gene region");
        //geneMatchPosSet.clear();
        //geneMatchPosSet.addAll(geneUpDown2kbPosSet);
        //geneMatchPosSet.retainAll(alignedPosSet);
        //int geneUpDown2kbMatchBaseNum = geneMatchPosSet.size();
        //System.out.println(mafname + " maf covers " + 100 * geneUpDown2kbMatchBaseNum / geneUpDown2kbBaseNum + "% bases of gene up and down stream 2kb region");
        // geneMatchPosSet.clear();
        //geneMatchPosSet.addAll(geneUpDown5kbPosSet);
        // geneMatchPosSet.retainAll(alignedPosSet);
        // int geneUpDown5kbMatchBaseNum = geneMatchPosSet.size();
        //System.out.println(mafname + " maf covers " + 100 * geneUpDown5kbMatchBaseNum / geneUpDown5kbBaseNum + "% bases of gene up and down stream 5kb region");

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
        this.gff3Dir = pLineList.get(0);
        this.mafDir = pLineList.get(1);

    }

    public static void main(String[] args) {
        alignedGeneRange le = new alignedGeneRange(args[0]);
        //  String gff3Dir="/Users/kanglipeng/Desktop/wheat.gff3";
        //String mafDir="/Users/kanglipeng/Desktop/wheatA-stiffbrome.maf";
        alignedGeneRange li = new alignedGeneRange(le.gff3Dir, le.mafDir);
        new alignedGeneRange(li, le.gff3Dir, le.mafDir);

    }

}
