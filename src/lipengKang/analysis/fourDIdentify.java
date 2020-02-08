/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;

import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashIntIntMaps;
import pgl.infra.range.RangeValStr;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
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
import lipengKang.analysis.ParseCDSFasta;

/**
 *
 * @author kanglipeng
 */
public class fourDIdentify {

    String mafDir = null;
    String gff3Dir = null;
    String CDSfasta =null;
    int phylipNum;

String outfile =null;
    HashMap<String, TIntArrayList> CDSPosMap = new HashMap<>();
    HashMap<Integer, TIntArrayList> alignedPosMap = new HashMap<>();
    TIntArrayList CDSPos = new TIntArrayList();
    ArrayList <String>CDSID = new ArrayList();
    ArrayList <Integer> alignedSegment = new ArrayList();
    TIntArrayList alignedPos = new TIntArrayList();
    public fourDIdentify(String parameterFileS) {
        this.parseParameters(parameterFileS);
    }
    public fourDIdentify(String gff3Dir, String mafDir) {
        getCDSPosMap(gff3Dir);
        getAlignedPosMap(mafDir);
    }  
    public String getoutfile(){
    return this.outfile;
    }
     public String getCDSfasta(){
    return this.CDSfasta;
    }
      public String getMaffile(){
    return this.mafDir;
    }
     public int getPhylipNum(){
    return this.phylipNum;
    } 
    public  ArrayList<String> getCDSID(){
    return this.CDSID;
    }
     public  ArrayList<Integer> getAlignedSegment(){
    return this.alignedSegment;
    }

    public HashMap<String, TIntArrayList> getCDSPosMap(String gff3Dir) {

        try {
            BufferedReader br = IOUtils.getTextReader(gff3Dir);
            String temp = null;
            String[] tem = null;
            String te[] = null;
            String CDS = null;

            while ((temp = br.readLine()) != null) {
                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);

                if (!tem[0].startsWith("chr")) {
                    continue;
                }
                if (!tem[0].contains("A")) {
                    continue;
                }
                if (!tem[2].startsWith("CDS")) {
                    continue;
                }
//split geneID
                te = tem[8].split(";");
                CDS = te[0].substring(3);
                int min = Integer.parseInt(tem[3]);
                int max = Integer.parseInt(tem[4]);
                CDSPos = new TIntArrayList();
                CDSPos.add(min);
                CDSPos.add(max);
                CDSID.add(CDS);
                CDSPosMap.put(CDS, CDSPos);
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("CDSPosMap established");
        return CDSPosMap;
//CDS position from gff3     TraesCS1A02G001900.6.cds1----->list(start,end)
    }
    
    public HashMap<Integer, TIntArrayList> getAlignedPosMap(String mafDir) {

        try{
            BufferedReader br = IOUtils.getTextReader(mafDir);
            String temp = null;
            String[] tem = null;
            int alignedSegmentNum = 0;

            while ((temp = br.readLine()) != null) {

                List<String> tList = KStringUtils.fastSplit(temp);
                List<String> tListNew = new ArrayList<>();
                if (!temp.startsWith("s wheatA.chr")){ 
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

                alignedSegmentNum++;
                int alignedSeqSize = tem[6].length();

               // if (tem[4].contains("+")) {
                    alignedPos = new TIntArrayList();

                    int startPos = Integer.parseInt(tem[2])+1;//maf produced by last aligner is 0-based, here we change to 1-based as CDS.fasta and gff3
                    int endPos = startPos + alignedSeqSize-1;
                    alignedPos.add(startPos);
                    alignedPos.add(endPos);

               /* } else {
                    alignedPos = new TIntArrayList();
                    int endPos = Integer.parseInt(tem[2]);
                    int startPos = endPos - alignedSeqSize+1;
                    alignedPos.add(startPos);
                    alignedPos.add(endPos);
                }*/
StringBuilder blockNum = new StringBuilder();
blockNum.append(tem[1].charAt(10)).append(alignedSegmentNum);
                alignedSegment.add(Integer.parseInt(blockNum.toString()));
                alignedPosMap.put(Integer.parseInt(blockNum.toString()), alignedPos);

        }

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
System.out.println("alignedPosMap established");
        return alignedPosMap;//postion map from maf   11------list(start,end) the first index is chromosome NUM , the second one is aligned block NUM
        
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
        this.CDSfasta=pLineList.get(2);
        this.outfile= pLineList.get(3);
        this.phylipNum=Integer.parseInt(pLineList.get(4));
        
    }

    public static void main(String[] args) {
        fourDIdentify la = new fourDIdentify(args[0]);
        fourDIdentify le = new fourDIdentify(la.gff3Dir, la.mafDir);
        new ParseCDSFasta(le,la);

    }

}
