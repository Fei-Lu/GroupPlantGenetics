/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import pgl.infra.utils.IOUtils;
import lipengKang.analysis.KStringUtils;
import pgl.infra.utils.PStringUtils;
import zhouyao.analysis.wheatHapMap.YaoIOUtils;
import lipengKang.analysis.fourDIdentify;
import static lipengKang.analysis.RandomArray.randomCommon;

/**
 *
 * @author kanglipeng
 */
public class ParseCDSFasta {

    HashMap<String, String> CDSMap = new HashMap<>();
    StringBuilder seq = new StringBuilder();
    ArrayList<String> mRnaList = new ArrayList<>();
    HashMap<String, TIntArrayList> fourDPos = new HashMap<>();
    HashMap<String, TIntArrayList> mRnaCDSPos = new HashMap<>();
    HashMap<String,Character > CDSstrandMap = new HashMap<>();
    TIntArrayList fourDList = new TIntArrayList();
    TIntArrayList CDSPosList = new TIntArrayList();
    ArrayList<String> CdsID = new ArrayList<>();
    HashMap<String, TIntArrayList> finalFourDMap = new HashMap<>();
    HashMap<String, TIntArrayList> fourDPosMap = new HashMap<>();
    HashMap<Integer, TIntArrayList> alignedFourDMap = new HashMap<>();
    HashMap<String, String> filteredCDSMap = new HashMap();
    ArrayList<String> filteredmRnaList = new ArrayList();
    HashMap<Integer, ArrayList<Integer>> chrFourDMap = new HashMap<>();
    HashMap<Integer, ArrayList<Integer>> alignedChrFourDMap = new HashMap<>();
    ArrayList<String> Cds = new ArrayList();
    ArrayList <Integer>fourDBlock = new ArrayList();
    fourDIdentify le = null;
    fourDIdentify la = null;

    public ParseCDSFasta(fourDIdentify le, fourDIdentify la) {
        this.la = la;
        this.le = le;
        this.getCDSSeq(la);
        this.parseFourDPos();
        this.mergeFourDPos();
        this.parseFourDPosMap(le);
        this.getAlignedChrFourDMap(le);
        this.extractMSAFourD(la);
    }

    public HashMap<String, String> getCDSSeq(fourDIdentify la) {
        try {
            BufferedReader br;
            String CDSfasta = la.getCDSfasta();
            if (CDSfasta.endsWith("gz")) {
                br = YaoIOUtils.getTextGzipReader(CDSfasta);

            } else {
                br = YaoIOUtils.getTextReader(CDSfasta);
            }
            String temp = null;
            String[] tem = null;

            while ((temp = br.readLine()) != null) {
                if (temp.startsWith(">") && temp.substring(9).startsWith("A")) {
                    List<String> tList = KStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                    mRnaList.add(tem[0].substring(1));
                    List<String> kList = KStringUtils.fastSplitComma(tem[3].substring(5));
                    char CDSstrand =tem[2].charAt(10);
                    for (int i = 0; i < kList.size(); i++) {
                        for (int j = 0; j < 2; j++) {
                            CDSPosList.add(Integer.parseInt(KStringUtils.fastSplitStrigula(kList.get(i)).get(j)));
//>TraesCS1A02G001900.6 gene=TraesCS1A02G001900 loc:chr1A(-)1159194-1162777 segs:1-46,47-249,250-447,448-625,626-796,797-828,829-871,872-939,940-978,979-1431
//1-46                        
                        }
                        if (CDSstrand=='+'){
                        CdsID.add(tem[0].substring(1) + ".cds" + (i + 1));
                        mRnaCDSPos.put(CdsID.get(CdsID.size() - 1), CDSPosList);
                        CDSstrandMap.put(CdsID.get(CdsID.size() - 1), CDSstrand);
                        //mRnaCDSPos  refers all cds position of its mRNA     TraesCS1A02G001900.6.cds1------>list(start,end)
                        CDSPosList = new TIntArrayList();}else{
                        CdsID.add(tem[0].substring(1) + ".cds" + (kList.size() - i));
                        mRnaCDSPos.put(CdsID.get(CdsID.size() - 1), CDSPosList);
                        CDSstrandMap.put(CdsID.get(CdsID.size() - 1), CDSstrand);
                         CDSPosList = new TIntArrayList();
                        }
                    }
                    if (temp.startsWith(">TraesCS1A02G000100.1")) {
                        seq = new StringBuilder();
                    } else {
                        CDSMap.put(mRnaList.get(mRnaList.size() - 2), seq.toString());
                        seq = new StringBuilder();
                    }
                } else {
                    seq.append(temp);
                    continue;
                }
            }
            CDSMap.put(mRnaList.get(mRnaList.size() - 1), seq.toString());
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("CDSMap established");
        return CDSMap;//CDS sequence from CDS.fasta   TraesCS1A02G000100.1------>mRNAseq 
    }

    public HashMap<String, TIntArrayList> parseFourDPos() {
        HashMap<String, String> CDSMap = this.CDSMap;
        ArrayList<String> mRnaList = this.mRnaList;
        String CDSSeq = null;
        ArrayList<String> fourDCodons = new ArrayList<String>(Arrays.asList("CT", "GT", "TC", "CC", "AC", "GC", "CG", "GG"));
        for (String i : mRnaList) {
            if (i.endsWith(".1")) {// filterd mrna except the first mrna
                filteredmRnaList.add(i);
                filteredCDSMap.put(i, CDSMap.get(i));
            }
        }
        for (int i = 0; i < filteredmRnaList.size(); i++) {
            CDSSeq = filteredCDSMap.get(filteredmRnaList.get(i));
            for (int j = 0; j < CDSSeq.length() - 2; j += 3) {
                if (fourDCodons.contains(CDSSeq.substring(j, j + 2))) {
                    fourDList.add(j + 3);
                }
            }
            fourDPos.put(filteredmRnaList.get(i), fourDList);
            fourDList = new TIntArrayList();
        }
        System.out.println("fourDPos established");
        return fourDPos;//four-fold degenerate sites position according to each mRNA     TraesCS1A02G000100.1------->(position--3,6,9...) 1-based    
    }

    public HashMap<String, TIntArrayList> mergeFourDPos() {
        HashMap<String, TIntArrayList> fourDPos = this.fourDPos;
        HashMap<String, TIntArrayList> mRnaCDSPos = this.mRnaCDSPos;
        ArrayList<String> CdsID = this.CdsID;
        ArrayList<String> filteredmRnaList = this.filteredmRnaList;
        TIntArrayList finalFourD = new TIntArrayList();
        TIntArrayList fourD = new TIntArrayList();
        TIntArrayList mRnaCDS = new TIntArrayList();
        String[] tep = null;

        for (int j = 0; j < filteredmRnaList.size(); j++) {
            for (int i = 0; i < CdsID.size(); i++) {
                List<String> tList = KStringUtils.fastSplitDot(CdsID.get(i));
                tep = tList.toArray(new String[tList.size()]);
                StringBuilder temp = new StringBuilder();
                temp = temp.append(tep[0]).append(".").append(tep[1]);
                if (temp.toString().equals(filteredmRnaList.get(j))) {
                    fourD = fourDPos.get(filteredmRnaList.get(j));
                    mRnaCDS = mRnaCDSPos.get(CdsID.get(i));
                    for (int k = 0; k < fourD.size(); k++) {
                        if (mRnaCDS.get(0) <= fourD.get(k) && fourD.get(k) <= mRnaCDS.get(1)) {
                            finalFourD.add(fourD.get(k) - mRnaCDS.get(0));
                        }
                    }
                    Cds.add(CdsID.get(i));
                    finalFourDMap.put(CdsID.get(i), finalFourD);

                    finalFourD = new TIntArrayList();
                }
            }
        }
        System.out.println("finalFourDMap established");
        return finalFourDMap;// TraesCS1A02G001900.6.cds1------>list(position,---2,5)     relative Positions according to each cds( 0-based)
    }


    /* public HashMap<String, TIntArrayList> parseFourDPosMap(fourDIdentify le) {
        HashMap<String, TIntArrayList> CDSPosMap = le.CDSPosMap;
        ArrayList<String> CDSID = le.getCDSID();
        HashMap<String, TIntArrayList> finalFourDMap = this.finalFourDMap;
        ArrayList<String> Cds = this.Cds;
        TIntArrayList fourDSite = new TIntArrayList();

        int CDSStartPos = 0;
        int fourDRelativePos = 0;
        for (int i = 0; i < Cds.size(); i++) {
            for (int j = 0; j < CDSID.size(); j++) {
                if (CDSID.get(j).equals(Cds.get(i))) {
                    CDSStartPos = CDSPosMap.get(CDSID.get(j)).get(0);
                    for (int k = 0; k < finalFourDMap.get(Cds.get(i)).size(); k++) {
                        fourDRelativePos = finalFourDMap.get(Cds.get(i)).get(k);
                        fourDSite.add(CDSStartPos + fourDRelativePos);
                    }
                    fourDPosMap.put(Cds.get(i), fourDSite);
                    fourDSite = new TIntArrayList();
                }
            }
        }
        System.out.println("fourDPosMap established");
        return fourDPosMap;//reference based position of fourfold degenerate site      TraesCS1A02G001900.6.cds1------>list(3,9...)
    }*/
    public HashMap<Integer, ArrayList<Integer>> parseFourDPosMap(fourDIdentify le) {
        HashMap<String, TIntArrayList> CDSPosMap = le.CDSPosMap;
        ArrayList<String> CDSID = le.getCDSID();
        HashMap<String, TIntArrayList> finalFourDMap = this.finalFourDMap;
        HashMap<String, Character> CDSstrandMap = this.CDSstrandMap;
        ArrayList<String> Cds = this.Cds;
        ArrayList<Integer> fourDSite = new ArrayList();
       
 
        int CDSStartPos = 0;
        int fourDRelativePos = 0;
        for (int l = 1; l <= 7; l++) {

            for (int i = 0; i < Cds.size(); i++) {
                for (int j = 0; j < CDSID.size(); j++) {

                    if (CDSID.get(j).substring(7).startsWith(String.valueOf(l))&&CDSID.get(j).equals(Cds.get(i))) {
                        if(CDSstrandMap.get(Cds.get(i))=='+'){
                        CDSStartPos = CDSPosMap.get(CDSID.get(j)).get(0);
                        for (int k = 0; k < finalFourDMap.get(Cds.get(i)).size(); k++) {
                            fourDRelativePos = finalFourDMap.get(Cds.get(i)).get(k);
                            fourDSite.add(CDSStartPos + fourDRelativePos);

                        }}else{
                            CDSStartPos = CDSPosMap.get(CDSID.get(j)).get(1);
                        for (int k = 0; k < finalFourDMap.get(Cds.get(i)).size(); k++) {
                            fourDRelativePos = finalFourDMap.get(Cds.get(i)).get(k);
                            fourDSite.add(CDSStartPos - fourDRelativePos);

                        }
                        
                        }
                     
                    }
                }
            }
            chrFourDMap.put(l, fourDSite);
            fourDSite = new ArrayList();
        }
        System.out.println("fourDPosMap established");//1-based
        return chrFourDMap;//reference based position of fourfold degenerate site      TraesCS1A02G001900.6.cds1------>list(3,9...)
    }
    //public HashMap<Integer, TIntArrayList> getRandomFourDPosMap(){
    //HashMap<String, TIntArrayList> fourDPosMap = this.fourDPosMap;
    // ArrayList<String> CdsID = this.CdsID;
    //  for(String i:CdsID){
    //if (KStringUtils.fastSplitDot(i).get(0).equals)

    // }
    //}
    public HashMap<Integer, ArrayList<Integer>>  getAlignedChrFourDMap(fourDIdentify le) {
       
        HashMap<Integer, ArrayList<Integer>> chrFourDMap = this.chrFourDMap;
        HashMap<Integer, TIntArrayList> alignedPosMap = le.alignedPosMap;
    ArrayList<Integer> alignedSegment = le.getAlignedSegment();
       ArrayList fourDSite = new ArrayList();
       ArrayList <Integer>alignedFourDSite = new ArrayList();
       
    for(int i=1;i<=7;i++){
        fourDSite=chrFourDMap.get(i);
        Collections.sort(fourDSite);
    for(int j=0;j<alignedSegment.size();j++){
        if(String.valueOf(alignedSegment.get(j)).startsWith(String.valueOf(i))){
         int blockStart=alignedPosMap.get(alignedSegment.get(j)).get(0);
        int blockEnd= alignedPosMap.get(alignedSegment.get(j)).get(1);
        for(int k=blockStart;k<=blockEnd;k++){
      if( Collections.binarySearch(fourDSite, k)>=0){
      alignedFourDSite.add(k-blockStart);
      
      }
        }
        alignedChrFourDMap.put(alignedSegment.get(j), alignedFourDSite);
        fourDBlock.add(alignedSegment.get(j));
        alignedFourDSite = new ArrayList();
        }
    
    }
    
    
    }

        return alignedChrFourDMap;
    }

   /* public HashMap<Integer, TIntArrayList> getAlignedFourDMap(fourDIdentify le) {
        HashMap<String, TIntArrayList> fourDPosMap = this.fourDPosMap;
        ArrayList<String> Cds = this.Cds;
        TIntArrayList alignedFourDSite = new TIntArrayList();
        ArrayList<Integer> alignedSegment = le.getAlignedSegment();
        HashMap<Integer, TIntArrayList> alignedPosMap = le.alignedPosMap;
        TIntArrayList fourDSites = new TIntArrayList();
        TIntArrayList alignedPos = new TIntArrayList();
        String[] tep = null;
        for (int i = 0; i < alignedSegment.size(); i++) {
            for (int j = 0; j < Cds.size(); j++) {
                List<String> tList = KStringUtils.fastSplitDot(Cds.get(j));
                tep = tList.toArray(new String[tList.size()]);
                if (alignedSegment.get(i).toString().startsWith(String.valueOf(tep[0].charAt(7)))) {
                    // if (CdsID.get(j).substring(8).startsWith("A")) {
                    fourDSites = fourDPosMap.get(Cds.get(j));
                    alignedPos = alignedPosMap.get(alignedSegment.get(i));
                    if (fourDSites == null || fourDSites.size() == 0) {
                        continue;
                    }
                    if (fourDSites.size() == 1) {
                        if (alignedPos.get(0) <= fourDSites.get(0) && fourDSites.get(0) <= alignedPos.get(1)) {
                            alignedFourDSite.add(fourDSites.get(0) - alignedPos.get(0));
                        }
                    }
                    if (fourDSites.size() >= 2) {
                        if (fourDSites.get(0) > alignedPos.get(1) || fourDSites.get(fourDSites.size() - 1) < alignedPos.get(0)) {
                            continue;
                        }
                        for (int k = 0; k < fourDSites.size(); k++) {
                            if (alignedPos.get(0) <= fourDSites.get(k) && fourDSites.get(k) <= alignedPos.get(1)) {
                                alignedFourDSite.add(fourDSites.get(k) - alignedPos.get(0));
                            }

                        }
                    }
                    // }
                }
            }

            alignedFourDMap.put(alignedSegment.get(i), alignedFourDSite);
            alignedFourDSite = new TIntArrayList();
        }
        System.out.println("alignedFourDMap established");
        return alignedFourDMap;//1------->(position--which base is fourD in each block)
    }*/

    public void extractMSAFourD(fourDIdentify la) {
        HashMap<Integer, ArrayList<Integer>> finalAlignedFourDMap = this.alignedChrFourDMap;
        ArrayList <Integer>fourDBlock = this.fourDBlock;
        // ArrayList<Integer> alignedSegment = le.getAlignedSegment();
        //ArrayList<Integer> alignedFourDList = new ArrayList();
        // HashMap<Integer, ArrayList<Integer>> finalAlignedFourDMap = new HashMap();

        //filter same position from different mrna
        /* for (int i : alignedSegment) {
            for (int j = 0; j < alignedFourDMap.get(i).size(); j++) {
                if (!alignedFourDList.contains(alignedFourDMap.get(i).get(j))) {
                    alignedFourDList.add(alignedFourDMap.get(i).get(j));

                }
            }

            finalAlignedFourDMap.put(i, alignedFourDList);
            alignedFourDList = new ArrayList();
        }*/
 /* for (int i : alignedSegment) {
          alignedFourDList = new ArrayList();
                alignedFourDSet = new HashSet();
           for (int j = 0; j < alignedFourDMap.get(i).size(); j++) {
               
               int x=alignedFourDMap.get(i).get(j);
               if(alignedFourDSet.add(x)){
               alignedFourDList.add(x);
               }
             
           }
        finalAlignedFourDMap.put(i, alignedFourDList);
       }*/
  int phylipNum=la.getPhylipNum();
        System.out.println("finalAlignedFourDMap established");
        BufferedWriter bw;
        BufferedReader br;
        String outfile = la.getoutfile();
        String mafDir = la.getMaffile();
        if (mafDir.endsWith("gz")) {
            br = YaoIOUtils.getTextGzipReader(mafDir);
            bw = YaoIOUtils.getTextGzipWriter(outfile);
        } else {
            br = YaoIOUtils.getTextReader(mafDir);
            bw = YaoIOUtils.getTextWriter(outfile);
        }
        String temp = null;
        String[] tem = null;
        int alignedNum = 0;//take care

        int alignedSegmentNum = 0;
        StringBuilder wheatA = new StringBuilder();
        StringBuilder randomWheatAFourD = new StringBuilder();
        StringBuilder wheatD = new StringBuilder();
        StringBuilder randomWheatDFourD = new StringBuilder();
        StringBuilder wheatB = new StringBuilder();
        StringBuilder randomWheatBFourD = new StringBuilder();
        StringBuilder barley = new StringBuilder();
        StringBuilder randomBarleyFourD = new StringBuilder();
        StringBuilder stiffbrome = new StringBuilder();
        StringBuilder randomStiffbromeFourD = new StringBuilder();
        StringBuilder rice = new StringBuilder();
        StringBuilder randomRiceFourD = new StringBuilder();
        StringBuilder maize = new StringBuilder();
        StringBuilder randomMaizeFourD = new StringBuilder();
        StringBuilder sorghum = new StringBuilder();
        StringBuilder randomSorghumFourD = new StringBuilder();
        System.out.println("Stringbuilder build successfully");
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("a")) {
                    alignedSegmentNum++;
                    StringBuilder blockNum = new StringBuilder();
                    
                    temp = br.readLine();
                    if (temp.startsWith("s wheatA")) {
                        blockNum.append(temp.charAt(12)).append(alignedSegmentNum);
                    alignedNum = Integer.parseInt(blockNum.toString());
                   /* Collections.sort(fourDBlock);
                    if( Collections.binarySearch(fourDBlock, alignedNum)>=0){*/

                        List<String> aList = KStringUtils.fastSplit(temp);
                        List<String> aListNew = new ArrayList<>();
                        for (int i = 0; i < aList.size(); i++) {

                            if (aList.get(i) != null && !aList.get(i).equals("")) {

                                aListNew.add(aList.get(i));
                            }
                        }
                        tem = aListNew.toArray(new String[aListNew.size()]);
                        for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                            wheatA.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                            
                        }
                    /*}else{
                    continue;

                    }*/
                    }
                    temp = br.readLine();

                    //wheatA-wheatD-
                    if (temp.startsWith("s wheatD")) {
                        List<String> bList = KStringUtils.fastSplit(temp);
                        List<String> bListNew = new ArrayList<>();
                        for (int i = 0; i < bList.size(); i++) {

                            if (bList.get(i) != null && !bList.get(i).equals("")) {

                                bListNew.add(bList.get(i));
                            }
                        }
                        tem = bListNew.toArray(new String[bListNew.size()]);
                        for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                            wheatD.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                        }
                        temp = br.readLine();
                        //wheatA-wheatD-wheatB-
                        if (temp.startsWith("s wheatB")) {
                            List<String> aList = KStringUtils.fastSplit(temp);
                            List<String> aListNew = new ArrayList<>();
                            for (int i = 0; i < aList.size(); i++) {

                                if (aList.get(i) != null && !aList.get(i).equals("")) {

                                    aListNew.add(aList.get(i));
                                }
                            }
                            tem = aListNew.toArray(new String[aListNew.size()]);
                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                wheatB.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                            }
                            temp = br.readLine();
                            //wheatA-wheatD-wheatB-barley-   
                            if (temp.startsWith("s barley")) {
                                List<String> cList = KStringUtils.fastSplit(temp);
                                List<String> cListNew = new ArrayList<>();
                                for (int i = 0; i < cList.size(); i++) {

                                    if (cList.get(i) != null && !cList.get(i).equals("")) {

                                        cListNew.add(cList.get(i));
                                    }
                                }
                                tem = cListNew.toArray(new String[cListNew.size()]);
                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                    barley.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                }
                                temp = br.readLine();
                                //wheatA-wheatD-wheatB-barley-stiffbrome-
                                if (temp.startsWith("s stiffbrome")) {
                                    List<String> dList = KStringUtils.fastSplit(temp);
                                    List<String> dListNew = new ArrayList<>();
                                    for (int i = 0; i < dList.size(); i++) {

                                        if (dList.get(i) != null && !dList.get(i).equals("")) {

                                            dListNew.add(dList.get(i));
                                        }
                                    }
                                    tem = dListNew.toArray(new String[dListNew.size()]);
                                    for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                        stiffbrome.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                    }
                                    temp = br.readLine();
                                    //wheatA-wheatD-wheatB-barley-stiffbrome-rice-
                                    if (temp.startsWith("s rice")) {
                                        List<String> eList = KStringUtils.fastSplit(temp);
                                        List<String> eListNew = new ArrayList<>();
                                        for (int i = 0; i < eList.size(); i++) {

                                            if (eList.get(i) != null && !eList.get(i).equals("")) {

                                                eListNew.add(eList.get(i));
                                            }
                                        }
                                        tem = eListNew.toArray(new String[eListNew.size()]);
                                        for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                            rice.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                        }
                                        temp = br.readLine();
                                        //wheatA-wheatD-wheatB-barley-stiffbrome-rice-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatD-wheatB-barley-stiffbrome-rice-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatD-wheatB-barley-stiffbrome-rice-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }

                                        } else {//wheatA-wheatD-wheatB-barley-stiffbrome-rice-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    maize.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatD-wheatB-barley-stiffbrome-rice

                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    maize.append('-');
                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        }

                                    } else {//wheatA-wheatD-wheatB-barley-stiffbrome-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);

                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatD-wheatB-barley-stiffbrome-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    rice.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatD-wheatB-barley-stiffbrome-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    rice.append('-');
                                                    sorghum.append('-');
                                                }
                                                continue;
                                            }

                                        } else {////wheatA-wheatD-wheatB-barley-stiffbrome-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {
                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);

                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    rice.append('-');
                                                    maize.append('-');
                                                }

                                                continue;

                                            } else {
                                                //wheatA-wheatD-wheatB-barley-stiffbrome
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    rice.append('-');
                                                    maize.append('-');
                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        }
                                    }

                                } else {//wheatA-wheatD-wheatB-barley-rice-
                                    if (temp.startsWith("s rice")) {
                                        List<String> eList = KStringUtils.fastSplit(temp);
                                        List<String> eListNew = new ArrayList<>();
                                        for (int i = 0; i < eList.size(); i++) {

                                            if (eList.get(i) != null && !eList.get(i).equals("")) {

                                                eListNew.add(eList.get(i));
                                            }
                                        }
                                        tem = eListNew.toArray(new String[eListNew.size()]);

                                        for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                            rice.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                        }
                                        temp = br.readLine();
                                        //wheatA-wheatD-wheatB-barley-rice-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);

                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatD-wheatB-barley-rice-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    stiffbrome.append('-');
                                                }

                                                continue;

                                            } else {    //wheatA-wheatD-wheatB-barley-rice-maize

                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    stiffbrome.append('-');
                                                    sorghum.append('-');
                                                }
                                                continue;
                                            }

                                        } else {//wheatA-wheatD-wheatB-barley-rice-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    stiffbrome.append('-');
                                                    maize.append('-');
                                                }

                                                continue;

                                            } else {//wheatA-wheatD-wheatB-barley-rice
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    stiffbrome.append('-');
                                                    maize.append('-');
                                                    sorghum.append('-');
                                                }
                                                continue;
                                            }
                                        }
                                    } else {////wheatA-wheatD-wheatB-barley-mazie-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatD-wheatB-barley-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    stiffbrome.append('-');

                                                    rice.append('-');
                                                }

                                                continue;

                                            } else {
                                                //wheatA-wheatD-wheatB-barley-mazie

                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    stiffbrome.append('-');
                                                    rice.append('-');
                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        } else {// //wheatA-wheatD-wheatB-barley-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    stiffbrome.append('-');
                                                    rice.append('-');
                                                    maize.append('-');
                                                }

                                                continue;

                                            } else {
                                                //wheatA-wheatD-wheatB-barley
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    stiffbrome.append('-');
                                                    rice.append('-');
                                                    maize.append('-');
                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        }
                                    }
                                }

                            } else {// //wheatA-wheatD-wheatB-stiffbrome-
                                if (temp.startsWith("s stiffbrome")) {
                                    List<String> dList = KStringUtils.fastSplit(temp);
                                    List<String> dListNew = new ArrayList<>();
                                    for (int i = 0; i < dList.size(); i++) {

                                        if (dList.get(i) != null && !dList.get(i).equals("")) {

                                            dListNew.add(dList.get(i));
                                        }
                                    }
                                    tem = dListNew.toArray(new String[dListNew.size()]);
                                    for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                        stiffbrome.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                    }
                                    temp = br.readLine();
                                    //wheatA-wheatD-wheatB-stiffbrome-rice-
                                    if (temp.startsWith("s rice")) {
                                        List<String> eList = KStringUtils.fastSplit(temp);
                                        List<String> eListNew = new ArrayList<>();
                                        for (int i = 0; i < eList.size(); i++) {

                                            if (eList.get(i) != null && !eList.get(i).equals("")) {

                                                eListNew.add(eList.get(i));
                                            }
                                        }
                                        tem = eListNew.toArray(new String[eListNew.size()]);
                                        for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                            rice.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                        }
                                        temp = br.readLine();
                                        //wheatA-wheatD-wheatB-stiffbrome-rice-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatD-wheatB-stiffbrome-rice-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    barley.append('-');
                                                }

                                                continue;

                                            } else {
                                                //wheatA-wheatD-wheatB-stiffbrome-rice-maize

                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    barley.append('-');
                                                    sorghum.append('-');
                                                }
                                                continue;
                                            }
                                        } else {  //wheatA-wheatD-wheatB-stiffbrome-rice-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    barley.append('-');
                                                    maize.append('-');
                                                }

                                                continue;

                                            } else {
                                                //wheatA-wheatD-wheatB-stiffbrome-rice
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    barley.append('-');
                                                    maize.append('-');
                                                    sorghum.append('-');
                                                }
                                                continue;
                                            }

                                        }
                                    } else { //wheatA-wheatD-wheatB-stiffbrome-maize
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatD-wheatB-stiffbrome-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    barley.append('-');
                                                    rice.append('-');
                                                }

                                                continue;

                                            } else {
                                                //wheatA-wheatD-wheatB-stiffbrome-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    barley.append('-');
                                                    rice.append('-');
                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }

                                        } else { //wheatA-wheatD-wheatB-stiffbrome-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    barley.append('-');
                                                    rice.append('-');
                                                    maize.append('-');
                                                }

                                                continue;

                                            } else {//wheatA-wheatD-wheatB-stiffbrome
                                                //!!!!!
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    barley.append('-');
                                                    rice.append('-');
                                                    maize.append('-');
                                                    sorghum.append('-');
                                                }
                                                continue;
                                            }

                                        }
                                    }
                                } else {//wheatA-wheatD-wheatB-rice-
                                    if (temp.startsWith("s rice")) {
                                        List<String> eList = KStringUtils.fastSplit(temp);
                                        List<String> eListNew = new ArrayList<>();
                                        for (int i = 0; i < eList.size(); i++) {

                                            if (eList.get(i) != null && !eList.get(i).equals("")) {

                                                eListNew.add(eList.get(i));
                                            }
                                        }
                                        tem = eListNew.toArray(new String[eListNew.size()]);
                                        for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                            rice.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                        }
                                        temp = br.readLine();
                                        //wheatA-wheatD-wheatB-rice-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatD-wheatB-rice-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    barley.append('-');
                                                    stiffbrome.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatD-wheatB-rice-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    barley.append('-');
                                                    stiffbrome.append('-');
                                                    sorghum.append('-');
                                                }
                                                continue;
                                            }

                                        } else {//wheatA-wheatD-wheatB-rice-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    barley.append('-');
                                                    stiffbrome.append('-');
                                                    maize.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatD-wheatB-rice
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    barley.append('-');
                                                    stiffbrome.append('-');
                                                    maize.append('-');
                                                    sorghum.append('-');
                                                }
                                                continue;
                                            }
                                        }
                                    } else { //wheatA-wheatD-wheatB-mazie-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {
                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatD-wheatB-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    barley.append('-');
                                                    stiffbrome.append('-');
                                                    rice.append('-');
                                                }

                                                continue;

                                            } else {
                                                //wheatA-wheatD-wheatB-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    barley.append('-');
                                                    stiffbrome.append('-');
                                                    rice.append('-');
                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }

                                        } else {//wheatA-wheatD-wheatB-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    barley.append('-');
                                                    stiffbrome.append('-');
                                                    rice.append('-');
                                                    maize.append('-');
                                                }

                                                continue;

                                            } else {//wheatA-wheatD-wheatB
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    barley.append('-');
                                                    stiffbrome.append('-');
                                                    rice.append('-');
                                                    maize.append('-');
                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        }
                                    }
                                }
                            }

                        } else {

                            //wheatA-wheatD-barley-
                            if (temp.startsWith("s barley")) {
                                List<String> cList = KStringUtils.fastSplit(temp);
                                List<String> cListNew = new ArrayList<>();
                                for (int i = 0; i < cList.size(); i++) {

                                    if (cList.get(i) != null && !cList.get(i).equals("")) {

                                        cListNew.add(cList.get(i));
                                    }
                                }
                                tem = cListNew.toArray(new String[cListNew.size()]);
                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                    barley.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                }
                                temp = br.readLine();
                                //wheatA-wheatD-barley-stiffbrome-
                                if (temp.startsWith("s stiffbrome")) {
                                    List<String> dList = KStringUtils.fastSplit(temp);
                                    List<String> dListNew = new ArrayList<>();
                                    for (int i = 0; i < dList.size(); i++) {

                                        if (dList.get(i) != null && !dList.get(i).equals("")) {

                                            dListNew.add(dList.get(i));
                                        }
                                    }
                                    tem = dListNew.toArray(new String[dListNew.size()]);
                                    for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                        stiffbrome.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                    }
                                    temp = br.readLine();
                                    //wheatA-wheatD-barley-stiffbrome-rice-
                                    if (temp.startsWith("s rice")) {
                                        List<String> eList = KStringUtils.fastSplit(temp);
                                        List<String> eListNew = new ArrayList<>();
                                        for (int i = 0; i < eList.size(); i++) {

                                            if (eList.get(i) != null && !eList.get(i).equals("")) {

                                                eListNew.add(eList.get(i));
                                            }
                                        }
                                        tem = eListNew.toArray(new String[eListNew.size()]);
                                        for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                            rice.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                        }
                                        temp = br.readLine();
                                        //wheatA-wheatD-barley-stiffbrome-rice-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatD-barley-stiffbrome-rice-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    wheatB.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatD-barley-stiffbrome-rice-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    wheatB.append('-');
                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        } else {//wheatA-wheatD-barley-stiffbrome-rice-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    wheatB.append('-');
                                                    maize.append('-');
                                                }

                                                continue;

                                            } else {//wheatA-wheatD-barley-stiffbrome-rice
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    wheatB.append('-');
                                                    maize.append('-');
                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        }

                                    } else {//wheatA-wheatD-barley-stiffbrome-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatD-barley-stiffbrome-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    wheatB.append('-');
                                                    rice.append('-');
                                                }

                                                continue;

                                            } else {    //wheatA-wheatD-barley-stiffbrome-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    wheatB.append('-');
                                                    rice.append('-');
                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        } else {//wheatA-wheatD-barley-stiffbrome-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    wheatB.append('-');
                                                    rice.append('-');
                                                    maize.append('-');
                                                }

                                                continue;

                                            } else {
                                                //wheatA-wheatD-barley-stiffbrome
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    wheatB.append('-');
                                                    rice.append('-');
                                                    maize.append('-');
                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        }
                                    }

                                } else { //wheatA-wheatD-barley-rice-
                                    if (temp.startsWith("s rice")) {
                                        List<String> eList = KStringUtils.fastSplit(temp);
                                        List<String> eListNew = new ArrayList<>();
                                        for (int i = 0; i < eList.size(); i++) {

                                            if (eList.get(i) != null && !eList.get(i).equals("")) {

                                                eListNew.add(eList.get(i));
                                            }
                                        }
                                        tem = eListNew.toArray(new String[eListNew.size()]);
                                        for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                            rice.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                        }
                                        temp = br.readLine();
                                        //wheatA-wheatD-barley-rice-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatD-barley-rice-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    wheatB.append('-');
                                                    stiffbrome.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatD-barley-rice-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    wheatB.append('-');
                                                    stiffbrome.append('-');
                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        } else {//wheatA-wheatD-barley-rice-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    wheatB.append('-');
                                                    stiffbrome.append('-');
                                                    maize.append('-');
                                                }

                                                continue;

                                            } else {//wheatA-wheatD-barley-rice
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    wheatB.append('-');
                                                    stiffbrome.append('-');
                                                    maize.append('-');
                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        }
                                    } else {//wheatA-wheatD-barley-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatD-barley-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    wheatB.append('-');
                                                    stiffbrome.append('-');
                                                    rice.append('-');
                                                }

                                                continue;

                                            } else {
                                                //wheatA-wheatD-barley-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    wheatB.append('-');
                                                    stiffbrome.append('-');
                                                    rice.append('-');
                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }

                                        } else {  //wheatA-wheatD-barley-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatB.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');

                                                    maize.append('-');
                                                }

                                                continue;

                                            } else {
                                                //wheatA-wheatD-barley
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatB.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');

                                                    maize.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        }
                                    }
                                }

                            } else {//wheatA-wheatD-stiffbrome-
                                if (temp.startsWith("s stiffbrome")) {
                                    List<String> dList = KStringUtils.fastSplit(temp);
                                    List<String> dListNew = new ArrayList<>();
                                    for (int i = 0; i < dList.size(); i++) {

                                        if (dList.get(i) != null && !dList.get(i).equals("")) {

                                            dListNew.add(dList.get(i));
                                        }
                                    }
                                    tem = dListNew.toArray(new String[dListNew.size()]);
                                    for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                        stiffbrome.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                    }
                                    temp = br.readLine();
                                    //wheatA-wheatD-stiffbrome-rice-
                                    if (temp.startsWith("s rice")) {
                                        List<String> eList = KStringUtils.fastSplit(temp);
                                        List<String> eListNew = new ArrayList<>();
                                        for (int i = 0; i < eList.size(); i++) {

                                            if (eList.get(i) != null && !eList.get(i).equals("")) {

                                                eListNew.add(eList.get(i));
                                            }
                                        }
                                        tem = eListNew.toArray(new String[eListNew.size()]);
                                        for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                            rice.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                        }
                                        temp = br.readLine();
                                        //wheatA-wheatD-stiffbrome-rice-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatD-stiffbrome-rice-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatB.append('-');

                                                    barley.append('-');
                                                }

                                                continue;

                                            } else {
                                                //wheatA-wheatD-stiffbrome-rice-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        } else {//wheatA-wheatD-stiffbrome-rice-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    maize.append('-');
                                                }

                                                continue;

                                            } else {
                                                //wheatA-wheatD-stiffbrome-rice
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    maize.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        }
                                    } else { //wheatA-wheatD-stiffbrome-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatD-stiffbrome-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    rice.append('-');
                                                }

                                                continue;

                                            } else {
                                                //wheatA-wheatD-stiffbrome-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    rice.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        } else {//wheatA-wheatD-stiffbrome-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    wheatB.append('-');
                                                    barley.append('-');
                                                    rice.append('-');
                                                    maize.append('-');
                                                }

                                                continue;

                                            } else {//wheatA-wheatD-stiffbrome
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    wheatB.append('-');
                                                    barley.append('-');
                                                    rice.append('-');
                                                    maize.append('-');
                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        }
                                    }

                                } else {//wheatA-wheatD-rice-
                                    if (temp.startsWith("s rice")) {
                                        List<String> eList = KStringUtils.fastSplit(temp);
                                        List<String> eListNew = new ArrayList<>();
                                        for (int i = 0; i < eList.size(); i++) {

                                            if (eList.get(i) != null && !eList.get(i).equals("")) {

                                                eListNew.add(eList.get(i));
                                            }
                                        }
                                        tem = eListNew.toArray(new String[eListNew.size()]);
                                        for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                            rice.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                        }
                                        temp = br.readLine();
                                        //wheatA-wheatD-rice-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatD-rice-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatD-rice-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatB.append('-');

                                                    stiffbrome.append('-');
                                                    barley.append('-');
                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }

                                        } else {//wheatA-wheatD-rice-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');

                                                    maize.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatD-rice
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');

                                                    maize.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        }
                                    } else {//wheatA-wheatD-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatD-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatD-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }

                                        } else { //wheatA-wheatD-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');

                                                    maize.append('-');
                                                }

                                                continue;

                                            } else {//wheatA-wheatD
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');

                                                    maize.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        }
                                    }
                                }
                            }
                        }

                    } else //wheatA-wheatB-
                    {
                        if (temp.startsWith("s wheatB")) {
                            List<String> aList = KStringUtils.fastSplit(temp);
                            List<String> aListNew = new ArrayList<>();
                            for (int i = 0; i < aList.size(); i++) {

                                if (aList.get(i) != null && !aList.get(i).equals("")) {

                                    aListNew.add(aList.get(i));
                                }
                            }
                            tem = aListNew.toArray(new String[aListNew.size()]);
                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                wheatB.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                            }
                            temp = br.readLine();
                            //wheatA-wheatB-barley-
                            if (temp.startsWith("s barley")) {
                                List<String> cList = KStringUtils.fastSplit(temp);
                                List<String> cListNew = new ArrayList<>();
                                for (int i = 0; i < cList.size(); i++) {

                                    if (cList.get(i) != null && !cList.get(i).equals("")) {

                                        cListNew.add(cList.get(i));
                                    }
                                }
                                tem = cListNew.toArray(new String[cListNew.size()]);
                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                    barley.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                }
                                temp = br.readLine();
                                //wheatA-wheatB-barley-stiffbrome-
                                if (temp.startsWith("s stiffbrome")) {
                                    List<String> dList = KStringUtils.fastSplit(temp);
                                    List<String> dListNew = new ArrayList<>();
                                    for (int i = 0; i < dList.size(); i++) {

                                        if (dList.get(i) != null && !dList.get(i).equals("")) {

                                            dListNew.add(dList.get(i));
                                        }
                                    }
                                    tem = dListNew.toArray(new String[dListNew.size()]);
                                    for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                        stiffbrome.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                    }
                                    temp = br.readLine();
                                    //wheatA-wheatB-barley-stiffbrome-rice-
                                    if (temp.startsWith("s rice")) {
                                        List<String> eList = KStringUtils.fastSplit(temp);
                                        List<String> eListNew = new ArrayList<>();
                                        for (int i = 0; i < eList.size(); i++) {

                                            if (eList.get(i) != null && !eList.get(i).equals("")) {

                                                eListNew.add(eList.get(i));
                                            }
                                        }
                                        tem = eListNew.toArray(new String[eListNew.size()]);
                                        for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                            rice.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                        }
                                        temp = br.readLine();
                                        //wheatA-wheatB-barley-stiffbrome-rice-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatB-barley-stiffbrome-rice-maize-sorghum

                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatB-barley-stiffbrome-rice-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }

                                        } else {
                                            //wheatA-wheatB-barley-stiffbrome-rice-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    maize.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatB-barley-stiffbrome-rice
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    maize.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        }

                                    } else {////wheatA-wheatB-barley-stiffbrome-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatB-barley-stiffbrome-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    rice.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatB-barley-stiffbrome-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    rice.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }

                                        } else {////wheatA-wheatB-barley-stiffbrome-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    rice.append('-');

                                                    maize.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatB-barley-stiffbrome
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    rice.append('-');

                                                    maize.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        }
                                    }

                                } else {
                                    // //wheatA-wheatB-barley-rice-
                                    if (temp.startsWith("s rice")) {
                                        List<String> eList = KStringUtils.fastSplit(temp);
                                        List<String> eListNew = new ArrayList<>();
                                        for (int i = 0; i < eList.size(); i++) {

                                            if (eList.get(i) != null && !eList.get(i).equals("")) {

                                                eListNew.add(eList.get(i));
                                            }
                                        }
                                        tem = eListNew.toArray(new String[eListNew.size()]);
                                        for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                            rice.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                        }
                                        temp = br.readLine();
                                        //wheatA-wheatB-barley-rice-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatB-barley-rice-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    stiffbrome.append('-');
                                                }
                                                continue;
                                            } else {
                                                //wheatA-wheatB-barley-rice-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    stiffbrome.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }

                                        } else {
                                            //wheatA-wheatB-barley-rice-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    stiffbrome.append('-');

                                                    maize.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatB-barley-rice
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    stiffbrome.append('-');

                                                    maize.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        }

                                    } else { //wheatA-wheatB-barley-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatB-barley-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatB-barley-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }

                                        } else {// //wheatA-wheatB-barley-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');

                                                    maize.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatB-barley
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');

                                                    maize.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        }
                                    }
                                }

                            } else {////wheatA-wheatB-stiffbrome
                                if (temp.startsWith("s stiffbrome")) {
                                    List<String> dList = KStringUtils.fastSplit(temp);
                                    List<String> dListNew = new ArrayList<>();
                                    for (int i = 0; i < dList.size(); i++) {

                                        if (dList.get(i) != null && !dList.get(i).equals("")) {

                                            dListNew.add(dList.get(i));
                                        }
                                    }
                                    tem = dListNew.toArray(new String[dListNew.size()]);
                                    for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                        stiffbrome.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                    }
                                    temp = br.readLine();
                                    //wheatA-wheatB-stiffbrome-rice-
                                    if (temp.startsWith("s rice")) {
                                        List<String> eList = KStringUtils.fastSplit(temp);
                                        List<String> eListNew = new ArrayList<>();
                                        for (int i = 0; i < eList.size(); i++) {

                                            if (eList.get(i) != null && !eList.get(i).equals("")) {

                                                eListNew.add(eList.get(i));
                                            }
                                        }
                                        tem = eListNew.toArray(new String[eListNew.size()]);
                                        for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                            rice.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                        }
                                        temp = br.readLine();
                                        //wheatA-wheatB-stiffbrome-rice-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatB-stiffbrome-rice-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    wheatD.append('-');
                                                    barley.append('-');
                                                    
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatB-stiffbrome-rice-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    barley.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        } else {
                                            ////wheatA-wheatB-stiffbrome-rice-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    barley.append('-');

                                                    maize.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatB-stiffbrome-rice
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    barley.append('-');

                                                    maize.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        }
                                    } else {
                                        // //wheatA-wheatB-stiffbrome-mazie-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatB-stiffbrome-mazie-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    barley.append('-');

                                                    rice.append('-');
                                                }

                                                continue;
                                            } else {//wheatA-wheatB-stiffbrome-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    barley.append('-');

                                                    rice.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        } else {////wheatA-wheatB-stiffbrome-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');
                                                    barley.append('-');
                                                    rice.append('-');

                                                    maize.append('-');
                                                }

                                            } else {
                                                //wheatA-wheatB-stiffbrome 
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    barley.append('-');

                                                    rice.append('-');

                                                    maize.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        }
                                    }
                                } else {// //wheatA-wheatB-rice-
                                    if (temp.startsWith("s rice")) {
                                        List<String> eList = KStringUtils.fastSplit(temp);
                                        List<String> eListNew = new ArrayList<>();
                                        for (int i = 0; i < eList.size(); i++) {

                                            if (eList.get(i) != null && !eList.get(i).equals("")) {

                                                eListNew.add(eList.get(i));
                                            }
                                        }
                                        tem = eListNew.toArray(new String[eListNew.size()]);
                                        for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                            rice.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                        }
                                        temp = br.readLine();
                                        ////wheatA-wheatB-rice-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            ////wheatA-wheatB-rice-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');
                                                }

                                                continue;
                                            } else {
                                                ////wheatA-wheatB-rice-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }

                                        } else {//wheatA-wheatB-rice-sorghum

                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');

                                                    maize.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatB-rice
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');

                                                    maize.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        }

                                    } else {//wheatA-wheatB-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-wheatB-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatB-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        } else {  //wheatA-wheatB-sorghum

                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');
                                                         barley.append('-');
                                                    stiffbrome.append('-');

                                                    rice.append('-');

                                                    maize.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-wheatB
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');

                                                    sorghum.append('-');
                                                    maize.append('-');
                                                }

                                                continue;
                                            }
                                        }
                                    }
                                }
                            }

                        } else {
                            //wheatA-barley
                            if (temp.startsWith("s barley")) {
                                List<String> cList = KStringUtils.fastSplit(temp);
                                List<String> cListNew = new ArrayList<>();
                                for (int i = 0; i < cList.size(); i++) {

                                    if (cList.get(i) != null && !cList.get(i).equals("")) {

                                        cListNew.add(cList.get(i));
                                    }
                                }
                                tem = cListNew.toArray(new String[cListNew.size()]);
                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                    barley.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                }
                                temp = br.readLine();
                                //wheatA-barley-stiffbrome-
                                if (temp.startsWith("s stiffbrome")) {
                                    List<String> dList = KStringUtils.fastSplit(temp);
                                    List<String> dListNew = new ArrayList<>();
                                    for (int i = 0; i < dList.size(); i++) {

                                        if (dList.get(i) != null && !dList.get(i).equals("")) {

                                            dListNew.add(dList.get(i));
                                        }
                                    }
                                    tem = dListNew.toArray(new String[dListNew.size()]);
                                    for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                        stiffbrome.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                    }
                                    temp = br.readLine();
                                    //wheatA-barley-stiffbrome-rice-
                                    if (temp.startsWith("s rice")) {
                                        List<String> eList = KStringUtils.fastSplit(temp);
                                        List<String> eListNew = new ArrayList<>();
                                        for (int i = 0; i < eList.size(); i++) {

                                            if (eList.get(i) != null && !eList.get(i).equals("")) {

                                                eListNew.add(eList.get(i));
                                            }
                                        }
                                        tem = eListNew.toArray(new String[eListNew.size()]);
                                        for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                            rice.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                        }
                                        temp = br.readLine();
                                        //wheatA-barley-stiffbrome-rice-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-barley-stiffbrome-rice-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    wheatB.append('-');
                                                }

                                                continue;

                                            } else {
                                                //wheatA-barley-stiffbrome-rice-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    wheatD.append('-');
                                                    wheatB.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        } else {
                                            //wheatA-barley-stiffbrome-rice-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    maize.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-barley-stiffbrome-rice
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    maize.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        }

                                    } else {
                                        //wheatA-barley-stiffbrome-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-barley-stiffbrome-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    rice.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-barley-stiffbrome-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    rice.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }

                                        } else {
                                            // //wheatA-barley-stiffbrome-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    rice.append('-');

                                                    maize.append('-');
                                                }

                                                continue;
                                            } else {

                                                //wheatA-barley-stiffbrome
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    rice.append('-');

                                                    maize.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        }
                                    }
                                } else {
                                    //wheatA-barley-rice-
                                    if (temp.startsWith("s rice")) {
                                        List<String> eList = KStringUtils.fastSplit(temp);
                                        List<String> eListNew = new ArrayList<>();
                                        for (int i = 0; i < eList.size(); i++) {

                                            if (eList.get(i) != null && !eList.get(i).equals("")) {

                                                eListNew.add(eList.get(i));
                                            }
                                        }
                                        tem = eListNew.toArray(new String[eListNew.size()]);
                                        for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                            rice.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                        }
                                        temp = br.readLine();
                                        // //wheatA-barley-rice-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            // //wheatA-barley-rice-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    stiffbrome.append('-');
                                                }

                                                continue;
                                            } else {

                                                //wheatA-barley-rice-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    stiffbrome.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        } else {
                                            ////wheatA-barley-rice-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    stiffbrome.append('-');

                                                    maize.append('-');
                                                }

                                                continue;
                                            } else {      ////wheatA-barley-rice
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    maize.append('-');
                                                    stiffbrome.append('-');
                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        }
                                    } else {
                                        //wheatA-barley-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-barley-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-barley-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        } else {
                                            //wheatA-barley-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');

                                                    maize.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-barley
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');

                                                    maize.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        }
                                    }
                                }

                            } else {
                                //wheatA-stiffbrome-
                                if (temp.startsWith("s stiffbrome")) {
                                    List<String> dList = KStringUtils.fastSplit(temp);
                                    List<String> dListNew = new ArrayList<>();
                                    for (int i = 0; i < dList.size(); i++) {

                                        if (dList.get(i) != null && !dList.get(i).equals("")) {

                                            dListNew.add(dList.get(i));
                                        }
                                    }
                                    tem = dListNew.toArray(new String[dListNew.size()]);
                                    for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                        stiffbrome.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                    }
                                    temp = br.readLine();
                                    //wheatA-stiffbrome-rice-
                                    if (temp.startsWith("s rice")) {
                                        List<String> eList = KStringUtils.fastSplit(temp);
                                        List<String> eListNew = new ArrayList<>();
                                        for (int i = 0; i < eList.size(); i++) {

                                            if (eList.get(i) != null && !eList.get(i).equals("")) {

                                                eListNew.add(eList.get(i));
                                            }
                                        }
                                        tem = eListNew.toArray(new String[eListNew.size()]);
                                        for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                            rice.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                        }
                                        temp = br.readLine();
                                        //wheatA-stiffbrome-rice-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-stiffbrome-rice-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    barley.append('-');
                                                }

                                                continue;
                                            } else {

                                                //wheatA-stiffbrome-rice-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        } else {
                                            //wheatA-stiffbrome-rice-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    maize.append('-');
                                                }

                                                continue;
                                            } else {
                                                ////wheatA-stiffbrome-rice
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    maize.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        }

                                    } else {
                                        //wheatA-stiffbrome-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }

                                            temp = br.readLine();
                                            //wheatA-stiffbrome-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    rice.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-stiffbrome-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    wheatD.append('-');
                                                    wheatB.append('-');

                                                    barley.append('-');
                                                     rice.append('-');
                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        } else {
                                            ////wheatA-stiffbrome-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    rice.append('-');

                                                    maize.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-stiffbrome
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    rice.append('-');

                                                    maize.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;

                                            }
                                        }
                                    }

                                } else {
                                    //wheatA-rice-
                                    if (temp.startsWith("s rice")) {
                                        List<String> eList = KStringUtils.fastSplit(temp);
                                        List<String> eListNew = new ArrayList<>();
                                        for (int i = 0; i < eList.size(); i++) {

                                            if (eList.get(i) != null && !eList.get(i).equals("")) {

                                                eListNew.add(eList.get(i));
                                            }
                                        }
                                        tem = eListNew.toArray(new String[eListNew.size()]);
                                        for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                            rice.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                        }
                                        temp = br.readLine();
                                        //wheatA-rice-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //wheatA-rice-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-rice-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }

                                        } else {
                                            //wheatA-rice-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);

                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');

                                                    maize.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-rice
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');

                                                    maize.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        }

                                    } else {
                                        //whestA-maize-
                                        if (temp.startsWith("s maize")) {
                                            List<String> fList = KStringUtils.fastSplit(temp);
                                            List<String> fListNew = new ArrayList<>();
                                            for (int i = 0; i < fList.size(); i++) {

                                                if (fList.get(i) != null && !fList.get(i).equals("")) {

                                                    fListNew.add(fList.get(i));
                                                }
                                            }
                                            tem = fListNew.toArray(new String[fListNew.size()]);
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                maize.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                            }
                                            temp = br.readLine();
                                            //whestA-maize-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));
                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');
                                                }

                                                continue;
                                            } else {
                                                //wheatA-maize
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            }
                                        } else {
                                            //wheatA-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);
                                                for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {
                                                    sorghum.append(tem[6].charAt(finalAlignedFourDMap.get(alignedNum).get(j)));

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');

                                                    maize.append('-');
                                                }

                                                continue;
                                            }else{
                                            for (int j = 0; j < finalAlignedFourDMap.get(alignedNum).size(); j++) {

                                                    wheatD.append('-');

                                                    wheatB.append('-');

                                                    barley.append('-');

                                                    stiffbrome.append('-');

                                                    rice.append('-');
                                                    maize.append('-');

                                                    sorghum.append('-');
                                                }

                                                continue;
                                            
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                }
            }
            System.out.println("all cut position established, starting random cutting "+phylipNum +" positions");
            int wheatAFourDLength=wheatA.toString().length();
            System.out.println(wheatAFourDLength);
            System.out.println(wheatD.toString().length());
            System.out.println(wheatB.toString().length());
            System.out.println(barley.toString().length());
            System.out.println(stiffbrome.toString().length());
            System.out.println(rice.toString().length());
            System.out.println(maize.toString().length());
            System.out.println(sorghum.toString().length());
            int[] randomPos = randomCommon(0, wheatAFourDLength-1, phylipNum);
             System.out.println("randomPos established"+" "+randomPos.length);
             
            for (int i=0;i<randomPos.length;i++) {
        
                // if (a != '-' && b != '-' && c != '-' && d != '-' && e != '-' && f != '-' && g != '-' && h != '-') {
                randomWheatAFourD.append(wheatA.charAt(randomPos[i]));
                randomWheatDFourD.append(wheatD.charAt(randomPos[i]));
                randomWheatBFourD.append(wheatB.charAt(randomPos[i]));
                randomBarleyFourD.append(barley.charAt(randomPos[i]));
                randomStiffbromeFourD.append(stiffbrome.charAt(randomPos[i]));
                randomRiceFourD.append(rice.charAt(randomPos[i]));
                randomMaizeFourD.append(maize.charAt(randomPos[i]));
                randomSorghumFourD.append(sorghum.charAt(randomPos[i]));
                // }
            }
            System.out.println("random Array build successfully");
            bw.write("8 ");
            bw.write(String.valueOf(phylipNum));
            bw.newLine();
            bw.write("wheatA        ");
 System.out.println("wheatA title build successfully");
            bw.write(randomWheatAFourD.toString());
            bw.flush();
    System.out.println("randomwheatA write successfully");        
            bw.newLine();
            bw.write("wheatD        ");

            bw.write(randomWheatDFourD.toString());
              bw.flush();
            bw.newLine();
            bw.write("wheatB        ");

            bw.write(randomWheatBFourD.toString());
              bw.flush();
            bw.newLine();
            bw.write("barley        ");

            bw.write(randomBarleyFourD.toString());
              bw.flush();
            bw.newLine();
            bw.write("stiffbrome    ");

            bw.write(randomStiffbromeFourD.toString());
              bw.flush();
            bw.newLine();
            bw.write("rice          ");

            bw.write(randomRiceFourD.toString());
              bw.flush();
            bw.newLine();
            bw.write("maize         ");

            bw.write(randomMaizeFourD.toString());
              bw.flush();
            bw.newLine();
            bw.write("sorghum       ");

            bw.write(randomSorghumFourD.toString());
       
            bw.flush();
            bw.close();

        } catch (Exception e) {
            System.out.println("Error in Reading and Write MSA");
             e.printStackTrace(); 
        }
    }
}
