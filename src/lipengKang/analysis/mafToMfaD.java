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
import  lipengKang.analysis.KStringUtils;
import zhouyao.analysis.wheatHapMap.YaoIOUtils;

/**
 *
 * @author kanglipeng
 */
//put all blocks from maf to one fasta
public class mafToMfaD {

    HashMap<Integer, TIntArrayList> alignedPosMap = new HashMap<>();
    ArrayList<Integer> alignedSegment = new ArrayList();
    TIntArrayList alignedPos = new TIntArrayList();

    public mafToMfaD(String mafDir, String wheatADir, String wheatDDir, String wheatBDir, String riceDir, String sorghumDir, String maizeDir, String barleyDir, String stiffbromeDir) {
        parseMaf(mafDir, wheatADir, wheatDDir, wheatBDir, riceDir, sorghumDir, maizeDir, barleyDir, stiffbromeDir);
    }

    public void parseMaf(String mafDir, String wheatADir, String wheatDDir, String wheatBDir, String riceDir, String sorghumDir, String maizeDir, String barleyDir, String stiffbromeDir) {
        BufferedWriter bb;
        BufferedWriter ba;
        BufferedWriter bc;
        BufferedWriter bd;
        BufferedWriter be;
        BufferedWriter bf;
        BufferedWriter bg;
        BufferedWriter bh;

        BufferedReader br;

        br = YaoIOUtils.getTextReader(mafDir);
        ba = YaoIOUtils.getTextWriter(wheatADir);
        bb = YaoIOUtils.getTextWriter(wheatDDir);
        bc = YaoIOUtils.getTextWriter(wheatBDir);
        bd = YaoIOUtils.getTextWriter(barleyDir);
        be = YaoIOUtils.getTextWriter(stiffbromeDir);
        bf = YaoIOUtils.getTextWriter(riceDir);
        bg = YaoIOUtils.getTextWriter(maizeDir);
        bh = YaoIOUtils.getTextWriter(sorghumDir);

        String temp = null;
        String[] tem = null;

        try {
            ba.write(">wheatA");
            ba.newLine();
            bb.write(">wheatD");
            bb.newLine();
            bc.write(">wheatB");
            bc.newLine();
            bd.write(">barley");
            bd.newLine();
            be.write(">stiffbrome");
            be.newLine();
            bf.write(">rice");
            bf.newLine();
            bg.write(">maize");
            bg.newLine();
            bh.write(">sorghum");
            bh.newLine();
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("a")) {
                    StringBuilder wheatD = new StringBuilder();

                    StringBuilder wheatA = new StringBuilder();

                    StringBuilder wheatB = new StringBuilder();

                    StringBuilder barley = new StringBuilder();

                    StringBuilder stiffbrome = new StringBuilder();

                    StringBuilder rice = new StringBuilder();

                    StringBuilder maize = new StringBuilder();

                    StringBuilder sorghum = new StringBuilder();
                    temp = br.readLine();
                    if (temp.startsWith("s wheatD")) {
                        List<String> aList = KStringUtils.fastSplit(temp);
                        List<String> aListNew = new ArrayList<>();
                        for (int i = 0; i < aList.size(); i++) {

                            if (aList.get(i) != null && !aList.get(i).equals("")) {

                                aListNew.add(aList.get(i));
                            }
                        }
                        tem = aListNew.toArray(new String[aListNew.size()]);
                        wheatD.append(tem[6]);

                    }
                    temp = br.readLine();

                    //wheatA-wheatD-
                    if (temp.startsWith("s wheatA")) {
                        List<String> bList = KStringUtils.fastSplit(temp);
                        List<String> bListNew = new ArrayList<>();
                        for (int i = 0; i < bList.size(); i++) {

                            if (bList.get(i) != null && !bList.get(i).equals("")) {

                                bListNew.add(bList.get(i));
                            }
                        }
                        tem = bListNew.toArray(new String[bListNew.size()]);
                        wheatA.append(tem[6]);
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
                            wheatB.append(tem[6]);
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
                                barley.append(tem[6]);
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
                                    stiffbrome.append(tem[6]);
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
                                        rice.append(tem[6]);
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
                                            maize.append(tem[6]);
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
                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-wheatD-wheatB-barley-stiffbrome-rice-maize
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }

                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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
                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());

                                                continue;
                                            } else {
                                                //wheatA-wheatD-wheatB-barley-stiffbrome-rice
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());

                                                continue;
                                            } else {
                                                //wheatA-wheatD-wheatB-barley-stiffbrome-maize
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());

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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bf.write(rice.toString());

                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {
                                                //wheatA-wheatD-wheatB-barley-stiffbrome
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                        rice.append(tem[6]);
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {    //wheatA-wheatD-wheatB-barley-rice-maize
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {//wheatA-wheatD-wheatB-barley-rice
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                            maize.append(tem[6]);
                                            temp = br.readLine();
                                            //wheatA-wheatD-wheatB-barley-mazie-sorghum
                                            if (temp.startsWith("s sorghum")) {
                                                List<String> gList = KStringUtils.fastSplit(temp);
                                                List<String> gListNew = new ArrayList<>();
                                                for (int i = 0; i < gList.size(); i++) {

                                                    if (gList.get(i) != null && !gList.get(i).equals("")) {

                                                        gListNew.add(gList.get(i));
                                                    }
                                                }
                                                tem = gListNew.toArray(new String[gListNew.size()]);

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());

                                                continue;

                                            } else {
                                                //wheatA-wheatD-wheatB-barley-mazie
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());

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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {
                                                //wheatA-wheatD-wheatB-barley
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                    stiffbrome.append(tem[6]);
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

                                        rice.append(tem[6]);
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
                                            maize.append(tem[6]);
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
                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {
                                                //wheatA-wheatD-wheatB-stiffbrome-rice-maize
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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
                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {
                                                //wheatA-wheatD-wheatB-stiffbrome-rice
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());

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

                                            maize.append(tem[6]);
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
                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {
                                                //wheatA-wheatD-wheatB-stiffbrome-maize
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());

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
                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {//wheatA-wheatD-wheatB-stiffbrome
                                               //!!!!!
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());

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

                                        rice.append(tem[6]);
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());

                                                continue;
                                            } else {
                                                //wheatA-wheatD-wheatB-rice-maize
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }

                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-wheatD-wheatB-rice
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {
                                                //wheatA-wheatD-wheatB-maize
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());

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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {//wheatA-wheatD-wheatB
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());

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

                                barley.append(tem[6]);
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

                                    stiffbrome.append(tem[6]);
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

                                        rice.append(tem[6]);
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());

                                                continue;
                                            } else {
                                                //wheatA-wheatD-barley-stiffbrome-rice-maize
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());

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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());

                                                continue;

                                            } else {//wheatA-wheatD-barley-stiffbrome-rice
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());

                                                continue;

                                            } else {    //wheatA-wheatD-barley-stiffbrome-maize
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {
                                                //wheatA-wheatD-barley-stiffbrome
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                        rice.append(tem[6]);
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-wheatD-barley-rice-maize
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {//wheatA-wheatD-barley-rice
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {
                                                //wheatA-wheatD-barley-maize
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {
                                                //wheatA-wheatD-barley
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                    stiffbrome.append(tem[6]);
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

                                        rice.append(tem[6]);
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {
                                                //wheatA-wheatD-stiffbrome-rice-maize
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());

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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {
                                                //wheatA-wheatD-stiffbrome-rice
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {
                                                //wheatA-wheatD-stiffbrome-maize
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {//wheatA-wheatD-stiffbrome
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());

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

                                        rice.append(tem[6]);
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());

                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-wheatD-rice-maize
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());

                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }

                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-wheatD-rice
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());

                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-wheatD-maize
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());

                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }

                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {//wheatA-wheatD
                                                bb.write(wheatD.toString());
                                                ba.write(wheatA.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                            wheatB.append(tem[6]);
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

                                barley.append(tem[6]);
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

                                    stiffbrome.append(tem[6]);
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

                                        rice.append(tem[6]);
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());

                                                continue;
                                            } else {
                                                //wheatA-wheatB-barley-stiffbrome-rice-maize
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-wheatB-barley-stiffbrome-rice
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-wheatB-barley-stiffbrome-maize
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-wheatB-barley-stiffbrome
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                        rice.append(tem[6]);
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-wheatB-barley-rice-maize
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-wheatB-barley-rice
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-wheatB-barley-maize
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-wheatB-barley
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                    stiffbrome.append(tem[6]);
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

                                        rice.append(tem[6]);
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-wheatB-stiffbrome-rice-maize
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-wheatB-stiffbrome-rice
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {//wheatA-wheatB-stiffbrome-maize
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-wheatB-stiffbrome 
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                        rice.append(tem[6]);
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                ////wheatA-wheatB-rice-maize
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-wheatB-rice
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-wheatB-maize
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-wheatB
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                barley.append(tem[6]);
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

                                    stiffbrome.append(tem[6]);
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

                                        rice.append(tem[6]);
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;

                                            } else {
                                                //wheatA-barley-stiffbrome-rice-maize
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-barley-stiffbrome-rice
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-barley-stiffbrome-maize
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {

                                                //wheatA-barley-stiffbrome
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                               
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                        rice.append(tem[6]);
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {

                                                //wheatA-barley-rice-maize
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());

                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {      ////wheatA-barley-rice
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());

                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-barley-maize
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-barley
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                    stiffbrome.append(tem[6]);
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

                                        rice.append(tem[6]);
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {

                                                //wheatA-stiffbrome-rice-maize
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                ////wheatA-stiffbrome-rice
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-stiffbrome-maize
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            } else {
                                                //wheatA-stiffbrome
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                        rice.append(tem[6]);
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

                                            maize.append(tem[6]);
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());

                                                continue;
                                            } else {
                                                //wheatA-rice-maize
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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

                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());

                                                continue;
                                            } else {
                                                //wheatA-rice
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
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

                                            maize.append(tem[6]);
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
                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());

                                                continue;
                                            } else {
                                                //wheatA-maize
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                bh.write(sorghum.toString());
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
                                                sorghum.append(tem[6]);
                                                bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                                continue;
                                            }else{
                                             bb.write(wheatD.toString());
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatA.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    wheatB.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    barley.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    stiffbrome.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    rice.append('-');
                                                }
                                                for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    maize.append('-');
                                                }
                                                 for (int i = 0; i < wheatD.toString().length(); i++) {
                                                    sorghum.append('-');
                                                }
                                                ba.write(wheatA.toString());
                                                bc.write(wheatB.toString());
                                                bd.write(barley.toString());
                                                be.write(stiffbrome.toString());
                                                bf.write(rice.toString());
                                                bg.write(maize.toString());
                                                bh.write(sorghum.toString());
                                            
                                            
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                }
            }

            bb.flush();
            bb.close();
            ba.flush();
            ba.close();
            bc.flush();
            bc.close();
            bd.flush();
            bd.close();
            be.flush();
            be.close();
            bf.flush();
            bf.close();
            bg.flush();
            bg.close();
            bh.flush();
            bh.close();

        } catch (Exception e) {
            System.out.println("Error in Reading and Write MSA");
            e.printStackTrace();
        }

    }

    public static void main(String[] args) {
    String wheatADir =    "/data1/home/lipeng/result/CNEr/axt/D/gerp/wheatA.mfa";
        String wheatDDir = "/data1/home/lipeng/result/CNEr/axt/D/gerp/wheatD.mfa";
        String wheatBDir = "/data1/home/lipeng/result/CNEr/axt/D/gerp/wheatB.mfa";
        String barleyDir = "/data1/home/lipeng/result/CNEr/axt/D/gerp/barley.mfa";
        String riceDir =   "/data1/home/lipeng/result/CNEr/axt/D/gerp/rice.mfa";
        String stiffbromeDir = "/data1/home/lipeng/result/CNEr/axt/D/gerp/stiffbrome.mfa";
        String maizeDir =  "/data1/home/lipeng/result/CNEr/axt/D/gerp/maize.mfa";
        String sorghumDir ="/data1/home/lipeng/result/CNEr/axt/D/gerp/sorghum.mfa";
        String mafDir =    "/data1/home/lipeng/result/CNEr/axt/D/gerp/wheatD.maf";

       /*String wheatADir = "/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/wheatA.mfa";
        String wheatDDir = "/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/wheatD.mfa";
        String wheatBDir = "/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/wheatB.mfa";
        String barleyDir = "/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/barley.mfa";
        String riceDir = "/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/rice.mfa";
        String stiffbromeDir = "/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/stiffbrome.mfa";
        String maizeDir = "/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/maize.mfa";
        String sorghumDir = "/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/sorghum.mfa";
        String mafDir = "/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/test.maf";*/
        new mafToMfaD(mafDir, wheatADir, wheatDDir, wheatBDir, riceDir, sorghumDir, maizeDir, barleyDir, stiffbromeDir);
    }
}
