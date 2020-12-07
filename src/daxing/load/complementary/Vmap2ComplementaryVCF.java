package daxing.load.complementary;

import daxing.common.IOTool;
import daxing.common.RowTableTool;
import daxing.common.StringTool;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.stat.inference.TestUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

public class Vmap2ComplementaryVCF {

    /**
     * 每个block计算一个加性和显性下的t统计量
     */
    public static void calculateT(String inputFile, String triadsBlockInfoFile, String outFile){
        Map<String, String> taxonGroupMap= RowTableTool.getMap(triadsBlockInfoFile, 0, 2);
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            bw.write("TriadsBlockID\tGroup\tGenotypedHexaploidTaxaNum\tGenotypedPseudoTaxaNum" +
                    "\tSlightlyOrStrongly\tAdditiveOrDominance" +
                    "\tT_equalVariances" +
                    "\tp_equalVariances\tT_notEqualVariances\tp_notEqualVariances");
            bw.newLine();
            String line;
            List<String> headerList= PStringUtils.fastSplit(br.readLine(), "\t");
            TIntArrayList pseudoIndexList=new TIntArrayList();
            TIntArrayList landraceIndexList=new TIntArrayList();
            TIntArrayList cultivarIndexList=new TIntArrayList();
            for (int i = 3; i < headerList.size(); i++) {
                if (taxonGroupMap.get(headerList.get(i)).length()>2){
                    pseudoIndexList.add(i-3);
                }else if (taxonGroupMap.get(headerList.get(i)).equals("CL")){
                    cultivarIndexList.add(i-3);
                }else if (taxonGroupMap.get(headerList.get(i)).equals("LR")){
                    landraceIndexList.add(i-3);
                }
            }
            List<double[][]> loadList, loadLRList, loadCLList, loadPseudoList;
            /**
             * dim 1 2:slightly/strongly, additive/dominance
             */
            TDoubleArrayList[][] loadLR;
            TDoubleArrayList[][] loadCL;
            TDoubleArrayList[][] loadPseudo;


            String triadsBlockID, group, slightly_strongly, additive_dominance;
            double t_equalVariants, t_notEqualVariants, p_equalVariants, p_notEqualVariants;
            int taxonNum_hexaploid, taxonNum_pseudo;
            StringBuilder sb=new StringBuilder();

            while ((line=br.readLine())!=null){
                loadList=calculateLoad(line);
                triadsBlockID=PStringUtils.fastSplit(line).get(0);


                loadLRList=sublist(loadList, landraceIndexList);
                loadCLList=sublist(loadList, cultivarIndexList);
                loadPseudoList=sublist(loadList, pseudoIndexList);

                loadLR=new TDoubleArrayList[2][];
                loadCL=new TDoubleArrayList[2][];
                loadPseudo=new TDoubleArrayList[2][];

                for (int i = 0; i < loadLR.length; i++) {
                    loadLR[i]=new TDoubleArrayList[2];
                    for (int j = 0; j < loadLR[i].length; j++) {
                        loadLR[i][j]=new TDoubleArrayList();
                    }
                }

                for (int i = 0; i < loadCL.length; i++) {
                    loadCL[i]=new TDoubleArrayList[2];
                    for (int j = 0; j < loadCL[i].length; j++) {
                        loadCL[i][j]=new TDoubleArrayList();
                    }
                }

                for (int i = 0; i < loadPseudo.length; i++) {
                    loadPseudo[i]=new TDoubleArrayList[2];
                    for (int j = 0; j < loadPseudo[i].length; j++) {
                        loadPseudo[i][j]=new TDoubleArrayList();
                    }
                }
                for (int i = 0; i < loadLRList.size(); i++) {
                    for (int j = 0; j < loadLRList.get(i).length; j++) {
                        for (int k = 0; k < loadLRList.get(i)[j].length; k++) {
                            if (loadLRList.get(i)[j][k] < 0) continue;
                            loadLR[j][k].add(loadLRList.get(i)[j][k]);
                        }
                    }
                }

                for (int i = 0; i < loadCLList.size(); i++) {
                    for (int j = 0; j < loadCLList.get(i).length; j++) {
                        for (int k = 0; k < loadCLList.get(i)[j].length; k++) {
                            if (loadCLList.get(i)[j][k] < 0) continue;
                            loadCL[j][k].add(loadCLList.get(i)[j][k]);
                        }
                    }
                }

                for (int i = 0; i < loadPseudoList.size(); i++) {
                    for (int j = 0; j < loadPseudoList.get(i).length; j++) {
                        for (int k = 0; k < loadPseudoList.get(i)[j].length; k++) {
                            if (loadPseudoList.get(i)[j][k] < 0) continue;
                            loadPseudo[j][k].add(loadPseudoList.get(i)[j][k]);
                        }
                    }
                }

                for (int i = 0; i < loadPseudo.length; i++) {
                    group="LR";
                    slightly_strongly= i==0 ? "slightly" : "strongly";
                    for (int j = 0; j < loadLR[i].length; j++) {
                        sb.setLength(0);
                        taxonNum_hexaploid=loadLR[i][j].size();
                        taxonNum_pseudo=loadPseudo[i][j].size();
                        additive_dominance= j==0 ? "additive" : "dominance";

                        if (taxonNum_hexaploid < 5 || taxonNum_pseudo<5){
                            sb.append(triadsBlockID).append("\t").append(group).append("\t");
                            sb.append(taxonNum_hexaploid).append("\t").append(taxonNum_pseudo).append("\t");
                            sb.append(slightly_strongly).append("\t").append(additive_dominance).append("\t");
                            sb.append("NA").append("\t").append("NA").append("\t");
                            sb.append("NA").append("\t").append("NA");
                            bw.write(sb.toString());
                            bw.newLine();
                            continue;
                        }

                        t_equalVariants=TestUtils.homoscedasticT(loadLR[i][j].toArray(), loadPseudo[i][j].toArray());
                        p_equalVariants=TestUtils.homoscedasticTTest(loadLR[i][j].toArray(), loadPseudo[i][j].toArray());
                        t_notEqualVariants=TestUtils.t(loadLR[i][j].toArray(), loadPseudo[i][j].toArray());
                        p_notEqualVariants=TestUtils.tTest(loadLR[i][j].toArray(), loadPseudo[i][j].toArray());


                        sb.append(triadsBlockID).append("\t").append(group).append("\t");
                        sb.append(taxonNum_hexaploid).append("\t").append(taxonNum_pseudo).append("\t");
                        sb.append(slightly_strongly).append("\t").append(additive_dominance).append("\t");
                        sb.append(t_equalVariants).append("\t").append(p_equalVariants).append("\t");
                        sb.append(t_notEqualVariants).append("\t").append(p_notEqualVariants);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }

                for (int i = 0; i < loadPseudo.length; i++) {
                    group="CL";
                    slightly_strongly= i==0 ? "slightly" : "strongly";
                    for (int j = 0; j < loadCL[i].length; j++) {
                        sb.setLength(0);
                        taxonNum_hexaploid=loadCL[i][j].size();
                        taxonNum_pseudo=loadPseudo[i][j].size();
                        additive_dominance= j==0 ? "additive" : "dominance";

                        if (taxonNum_hexaploid < 5 || taxonNum_pseudo<5){
                            sb.append(triadsBlockID).append("\t").append(group).append("\t");
                            sb.append(taxonNum_hexaploid).append("\t").append(taxonNum_pseudo).append("\t");
                            sb.append(slightly_strongly).append("\t").append(additive_dominance).append("\t");
                            sb.append("NA").append("\t").append("NA").append("\t");
                            sb.append("NA").append("\t").append("NA");
                            bw.write(sb.toString());
                            bw.newLine();
                            continue;
                        }

                        t_equalVariants=TestUtils.homoscedasticT(loadCL[i][j].toArray(), loadPseudo[i][j].toArray());
                        t_notEqualVariants=TestUtils.t(loadCL[i][j].toArray(), loadPseudo[i][j].toArray());
                        p_equalVariants=TestUtils.homoscedasticTTest(loadCL[i][j].toArray(), loadPseudo[i][j].toArray());
                        p_notEqualVariants=TestUtils.tTest(loadCL[i][j].toArray(), loadPseudo[i][j].toArray());


                        sb.append(triadsBlockID).append("\t").append(group).append("\t");
                        sb.append(taxonNum_hexaploid).append("\t").append(taxonNum_pseudo).append("\t");
                        sb.append(slightly_strongly).append("\t").append(additive_dominance).append("\t");
                        sb.append(t_equalVariants).append("\t").append(p_equalVariants).append("\t");
                        sb.append(t_notEqualVariants).append("\t").append(p_notEqualVariants);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static List<double[][]> sublist(List<double[][]> loadList, TIntArrayList indexList){
        indexList.sort();
        List<double[][]> res=new ArrayList<>();
        for (int i = 0; i < indexList.size(); i++) {
            res.add(loadList.get(indexList.get(i)));
        }
        return res;
    }

    /**
     *
     * @param vcfLine
     * @return 返回所有个体的load, load 是二维数组，
     * [slightly/strongly][additive/dominance], -1 if NA
     */
    private static List<double[][]> calculateLoad(String vcfLine){
        List<double[][]> res=new ArrayList<>();
        List<String> temp=PStringUtils.fastSplit(vcfLine);
        List<String> tem, te;
        double[][] slightlyStronglyABDLoad; // dim 1 2: slightlyStrongly ABD
        double[][] slightlyStronglyAdditiveDominanceLoad; //dim 1 2: slightlyStrongly additiveDominance
        double additive, dominance;
        for (int i = 3; i < temp.size(); i++) {
            tem=PStringUtils.fastSplit(temp.get(i), "|");
            slightlyStronglyABDLoad=new double[2][];
            for (int j = 0; j < slightlyStronglyABDLoad.length; j++) {
                slightlyStronglyABDLoad[j]=new double[3];
                Arrays.fill(slightlyStronglyABDLoad[j], -1);
            }
            for (int j = 0; j <tem.size(); j++) {
                te=PStringUtils.fastSplit(tem.get(j), ",");
                for (int k = 0; k < te.size(); k++) {
                    if (!StringTool.isNumeric(te.get(k))) continue;
                    slightlyStronglyABDLoad[j][k]=Double.parseDouble(te.get(k));
                }
            }
            slightlyStronglyAdditiveDominanceLoad=new double[2][];
            for (int j = 0; j < slightlyStronglyAdditiveDominanceLoad.length; j++) {
                slightlyStronglyAdditiveDominanceLoad[j]=new double[2];
                Arrays.fill(slightlyStronglyAdditiveDominanceLoad[j], -1);
            }
            for (int j = 0; j < slightlyStronglyABDLoad.length; j++) {
                if (Arrays.stream(slightlyStronglyABDLoad[j]).min().getAsDouble()<0) continue;
                additive=Arrays.stream(slightlyStronglyABDLoad[j]).sum();
                dominance=Arrays.stream(slightlyStronglyABDLoad[j]).min().getAsDouble();
                slightlyStronglyAdditiveDominanceLoad[j][0]=additive;
                slightlyStronglyAdditiveDominanceLoad[j][1]=dominance;
            }
            res.add(slightlyStronglyAdditiveDominanceLoad);
        }
        return res;
    }
}
