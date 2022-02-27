package xiaohan.rareallele;

import daxing.load.ancestralSite.Standardization;
import pgl.infra.utils.PStringUtils;
import xiaohan.utils.IOUtils;
import xiaohan.utils.RowTable;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.List;


/**
 * @ author: yxh
 * @ created: 2022-02-21 : 10:32 AM
 */
public class Triads {

    HashMap<String, String> triadsMapAB;
    HashMap<String, String> triadsMapAD;
    RowTable ExpressionTable;
    RowTable RVCountTable;
    HashMap<String, Integer> expressGeneMap;
    HashMap<String, Integer> RVcountGeneMap;
    double ExpgeneA;
    double ExpgeneB;
    double ExpgeneD;
    int CountgeneA;
    int CountgeneB;
    int CountgeneD;
    String ABD;

    public Triads(String[] args) {
//        System.out.println("=================initiating triads==============");
//        this.getTriads(args[0]);
//        System.out.println("=================initiating express==============");
//        this.initExpressionTable(args[1]);
//        System.out.println("=================initiating count==============");
//        this.initRVCountTable(args[2]);
//        System.out.println("=================writing files==============");
//        this.outFile(args[3]);
        this.genopheno(args[0], args[1]);

    }

    public void genopheno(String infile, String outfile) {
        BufferedReader br = IOUtils.getInFile(new File(infile).getAbsolutePath());
        BufferedWriter bw = IOUtils.getOutFile(new File(outfile).getAbsolutePath());
        String temp = null;
        List<String> temps;
        double[] ABD = new double[3];
        String pattern = null;
        StringBuilder sb = new StringBuilder();
        try {
            while ((temp = br.readLine()) != null) {
                temps = PStringUtils.fastSplit(temp);
                for (int i = 0; i < ABD.length; i++) {
                    ABD[i] = Double.parseDouble(temps.get(i+1));
                }
                pattern = Standardization.getNearestPointIndex(ABD).getRegion();
                sb.setLength(0);
                sb.append(temps.get(0) + "\t" + pattern + "\n");
                bw.write(sb.toString());
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void initExpressionTable(String expression) {
        ExpressionTable = new RowTable(new File(expression).getAbsolutePath());
        expressGeneMap = new HashMap<>();
        for (int i = 0; i < ExpressionTable.getRowNumber(); i++) {
            expressGeneMap.put(String.valueOf(ExpressionTable.getCell(i, 0)), i);
        }
    }

    private void initRVCountTable(String RVcount) {
        RVCountTable = new RowTable(new File(RVcount).getAbsolutePath());
        RVcountGeneMap = new HashMap<>();
        for (int i = 0; i < RVCountTable.getRowNumber(); i++) {
            RVcountGeneMap.put(String.valueOf(RVCountTable.getCell(i, 0)), i);
        }
    }

    private void outFile(String outfile) {
        BufferedWriter bw = IOUtils.getOutFile(outfile);
        StringBuilder sb = new StringBuilder();
        DecimalFormat dec = new DecimalFormat("0.00000");
        try {
            for (int i = 1; i <= 360; i++) {
                String sample = "E" + PStringUtils.getNDigitNumber(3, i);
                for (String geneA : triadsMapAB.keySet()) {
                    String geneB = triadsMapAB.get(geneA);
                    String geneD = triadsMapAD.get(geneA);
                    ExpgeneA = getRVexpress(geneA, sample);
                    ExpgeneB = getRVexpress(geneB, sample);
                    ExpgeneD = getRVexpress(geneD, sample);
                    if (ExpgeneD == -1 || ExpgeneB == -1 || ExpgeneA == -1) continue;
                    CountgeneA = getRVcount(geneA, sample);
                    CountgeneB = getRVcount(geneB, sample);
                    CountgeneD = getRVcount(geneD, sample);
                    if (CountgeneA == -1 || CountgeneB == -1 || CountgeneD == -1)
                        continue;
                    ABD = String.valueOf(CountgeneA) + "_" + String.valueOf(CountgeneB) + "_" + String.valueOf(CountgeneD);
                    double all = ExpgeneA + ExpgeneB + ExpgeneD;
                    if (all == 0) continue;
                    sb.setLength(0);
//                    sb.append(ABD + "\t" + ExpgeneA + "\t" + ExpgeneB + "\t" + ExpgeneD + "\t" + ExpgeneA/all + "\t" + ExpgeneB/all + "\t" + ExpgeneD/all);
                    sb.append(ABD + "\t" + dec.format(ExpgeneA / all) + "\t" + dec.format(ExpgeneB / all) + "\t" + dec.format(ExpgeneD / all));
                    bw.write(sb.toString() + "\n");
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private double getRVexpress(String GeneName, String SampleName) {
        if (expressGeneMap.containsKey(GeneName) && ExpressionTable.getColumnIndex(SampleName) != -1) {
            return Double.parseDouble(String.valueOf(ExpressionTable.getCell(expressGeneMap.get(GeneName), ExpressionTable.getColumnIndex(SampleName))));
        } else return -1;
    }

    private int getRVcount(String GeneName, String SampleName) {
        if (!RVcountGeneMap.containsKey(GeneName) || RVCountTable.getColumnIndex(SampleName) == -1) {
            return -1;
        } else {
//            System.out.println(GeneName);
            int res = RVCountTable.getCellAsInteger(RVcountGeneMap.get(GeneName), RVCountTable.getColumnIndex(SampleName));
            if (res > 0) {
                return 1;
            } else {
                return 0;
            }
        }
    }

    private void getTriads(String infile) {
        BufferedReader br = xiaohan.utils.IOUtils.getInFile(new File(infile).getAbsolutePath());
        String temp = null;
        List<String> temps;
        triadsMapAB = new HashMap<>();
        triadsMapAD = new HashMap<>();
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("A")) continue;
                temps = PStringUtils.fastSplit(temp);
                triadsMapAB.put(temps.get(0), temps.get(1));
                triadsMapAD.put(temps.get(0), temps.get(2));
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}