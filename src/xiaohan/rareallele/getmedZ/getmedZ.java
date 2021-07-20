package xiaohan.rareallele.getmedZ;

import pgl.infra.utils.PStringUtils;
import xiaohan.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

/**
 * @ author: yxh
 * @ created: 2021-07-20 : 10:37 PM
 */
public class getmedZ {

    public void getmedZ() {
        String infile = "/data2/xiaohan/pheno/peerless/all_normalized_expression.txt";
//        String infile = "/Users/yxh/Documents/RareAllele/script/001pheno_processing/all_normalized_expression.txt";
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedReader br1 = IOUtils.getTextReader(infile);
        String outfile = "/data2/xiaohan/pheno/peerless/outliers_medz_picked.txt";
        String outfilecount = "/data2/xiaohan/pheno/peerless/outliers_medz_counts.txt";
//        String outfile = "/Users/yxh/Documents/RareAllele/script/001pheno_processing/outliers_medz_picked.txt";
        String temp = null;
        List<String> temps = null;
        HashSet<String> tissueSet = new HashSet<>();
        HashSet<String> geneSet = new HashSet<>();
        StringBuilder sb = new StringBuilder();
        int countline1 = 0;
        int countline = 0;
        int columnnumber = 0;
        String error = null;
        String error1 = null;
        //picked variables
        String gene = null;
        String name = null;
        String tissue = null;
        String max = null;
        try {
            String[] header = null;
            String headers = null;
            System.out.println("Starting initiation ----------------------------------");
            while ((temp = br.readLine()) != null) {
                countline++;
                if (countline % 5000 == 0) {
                    System.out.println(countline);
                }
                temps = PStringUtils.fastSplit(temp);
                if (temp.startsWith("Tissue")) {
                    header = new String[temps.size() - 2];
                    for (int i = 2; i < temps.size(); i++) {
                        header[i - 2] = temps.get(i);
                    }
                    sb.setLength(0);
                    sb.append("GENE\t");
                    for (int i = 0; i < header.length; i++) {
                        sb.append(header[i] + "\t");
                    }
                    headers = sb.toString().replaceAll("\\s+$", "");
                    columnnumber = temps.size() - 2;
                    System.out.println(columnnumber);
                    continue;
                }
                tissueSet.add(temps.get(0));
                geneSet.add(temps.get(1));
            }
            br.close();
            String[][][] expression = new String[geneSet.size()][columnnumber][tissueSet.size()];
            String[] tissues = tissueSet.toArray(new String[0]);
            String[] genes = geneSet.toArray(new String[0]);
            HashMap<String, Integer> tissuesIndex = new HashMap<>();
            HashMap<String, Integer> genesIndex = new HashMap<>();
            for (int i = 0; i < tissues.length; i++) {
                tissuesIndex.put(tissues[i], i);
                for (int j = 0; j < genes.length; j++) {
                    genesIndex.put(genes[j], j);
                    for (int k = 0; k < columnnumber; k++) {
                        expression[j][k][i] = "NA";
                    }
                }
            }
            System.out.println("Starting building expression file ----------------------------------");
            while ((temp = br1.readLine()) != null) {
                countline1++;
                if (countline1 % 5 == 0) {
//                    System.out.println(countline1);
                }
                temps = PStringUtils.fastSplit(temp);
                if (temp.startsWith("Tissue")) continue;
                for (int i = 0; i < columnnumber; i++) {
                    expression[genesIndex.get(temps.get(1))][i][tissuesIndex.get(temps.get(0))] = temps.get(i + 2);
                }
            }
            System.out.println("Finished building expression file !");
            br1.close();
            String[][] expressionmedian = new String[geneSet.size()][columnnumber];
            String[][] expressionN = new String[geneSet.size()][columnnumber];
            String[] array = new String[tissues.length];
            int len = -1;
            double median = -1;
            System.out.println("Starting calculating n and median ----------------------------------");
            for (int i = 0; i < genes.length; i++) {
                for (int j = 0; j < columnnumber; j++) {
                    array = expression[i][j];
                    expressionN[i][j] = medZUtils.medn(array);
                    System.out.println(i + "\t" + j);
                    System.out.println(genes[i] + "\t" + header[j]);
                    expressionmedian[i][j] = String.valueOf(medZUtils.medmedian(array));
                }
            }
            // writing medz picked file
            BufferedWriter bw = IOUtils.getTextWriter(outfile);
            System.out.println("Starting picking outliers ----------------------------------");
            for (int i = 0; i < genes.length; i++) {
                sb.setLength(0);
                sb.append(genes[i] + "\t");
                String[] outliers = medZUtils.pickoutlier(expressionmedian[i]);
//                System.out.println(genes[i] + "_" + expressionmedian[i][0] +"_"+ outliers[0] +"_"+ outliers[1]);
                sb.append(header[Integer.parseInt(outliers[1])] + "\t");
                sb.append(expressionN[i][Integer.parseInt(outliers[1])] + "\t");
                sb.append(outliers[0]);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            // writing medz counts file
            BufferedWriter bwcount = IOUtils.getTextWriter(outfilecount);
            bwcount.write(headers + "\n");
            for (int i = 0; i < genes.length; i++) {
                sb.setLength(0);
                sb.append(genes[i] + "\t");
                for (int j = 0; j < columnnumber; j++) {
                    sb.append(expressionN[i][j] + "\t");
                }
                bwcount.write(sb.toString().replaceAll("\\s+$", "") + "\n");
            }
            bwcount.flush();
            bwcount.close();
            // writing single counts file and picked file
            String[] array1 = new String[columnnumber];
            String[] outliers = new String[2];
            BufferedWriter[] bw1 = new BufferedWriter[tissues.length];
            BufferedWriter[] bw2 = new BufferedWriter[tissues.length];
            for (int i = 0; i < tissues.length; i++) {
                bw1[i] = IOUtils.getTextWriter(new File("/data2/xiaohan/pheno/peerless/" + tissues[i] + "_outliers_counts.txt").getAbsolutePath());
                bw2[i] = IOUtils.getTextWriter(new File("/data2/xiaohan/pheno/peerless/" + tissues[i] + "_outliers_singlez_picked.txt").getAbsolutePath());
                bw1[i].write(headers + "\n");
                bw2[i].write("GENE\tIND\tDFS\tZ\n");
            }
            for (int i = 0; i < genes.length; i++) {
                for (int j = 0; j < tissues.length; j++) {
                    gene = genes[i];
                    sb.setLength(0);
                    sb.append(gene + "\t");
                    for (int k = 0; k < columnnumber; k++) {
                        array1[k] = expression[i][k][j];
                        if (array1[k].equals("NA")) {
                            sb.append("0\t");
                        } else {
                            sb.append("1\t");
                        }
                    }
                    outliers = medZUtils.pickoutlier(array1);
                    if (outliers[0].equals("NA") && outliers[1].equals("NA")) continue;
                    name = header[Integer.parseInt(outliers[1])];
                    tissue = tissues[j];
                    max = expression[i][Integer.parseInt(outliers[1])][j];
                    bw1[j].write(sb.toString().replaceAll("\\s+$", "") + "\n");
                    bw2[j].write(gene + "\t" + name + "\t" + tissue + "\t" + max + "\n");
                }
            }
            for (int i = 0; i < tissues.length; i++) {
                bw1[i].flush();
                bw1[i].close();
                bw2[i].flush();
                bw2[i].close();
            }
        } catch (Exception e) {
            e.printStackTrace();
//            System.out.println(countline1);
//            System.out.println(countline);
//            System.out.println(columnnumber);
            System.out.println(error);
            System.out.println(error1);
        }
    }
}
