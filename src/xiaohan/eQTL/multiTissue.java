package xiaohan.eQTL;

import xiaohan.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class multiTissue {

    public multiTissue(String[] args) {
        this.getTissueSharing(args);
    }

    public void getTissueSharing(String[] args) {
        try {
            BufferedReader br[] = new BufferedReader[args.length];
            String[] tissues = new String[args.length];
            String subgenome = null;
            for (int i = 0; i < tissues.length; i++) {
                tissues[i] = args[i].split("/")[args[i].split("/").length-2];
                subgenome = args[i].split("/")[args[i].split("/").length-1].split("\\.")[0];
            }
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < tissues.length; i++) {
                sb.append(tissues[i]+"_");
            }
            String outputDir = "/data2/xiaohan/tensorResult/tissueShared";
            for (int i = 0; i < args.length; i++) {
                if (args[i].endsWith("gz")) {
                    br[i] = IOUtils.getTextGzipReader(args[i]);
                } else {
                    br[i] = IOUtils.getTextReader(args[i]);
                }
            }

            String temp = null;
            HashSet<String> totalnameSet = new HashSet<>();
            HashMap<String, Integer> nameNumber = new HashMap<>();
            HashMap<String,Integer> geneIndex = new HashMap<>();
            for (int i = 0; i < args.length; i++) {
                while ((temp = br[i].readLine()) != null) {
                    if (!totalnameSet.contains(temp)) {
                        totalnameSet.add(temp);
                    }
                }
                br[i].close();
            }

            String[] geneNames = totalnameSet.toArray(new String[0]);
            Arrays.sort(geneNames);
            for (int i = 0; i < geneNames.length; i++) {
                geneIndex.put(geneNames[i],i);
            }

            int[][] geneMatrix = new int[geneNames.length][args.length];
            for (int i = 0; i < geneNames.length; i++) {
                for (int j = 0; j < args.length; j++) {
                    geneMatrix[i][j] = 0;
                }
            }

            for (int i = 0; i < args.length; i++) {
                if (args[i].endsWith("gz")) {
                    br[i] = IOUtils.getTextGzipReader(args[i]);
                } else {
                    br[i] = IOUtils.getTextReader(args[i]);
                }
            }

            for (int i = 0; i <args.length;  i++) {
                while ((temp = br[i].readLine()) != null) {
                    int index = geneIndex.get(temp);
                    geneMatrix[index][i] = 1;
                }
                br[i].close();
            }

            for (int i = 0; i < geneNames.length; i++) {
                String gene = geneNames[i];
                int num = 0;
                for (int j = 0; j < args.length; j++) {
                    num += geneMatrix[i][j];
                }
                nameNumber.put(gene, num);
            }

            BufferedWriter bw[] = new BufferedWriter[args.length];
            for (int i = 0; i < args.length; i++) {
                int num = i + 1;
                bw[i] = IOUtils.getTextWriter(new File(outputDir, "SharedBy_" + num + "tissue_"+sb.toString()+subgenome+".txt").getAbsolutePath());
            }

            for (int i = 0; i < geneNames.length; i++) {
                String gene = geneNames[i];
                int index = nameNumber.get(gene) - 1;
                bw[index].write(gene + "\n");
            }

            for (int i = 0; i < args.length; i++) {
                bw[i].flush();
                bw[i].close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        new multiTissue(args);
    }
}
