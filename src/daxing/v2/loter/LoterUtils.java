package daxing.v2.loter;

import daxing.common.factors.HexaploidBySubcontinent;
import daxing.common.table.RowTableTool;
import daxing.common.utiles.ArrayTool;
import daxing.common.utiles.IOTool;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

/**
 * this class contain statics methods calculating individual local ancestry size and individual site percent
 */
public class LoterUtils {
    /**
     * AB: WE DE FTT
     * D: ATF ATN
     */
    private enum AncestryNumber{
        A(3), B(3), D(2);

        int ancestryNum;

        AncestryNumber(int ancestryNum){
            this.ancestryNum=ancestryNum;
        }

        public int getAncestryNum() {
            return ancestryNum;
        }
    }

    public enum GroupInTaxaInfo{
        Subspecies_by6_TreeValidated(10, "LRCL"), GroupbyContinent(36,"GroupbyContinent"), IfAnalysisByPloidy(11,"Hexaploid");

        int index;
        String meaning;

        GroupInTaxaInfo(int index, String meaning){
            this.index=index;
            this.meaning = meaning;
        }

        public int getIndex() {
            return index;
        }

        public List<String> getGroup(){
            List<String> groupList = new ArrayList<>();
            if (this==Subspecies_by6_TreeValidated){
                groupList.add("Landrace");
                groupList.add("Cultivar");
            }else if (this == GroupbyContinent){
                for (HexaploidBySubcontinent hexaploidBySubcontinent: HexaploidBySubcontinent.values()){
                    groupList.add(hexaploidBySubcontinent.name());
                }
            } else if (this == IfAnalysisByPloidy){
                groupList.add("Hexaploid");
            }
            return groupList;
        }
    }

    public static void calculateIndividualSitePercent(String inputDir, String genotypeDir, String taxaInfo, GroupInTaxaInfo groupInTaxaInfo,
                                                      String outDir){
        List<File> laiFiles = IOTool.getFileListInDirEndsWith(inputDir, ".gz");
        List<File> vcfFiles = IOTool.getFileListInDirEndsWith(genotypeDir, "LRCL.recode.vcf.gz");
        String[] outNames =
                laiFiles.stream().map(File::getName).map(s -> s.replaceAll(".txt.gz",".ancestry.txt.gz")).toArray(String[]::new);
        Map<String, String> taxaGroupMap = RowTableTool.getMap(taxaInfo, 0, groupInTaxaInfo.getIndex());
        List<String> groupsList = groupInTaxaInfo.getGroup();
        IntStream.range(0, laiFiles.size()).forEach(e->{
            long start = System.nanoTime();
            try (BufferedReader br = IOTool.getReader(laiFiles.get(e));
                 BufferedReader brGenotype = IOTool.getReader(vcfFiles.get(e));
                 BufferedWriter bw =IOTool.getWriter(new File(outDir, outNames[e]))) {
                int chrID =Integer.parseInt(laiFiles.get(e).getName().substring(3,6));
                String chrSub = RefV1Utils.getChromosome(chrID, 1).substring(1,2);
                AncestryNumber ancestryNumber= AncestryNumber.valueOf(chrSub);
                String line;
                List<String> temp;
                List<TIntArrayList> individualSNPancestry = new ArrayList<>();
                List<String> taxaList= new ArrayList<>();
                TIntArrayList posList = new TIntArrayList();
                TIntArrayList posAncestry;
                boolean flag = true;
                while ((line=br.readLine())!=null){
                    temp =PStringUtils.fastSplit(line, " ");
                    posAncestry = new TIntArrayList();
                    for (String s: temp){
                        posAncestry.add(Integer.parseInt(s));
                    }
                    if (flag){
                        individualSNPancestry.add(posAncestry);
                        flag=false;
                    }else {
                        flag =true;
                    }
                }
                while ((line=brGenotype.readLine()).startsWith("##")){}
                temp = PStringUtils.fastSplit(line);
                for (int i = 9; i < temp.size(); i++) {
                    taxaList.add(temp.get(i));
                }
                while ((line=brGenotype.readLine())!=null){
                    temp = PStringUtils.fastSplit(line.substring(0, 30));
                    posList.add(Integer.parseInt(temp.get(1)));
                }
                double[][] ancestryProportion= new double[ancestryNumber.getAncestryNum()][];
                for (int i = 0; i < ancestryProportion.length; i++) {
                    ancestryProportion[i] = new double[individualSNPancestry.get(0).size()];
                    Arrays.fill(ancestryProportion[i], -1);
                }
                StringBuilder sb = new StringBuilder();
                sb.append("Chr\tPos\tGroup\tAncestry\tAncestryProportion");
                bw.write(sb.toString());
                bw.newLine();
                int[] count;
                NumberFormat numberFormat = NumberFormat.getInstance();
                numberFormat.setGroupingUsed(false);
                numberFormat.setMaximumFractionDigits(5);
                for (int i = 0; i < individualSNPancestry.get(0).size(); i++) {
                    for (int l = 0; l < groupsList.size(); l++) {
                        count = new int[ancestryNumber.getAncestryNum()];
                        for (int j = 0; j < individualSNPancestry.size(); j++) {
                            if (taxaGroupMap.get(taxaList.get(j)).equals(groupsList.get(l))){
                                for (int k = 0; k < count.length; k++) {
                                    if (individualSNPancestry.get(j).get(i) == k){
                                        count[k]++;
                                    }
                                }
                            }
                        }
                        double[] proportion = ArrayTool.getElementPercent(count);
                        for (int j = 0; j < proportion.length; j++) {
                            sb.setLength(0);
                            sb.append(chrID).append("\t").append(posList.get(i)).append("\t");
                            sb.append(groupsList.get(l)).append("\t");
                            sb.append(j).append("\t").append(numberFormat.format(proportion[j]));
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                    }
                }
                System.out.println(new File(outDir, outNames[e]).getName() + " completed in "+ Benchmark.getTimeSpanMinutes(start)+
                        " minutes");
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }
}
