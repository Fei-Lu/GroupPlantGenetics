package daxing.load.complementary;

import com.ibm.icu.text.NumberFormat;
import daxing.common.IOTool;
import gnu.trove.list.array.TDoubleArrayList;
import org.apache.commons.lang3.math.NumberUtils;
import pgl.PGLConstraints;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.*;

/**
 * this class is for evaluating the W (one sample) statistic at the traids level
 */
public class WilcoxonRankUtils {

    // individualsLoadList contains hexaploid and pseudohexaploid
    // dim1: AdditiveOrDominance; dim2: SlightlyOrStrongly; dim3: triads
    // triads num is 4, T000072,T001133,T003691,T005781
    List<double[][][]> individualsLoadList;
    List<Group> groupList;
    List<String> taxaNameList;

    public WilcoxonRankUtils(String triadsLoadHexaploidPseudoFile){
        individualsLoadList=new ArrayList<>();
        groupList=new ArrayList<>();
        taxaNameList=new ArrayList<>();
        try (BufferedReader br = IOTool.getReader(triadsLoadHexaploidPseudoFile)) {
            String line, taxaName;
            List<String> temp;
            br.readLine();
            double[][][] load=null;
            Group group;
            AdditiveOrDominance additiveOrDominance;
            SlightlyOrStrongly slightlyOrStrongly;
            int flag=0;
            while ((line=br.readLine())!=null){
                flag++;
                if (flag%4==1){
                    load=new double[AdditiveOrDominance.values().length][][];
                    for (int i = 0; i < load.length; i++) {
                        load[i]=new double[SlightlyOrStrongly.values().length][];
                        for (int j = 0; j < load[i].length; j++) {
                            load[i][j]=new double[4];
                            Arrays.fill(load[i][j], -1);
                        }
                    }
                }
                temp= PStringUtils.fastSplit(line);
                taxaName=temp.get(0);
                group = temp.get(1).contains("_AT") ? Group.Pseudo : Group.Hexaploid;
                additiveOrDominance = AdditiveOrDominance.valueOf(temp.get(2).toUpperCase(Locale.ROOT));
                slightlyOrStrongly = SlightlyOrStrongly.valueOf(temp.get(3).toUpperCase());
                for (int i = 4; i < temp.size(); i++) {
                    if (!NumberUtils.isNumber(temp.get(i))) continue;
                    load[additiveOrDominance.index][slightlyOrStrongly.index][i-4]=NumberUtils.toDouble(temp.get(i));
                }
                if (flag%4==0){
                    individualsLoadList.add(load);
                    groupList.add(group);
                    taxaNameList.add(taxaName);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private enum Group {
        Hexaploid, Pseudo
    }

    private enum AdditiveOrDominance{
        ADDITIVE(0),DOMINANCE(1);

        int index;
        AdditiveOrDominance(int index) {
            this.index=index;
        }
    }

    private enum SlightlyOrStrongly{
        SLIGHTLY(0), STRONGLY(1);

        int index;
        SlightlyOrStrongly(int index) {
            this.index=index;
        }
    }

    private List<double[][][]> getPseudoLoad(){
        List<double[][][]> pseudoLoad=new ArrayList<>();
        for (int i = 0; i < individualsLoadList.size(); i++) {
            if (groupList.get(i)!=Group.Pseudo) continue;
            pseudoLoad.add(individualsLoadList.get(i));
        }
        return pseudoLoad;
    }

    /**
     *
     * @return WilcoxonSignedRank of all hexaploid and pseudo
     * dim1: AdditiveOrDominance; dim2: SlightlyOrStrongly; dim3: triads
     */
    private List<double[][][]> calculateAllIndividualParallel(){
        List<Callable<double[][][]>> callableList=new ArrayList<>();
        List<double[][][]> pseudoLoadList=getPseudoLoad();
        for (double[][][] individualLoad:individualsLoadList){
            callableList.add(()->calculateWilcoxonSignedRank(individualLoad, pseudoLoadList));
        }
        ExecutorService executorService= Executors.newFixedThreadPool(PGLConstraints.parallelLevel);
        try {
            List<Future<double[][][]>> futureList=executorService.invokeAll(callableList);
            executorService.shutdown();
            executorService.awaitTermination(Long.MAX_VALUE, TimeUnit.MICROSECONDS);
            List<double[][][]> res=new ArrayList<>();
            for (Future<double[][][]> future : futureList) {
                res.add(future.get());
            }
            return res;
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (ExecutionException e) {
            e.printStackTrace();
        }
        return null;
    }

    public void writeWilcoxonSignedRank(String outFile){
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            StringBuilder sb=new StringBuilder();
            sb.append("TaxonID\tGroup\tAdditiveOrDominance\tSlightlyOrStrongly\tT000072\tT001133\tT003691\tT005781");
            bw.write(sb.toString());
            bw.newLine();
            List<double[][][]> wilcoxonSignedRank_AllIndis=this.calculateAllIndividualParallel();
            String taxaName;
            Group group;
            double[][][] wilcoxonSignedRank;
            AdditiveOrDominance additiveOrDominance;
            SlightlyOrStrongly slightlyOrStrongly;
            double[] triadsWilcoxonSignedRank;
            NumberFormat numberFormat=NumberFormat.getInstance();
            numberFormat.setGroupingUsed(false);
            numberFormat.setMinimumFractionDigits(3);
            for (int i = 0; i < wilcoxonSignedRank_AllIndis.size(); i++) {
                taxaName=taxaNameList.get(i);
                group=groupList.get(i);
                wilcoxonSignedRank=wilcoxonSignedRank_AllIndis.get(i);
                for (int j = 0; j < AdditiveOrDominance.values().length; j++) {
                    additiveOrDominance= j ==0 ? AdditiveOrDominance.ADDITIVE : AdditiveOrDominance.DOMINANCE;
                    for (int k = 0; k < SlightlyOrStrongly.values().length; k++) {
                        slightlyOrStrongly= k==0 ? SlightlyOrStrongly.SLIGHTLY : SlightlyOrStrongly.STRONGLY;
                        triadsWilcoxonSignedRank=wilcoxonSignedRank[j][k];
                        sb.setLength(0);
                        sb.append(taxaName).append("\t").append(group.name().toLowerCase()).append("\t");
                        sb.append(additiveOrDominance).append("\t").append(slightlyOrStrongly).append("\t");
                        for (int l = 0; l < 4; l++) {
                            if (Double.isNaN(triadsWilcoxonSignedRank[l])){
                                sb.append("NaN").append("\t");
                            }else {
                                sb.append(numberFormat.format(triadsWilcoxonSignedRank[l])).append("\t");
                            }
                        }
                        sb.deleteCharAt(sb.length()-1);
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

    /**
     *
     * @param hexaploidIndividualLoad individual load
     * @param pseudoLoadList all pseudo individuals
     * @return
     */
    private double[][][] calculateWilcoxonSignedRank(double[][][] hexaploidIndividualLoad,
                                                List<double[][][]> pseudoLoadList){
        double[][][] res=new double[AdditiveOrDominance.values().length][][];
        for (int i = 0; i < res.length; i++) {
            res[i]=new double[SlightlyOrStrongly.values().length][];
            for (int j = 0; j < res[i].length; j++) {
                res[i][j]=new double[4];
            }
        }
        double[][][][] transformedPseudoLoad=transformLoad(pseudoLoadList);
        double hexaploidLoad;
        double[] pseudoLoad;
        TDoubleArrayList pseudoLoadRemovedNAInfNaNList;
        Vmap2ComplementaryVCF.Statics statics= Vmap2ComplementaryVCF.Statics.WSR;
        for (int i = 0; i < AdditiveOrDominance.values().length; i++) {
            for (int j = 0; j < SlightlyOrStrongly.values().length; j++) {
                // triads num is 4
                for (int k = 0; k < 4; k++) {
                    hexaploidLoad=hexaploidIndividualLoad[i][j][k];
                    if (hexaploidLoad < 0) continue;
                    pseudoLoad=transformedPseudoLoad[i][j][k];
                    pseudoLoadRemovedNAInfNaNList=new TDoubleArrayList();
                    for (double v:pseudoLoad){
                        if (v < 0) continue;
                        if (Double.isNaN(v)) continue;
                        if (Double.isInfinite(v)) continue;
                        pseudoLoadRemovedNAInfNaNList.add(v);
                    }
                    if (pseudoLoadRemovedNAInfNaNList.size() < 2) continue;
                    res[i][j][k]=statics.getHexaploidStaticsValue(hexaploidLoad, pseudoLoadRemovedNAInfNaNList.toArray());
                }
            }
        }
        return res;
    }

    private double[][][][] transformLoad(List<double[][][]> pseudoLoadList){
        double[][][][] res=new double[AdditiveOrDominance.values().length][][][];
        for (int i = 0; i < res.length; i++) {
            res[i]=new double[SlightlyOrStrongly.values().length][][];
            for (int j = 0; j < res[i].length; j++) {
                res[i][j]=new double[4][];
                for (int k = 0; k < res[i][j].length; k++) {
                    res[i][j][k]=new double[pseudoLoadList.size()];
                }
            }
        }
        for (int i = 0; i < pseudoLoadList.size(); i++) {
            for (int j = 0; j < pseudoLoadList.get(i).length; j++) {
                for (int k = 0; k < pseudoLoadList.get(i)[j].length; k++) {
                    for (int l = 0; l < pseudoLoadList.get(i)[j][k].length; l++) {
                        res[j][k][l][i]=pseudoLoadList.get(i)[j][k][l];
                    }
                }
            }
        }
        return res;
    }
}
