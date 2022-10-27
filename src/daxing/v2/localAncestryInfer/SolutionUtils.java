package daxing.v2.localAncestryInfer;

import it.unimi.dsi.fastutil.ints.*;

import java.util.*;

public class SolutionUtils {

    /**
     *
     * @param srcGenotype the first dim is haplotype, the second dim is SNP position
     * @param queryGenotype
     * @param switchCostScore
     * @return mini cost score matrix, the first dim is haplotype, the second dim is SNP position
     */
    public static double[][] getMiniCostScore(double[][] srcGenotype, double[] queryGenotype,
                                              double switchCostScore){
        //        double switchCostScore= 1.5;
//        int[][] srcGenotype = {{0,1,0,1,0,1,0,0,0,0,1,1},
//                            {0,0,0,1,0,1,1,0,0,0,1,1},
//                            {0,0,1,0,1,0,0,0,1,0,1,1},
//                            {0,0,0,0,1,0,1,0,1,1,1,1},
//                            {1,1,0,0,0,0,1,1,1,1,0,0},
//                            {1,0,0,1,0,0,1,1,1,1,0,0}};
//        int[] queryGenotype =       {1,1,0,0,0,1,0,0,1,1,1,1};

//        long start = System.nanoTime();

        // distance
        double[][] distance = new double[srcGenotype.length][];
        for (int i = 0; i < distance.length; i++) {
            distance[i]= new double[srcGenotype[i].length];
            Arrays.fill(distance[i], -1);
        }
        for (int i = 0; i < distance.length; i++) {
            for (int j = 0; j < distance[i].length; j++) {
                distance[i][j]=Math.abs(srcGenotype[i][j]- queryGenotype[j]);
            }
        }

        // initialize mini cost score
        double[][] miniCost = new double[distance.length][];
        for (int i = 0; i < miniCost.length; i++) {
            miniCost[i] = new double[distance[0].length];
            Arrays.fill(miniCost[i], -1);
            miniCost[i][0] = distance[i][0];
        }


        // i is SNP position
        // j is haplotype index of source population
        // miniCost
        for (int i = 1; i < distance[0].length; i++) {

            // j-1 SNP位置，单倍型路径发生switch对应的最小Cost
            double miniCostSwitch=Double.MAX_VALUE;
            for (int j = 0; j < distance.length; j++) {
                miniCostSwitch = miniCost[j][i-1] < miniCostSwitch ? miniCost[j][i-1] : miniCostSwitch;
            }

            for (int j = 0; j < distance.length; j++) {
                // 最小cost路径对应当前haplotype
                if (miniCost[j][i-1] < miniCostSwitch+switchCostScore){
                    miniCost[j][i] = miniCost[j][i-1] + distance[j][i];
                }else {
                    // 最小cost路径对应转换单倍型
                    miniCost[j][i] = miniCostSwitch+switchCostScore+distance[j][i];
                }
            }
        }
//        System.out.println("calculate mini cost matrix take "+Benchmark.getTimeSpanSeconds(start)+ " seconds");
        return miniCost;
    }

    /**
     * 考虑单倍型所属的群体, 进行罚分
     * @param srcGenotype the first dim is haplotype, the second dim is SNP position
     * @param queryGenotype genotype
     * @param switchCostScore λ
     * @param srcIndiList 单倍型的样本名称
     * @param taxaSourceMap 样本对应群体的映射
     * @return mini cost score matrix, the first dim is haplotype, the second dim is SNP position
     */
    public static double[][] getMiniCostScore(double[][] srcGenotype, double[] queryGenotype,
                                              double switchCostScore,
                                              List<String> srcIndiList,
                                              Map<String, WindowSource.Source> taxaSourceMap){
//        long start = System.nanoTime();
        // distance
        double[][] distance = new double[srcGenotype.length][];
        for (int i = 0; i < distance.length; i++) {
            distance[i]= new double[srcGenotype[i].length];
            Arrays.fill(distance[i], -1);
        }
        for (int i = 0; i < distance.length; i++) {
            for (int j = 0; j < distance[i].length; j++) {
                distance[i][j]=Math.abs(srcGenotype[i][j]- queryGenotype[j]);
            }
        }

        // initialize mini cost score
        double[][] miniCost = new double[distance.length][];
        for (int i = 0; i < miniCost.length; i++) {
            miniCost[i] = new double[distance[0].length];
            Arrays.fill(miniCost[i], -1);
            miniCost[i][0] = distance[i][0];
        }


        // i is SNP position
        // j is haplotype index of source population
        // miniCost

        IntSet miniCostSwitchIndexSet;
        for (int i = 1; i < distance[0].length; i++) {

            // j-1 SNP位置，单倍型路径发生switch对应的最小Cost
            double miniCostSwitch=Double.MAX_VALUE;
            for (int j = 0; j < distance.length; j++) {
                miniCostSwitch = miniCost[j][i-1] < miniCostSwitch ? miniCost[j][i-1] : miniCostSwitch;
            }

            // j-1 SNP位置 单倍型路径发生switch对应的最小Cost index
            miniCostSwitchIndexSet = new IntOpenHashSet();
            for (int j = 0; j < distance.length; j++) {
                if (miniCostSwitch==miniCost[j][i-1]){
                    miniCostSwitchIndexSet.add(j);
                }
            }

            for (int j = 0; j < distance.length; j++) {
                // 最小cost路径对应当前haplotype
                if (miniCost[j][i-1] < miniCostSwitch+switchCostScore){
                    miniCost[j][i] = miniCost[j][i-1] + distance[j][i];
                }else {
                    // 最小cost路径对应转换单倍型

                    // i-1 SNP和当前SNP所属的单倍型是否来自同一群体
                    // 转换单倍型中i-1 SNP和当前SNP是否有在同一群体的案例
                    boolean ifExistAnyPairBelongToSamePop=false;
                    for (int index:miniCostSwitchIndexSet){
                        if(taxaSourceMap.get(srcIndiList.get(index)).equals(taxaSourceMap.get(srcIndiList.get(j)))){
                            ifExistAnyPairBelongToSamePop = true;
                            break;
                        }
                    }

                    // || taxaSourceMap.get(srcIndiList.get(j)).equals(WindowSource.Source.NONE)
                    if (ifExistAnyPairBelongToSamePop){
                        // i-1 SNP和当前SNP i所属的单倍型有来自同一群体 λ
                        // λ
                        miniCost[j][i] = miniCostSwitch+switchCostScore+distance[j][i];
                    }else {
                        // i-1 SNP和当前SNP所属的单倍型来自不同群体
                        // λ+1
                        miniCost[j][i] = miniCostSwitch+(switchCostScore+1)+distance[j][i];
                    }

                }
            }
        }

//        System.out.println("calculate mini cost matrix take "+Benchmark.getTimeSpanSeconds(start)+ " seconds");

        return miniCost;
    }

    /**
     *
     * @param miniCost the first dim is haplotype, the second dim is SNP position
     * @param switchCostScore λ
     * @return candidate solution, the first dim is mini cost sum haplotype index, the second dim is SNP position
     */
    public static IntSet[][] getCandidateSolution(double[][] miniCost, double switchCostScore){

//        long start = System.nanoTime();
        int haplotypeLen=miniCost[0].length;

        // miniCost score indexList
        double miniCostScore= Double.MAX_VALUE;
        IntList miniCostScoreIndexList = new IntArrayList();
        for (int i = 0; i < miniCost.length; i++) {
            miniCostScore = miniCost[i][haplotypeLen-1] < miniCostScore ? miniCost[i][haplotypeLen-1] :miniCostScore;
        }
        for (int i = 0; i < miniCost.length; i++) {
            if (miniCostScore==miniCost[i][haplotypeLen-1]){
                miniCostScoreIndexList.add(i);
            }
        }

        // new solution, the first dim is miniCostScoreList, the second dim is SNP
        IntSet[][] solution = new IntSet[miniCostScoreIndexList.size()][];
        for (int i = 0; i < solution.length; i++) {
            solution[i] = new IntSet[haplotypeLen];
            for (int j = 0; j < solution[i].length; j++) {
                solution[i][j] = new IntOpenHashSet();
            }
        }
        for (int i = 0; i < solution.length; i++) {
            solution[i][haplotypeLen-1].add(miniCostScoreIndexList.getInt(i));
        }


        // find all solution
        IntSet currentIndexSet;
        IntIterator tIntIterator;
        int index;
        for (int i = 0; i < solution.length; i++) {
            index = miniCostScoreIndexList.getInt(i);
            currentIndexSet = new IntOpenHashSet();
            currentIndexSet.add(index);
            for (int j = haplotypeLen - 1; j > 0; j--) {
                tIntIterator = currentIndexSet.iterator();
                while (tIntIterator.hasNext()){
                    index = tIntIterator.nextInt();

                    // 当前单倍型
                    if (miniCost[index][j-1] <= miniCost[index][j]){
                        solution[i][j-1].add(index);
                    }

                    // 转换单倍型
                    for (int k = 0; k < miniCost.length; k++) {
                        if (k==index) continue;
                        if ((miniCost[k][j-1]+switchCostScore) <= miniCost[index][j]){
                            solution[i][j-1].add(k);
                        }
                    }
                }
                currentIndexSet = solution[i][j-1];
            }
        }

//        System.out.println("calculate candidate solution take "+Benchmark.getTimeSpanSeconds(start)+" seconds");
        return solution;
    }

    public static EnumSet<WindowSource.Source>[][] getCandidateSourceSolutionEnumSet(double[][] miniCost,
                                                                                     double switchCostScore,
                                                                                     List<String> srcIndiList,
                                                                                     Map<String, WindowSource.Source> taxaSourceMap){
        int haplotypeLen=miniCost[0].length;

        IntSet[][] solution = GenotypeTable.getCandidateSolution(miniCost, switchCostScore);

        EnumSet<WindowSource.Source>[][] solutionSource= new EnumSet[solution.length][];
        for (int i = 0; i < solutionSource.length; i++) {
            solutionSource[i] = new EnumSet[haplotypeLen];
            for (int j = 0; j < solutionSource[i].length; j++) {
                solutionSource[i][j] = EnumSet.noneOf(WindowSource.Source.class);
            }
        }

        for (int i = 0; i < solutionSource.length; i++) {
            for (int j = 0; j < solutionSource[i].length; j++) {
                for (int ele:solution[i][j]){
                    solutionSource[i][j].add(taxaSourceMap.get(srcIndiList.get(ele)));
                }
            }
        }
        return solutionSource;
    }


    /**
     *
     * @param miniCost the first dim is haplotype, the second dim is SNP position
     * @param switchCostScore λ
     * @param srcIndiList 单倍型的样本名称
     * @param taxaSourceMap 样本对应群体的映射
     * @return candidate source solution, the first dim is mini cost sum haplotype index, the second dim is SNP position
     */
    public static List<WindowSource.Source>[][] getCandidateSourceSolutionList(double[][] miniCost, double switchCostScore,
                                                                               List<String> srcIndiList,
                                                                               Map<String, WindowSource.Source> taxaSourceMap){
//        long start = System.nanoTime();
        int haplotypeLen=miniCost[0].length;

        EnumSet<WindowSource.Source>[][] solutionSource= GenotypeTable.getCandidateSourceSolutionEnumSet(miniCost,
                switchCostScore, srcIndiList, taxaSourceMap);

        List<WindowSource.Source>[][] solutionSourceList = new List[solutionSource.length][];
        for (int i = 0; i < solutionSourceList.length; i++) {
            solutionSourceList[i] = new List[haplotypeLen];
            for (int j = 0; j < solutionSource[i].length; j++) {
                solutionSourceList[i][j]=new ArrayList<>(solutionSource[i][j]);
                Collections.sort(solutionSourceList[i][j]);
            }
        }


        // 从单倍型solution到Source solution的转变过程中, 会有重复的solution
        Set<List<List<WindowSource.Source>>> set= new HashSet<>();
        List<List<WindowSource.Source>> list;
        for (int i = 0; i < solutionSourceList.length; i++) {
            list = new ArrayList<>();
            for (int j = 0; j < solutionSourceList[i].length; j++) {
                list.add(solutionSourceList[i][j]);
            }
            set.add(list);
        }

        List<List<List<WindowSource.Source>>> listRemovedDup= new ArrayList<>(set);
        List<WindowSource.Source>[][] solutionSourceListRemovedDup= new List[listRemovedDup.size()][];
        for (int i = 0; i < solutionSourceListRemovedDup.length; i++) {
            solutionSourceListRemovedDup[i] = new List[haplotypeLen];
            for (int j = 0; j < solutionSourceListRemovedDup[i].length; j++) {
                solutionSourceListRemovedDup[i][j] = listRemovedDup.get(i).get(j);
                Collections.sort(solutionSourceListRemovedDup[i][j]);
            }
        }

//        System.out.println("calculate candidate solution take "+Benchmark.getTimeSpanSeconds(start)+" seconds");
        return solutionSourceListRemovedDup;
    }


}
