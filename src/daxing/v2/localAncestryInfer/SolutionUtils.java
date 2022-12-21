package daxing.v2.localAncestryInfer;

import it.unimi.dsi.fastutil.ints.*;
import java.util.*;
import java.util.List;

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
     * @param srcGenotype
     * @param queryGenotype
     * @param switchCostScore
     * @param srcIndiList
     * @param taxaSourceMap
     * @return candidate solution, it consists of multiple solutions with the same mini cost score
     * dim1 is mini cost score index, dim2 is solution
     * Solution consists of multiple groups, every three numbers as a group, representing a tract
     * the first number is source population index, equal WindowSource.Source.index()
     * the second and third number is start(inclusive) position and end(inclusive) position
     */
    public static IntList[] getCandidateSolution(double[][] srcGenotype, double[] queryGenotype, double switchCostScore, List<String> srcIndiList,
                                                 Map<String, WindowSource.Source> taxaSourceMap){

        double[][] miniCost = SolutionUtils.getMiniCostScore(srcGenotype, queryGenotype, switchCostScore);

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
        IntList[] solutions = new IntList[miniCostScoreIndexList.size()];

        WindowSource.Source source;
        for (int i = 0; i < solutions.length; i++) {
            solutions[i] = new IntArrayList();
            source = taxaSourceMap.get(srcIndiList.get(miniCostScoreIndexList.getInt(i)));
            solutions[i].add(SourceType.valueOf(source.name()).getSourceFeature());
            solutions[i].add(haplotypeLen-1);
            solutions[i].add(haplotypeLen-1);
        }

        // find all solution
        IntSet currentIndexSet, nextIndexSet;
        int currentSourceFeature, nextSourceFeature;
        IntIterator tIntIterator;
        int index, currentSolutionElementIndex;

        for (int i = 0; i < solutions.length; i++) {
            index = miniCostScoreIndexList.getInt(i);
            currentIndexSet = new IntOpenHashSet();
            currentIndexSet.add(index);

            currentSolutionElementIndex=0;
            for (int j = haplotypeLen - 1; j > 0; j--) {
                tIntIterator = currentIndexSet.iterator();
                nextIndexSet = new IntOpenHashSet();
                while (tIntIterator.hasNext()){

                    index = tIntIterator.nextInt();

                    // 当前单倍型
                    if (miniCost[index][j-1] <= miniCost[index][j]){
                        nextIndexSet.add(index);
                    }

                    // 转换单倍型
                    for (int k = 0; k < miniCost.length; k++) {
                        if (k==index) continue;
                        if ((miniCost[k][j-1]+switchCostScore) <= miniCost[index][j]){
                            nextIndexSet.add(k);
                        }
                    }
                }

                currentSourceFeature = SolutionUtils.getSourceFutureFrom(currentIndexSet, srcIndiList, taxaSourceMap);
                nextSourceFeature = SolutionUtils.getSourceFutureFrom(nextIndexSet, srcIndiList, taxaSourceMap);
                if (currentSourceFeature==nextSourceFeature){
                    solutions[i].set(currentSolutionElementIndex*3+2, j-1);
                }else {
                    solutions[i].add(nextSourceFeature);
                    solutions[i].add(j-1);
                    solutions[i].add(j-1);
                    currentSolutionElementIndex++;
                }

                currentIndexSet = nextIndexSet;
            }
        }
        return solutions;
    }



    private static int getSourceFutureFrom(IntSet sourceIndexSet, List<String> srcIndiList,
                                                                 Map<String, WindowSource.Source> taxaSourceMap){
        EnumSet<WindowSource.Source> sourceEnumSet = EnumSet.noneOf(WindowSource.Source.class);
        for (int index : sourceIndexSet){
            sourceEnumSet.add(taxaSourceMap.get(srcIndiList.get(index)));
        }

        return SourceType.getSourceFeature(sourceEnumSet);
    }

    public static int getMiniOptimalSolutionSize(IntList[] solutions){
        int[] size = SolutionUtils.getOptimalSolutionsSize(solutions);
        int mini = Integer.MAX_VALUE;
        for (int i = 0; i < size.length; i++) {
            mini = size[i] < mini ? size[i] : mini;
        }
        return mini;
    }

    public static int[] getOptimalSolutionsSize(IntList[] solutions){
        int[] sizes = new int[solutions.length];
        Arrays.fill(sizes, -1);
        int cumSize=1;
        EnumSet<WindowSource.Source> sources;
        for (int i = 0; i < solutions.length; i++) {
            for (int j = solutions[i].size()-1; j > 0; j=j-3) {
                sources = SourceType.getSourcesFrom(solutions[i].getInt(j-2));
                if (sources.size()==1) continue;
                cumSize*=sources.size();
            }
            sizes[i] = cumSize;
            cumSize = 1;
        }
        return sizes;
    }

    public static int getMiniSolutionEleCount(IntList[] solutions){
        int miniSolutionEleCount = Integer.MAX_VALUE;
        for (int i = 0; i < solutions.length; i++) {
            miniSolutionEleCount = solutions[i].size() < miniSolutionEleCount ? solutions[i].size() :
                    miniSolutionEleCount;

        }
        return miniSolutionEleCount;
    }

    /**
     *
     * @param srcGenotype
     * @param queryGenotype
     * @param switchCostScore
     * @param srcIndiList
     * @param taxaSourceMap
     * @return BiDirectionCandidateSourceSolution, Forward candidate source solution, Reverse candidate solution
     */
    public static EnumMap<Direction, IntList[]> getBiDirectionCandidateSourceSolution(double[][] srcGenotype,
                                                                                      double[] queryGenotype,
                                                                                      double switchCostScore,
                                                                                      List<String> srcIndiList,
                                                                                      Map<String, WindowSource.Source> taxaSourceMap,
                                                                                      int maxSolutionCount){

        IntList[] forwardCandidateSolutionCurrent = SolutionUtils.getCandidateSolution(srcGenotype, queryGenotype, switchCostScore,
                        srcIndiList, taxaSourceMap);
        IntList[] forwardCandidateSolutionNext =
                SolutionUtils.getCandidateSolution(srcGenotype, queryGenotype, switchCostScore+1,
                        srcIndiList, taxaSourceMap);

        int totalSolutionSizeCurrent = SolutionUtils.getMiniOptimalSolutionSize(forwardCandidateSolutionCurrent);
        int totalSolutionSizeNext = SolutionUtils.getMiniOptimalSolutionSize(forwardCandidateSolutionNext);

        //  totalSolutionSizeCurrent < 0 是因为 两个Int相乘的结果大于Int max
        if ((totalSolutionSizeCurrent > maxSolutionCount && totalSolutionSizeNext <= totalSolutionSizeCurrent/2) || totalSolutionSizeCurrent <= 0){
            return SolutionUtils.getBiDirectionCandidateSourceSolution(srcGenotype, queryGenotype, switchCostScore+1, srcIndiList,
                    taxaSourceMap, maxSolutionCount);
        }

        if (totalSolutionSizeCurrent > maxSolutionCount){
            return SolutionUtils.getBiDirectionCandidateSourceSolution(srcGenotype, queryGenotype, switchCostScore+1, srcIndiList,
                    taxaSourceMap, maxSolutionCount);
//            return new EnumMap<>(Solution.Direction.class);
        }

        IntList[] reverseCandidateSolution =
                SolutionUtils.getCandidateSolution(SolutionUtils.reverseSrcGenotype(srcGenotype),
                        SolutionUtils.reverseGenotype(queryGenotype),switchCostScore, srcIndiList, taxaSourceMap);
        EnumMap<Direction, IntList[]> enumMap = new EnumMap<>(Direction.class);
        enumMap.put(Direction.F, forwardCandidateSolutionCurrent);
        enumMap.put(Direction.R, reverseCandidateSolution);
        return enumMap;
    }

    public static int[] getTargetSourceCumLen(List<IntList> solutions){
        int[] cumLen = new int[solutions.size()];
        IntList singleSourceFeatureList = SourceType.getSingleSourceFeatureList();
        for (int i = 0; i < solutions.size(); i++) {
            for (int j = 0; j < solutions.get(i).size(); j=j+3) {
                if (singleSourceFeatureList.contains(solutions.get(i).getInt(j))){
                    cumLen[i]+=(solutions.get(i).getInt(j+1)) - (solutions.get(i).getInt(j+2));
                }
            }
        }
        return cumLen;
    }

    public static IntList coalescentForward(IntList[] solutions){
        int[] miniSolutionSizeArray = SolutionUtils.getOptimalSolutionsSize(solutions);
        int miniSolutionSize = SolutionUtils.getMiniOptimalSolutionSize(solutions);
        int miniSolutionEleCount = SolutionUtils.getMiniSolutionEleCount(solutions);
        int solutionEleCount;
        Set<IntList> solutionSet = new HashSet<>();
        for (int i = 0; i < solutions.length; i++) {
            if (miniSolutionSizeArray[i]!=miniSolutionSize) continue;
            if (solutions[i].size()!=miniSolutionEleCount) continue;
            solutionEleCount =  solutions[i].size();
            // filter Source is NONE
            if (solutionEleCount==3 && (solutions[i].getInt(0)==16)) continue;
            solutionSet.add(solutions[i]);
        }
        List<IntList> solutionList = new ArrayList<>(solutionSet);
        if (solutionList.size()==0) return new IntArrayList();
        int[] targetSourceCumLen = SolutionUtils.getTargetSourceCumLen(solutionList);
        int miniTargetSourceCumLen = Integer.MAX_VALUE;
        for (int i = 0; i < targetSourceCumLen.length; i++) {
            miniTargetSourceCumLen = targetSourceCumLen[i] < miniTargetSourceCumLen ? targetSourceCumLen[i] :
                    miniTargetSourceCumLen;
        }
        IntList miniTargetSourceCumLenSolution=null;
        for (int i = 0; i < targetSourceCumLen.length; i++) {
            if (targetSourceCumLen[i] == miniTargetSourceCumLen){
                miniTargetSourceCumLenSolution = solutionList.get(i);
                break;
            }
        }

        IntList solutionRes = new IntArrayList();
        for (int i = miniTargetSourceCumLenSolution.size()-1; i > 0; i=i-3) {
            solutionRes.add(miniTargetSourceCumLenSolution.getInt(i-2));
            solutionRes.add(miniTargetSourceCumLenSolution.getInt(i));
            solutionRes.add(miniTargetSourceCumLenSolution.getInt(i-1));
        }
        return solutionRes;
    }

    public static IntList coalescentReverse(IntList[] solutions){
        int[] miniSolutionSizeArray = SolutionUtils.getOptimalSolutionsSize(solutions);
        int miniSolutionSize = SolutionUtils.getMiniOptimalSolutionSize(solutions);
        int miniSolutionEleCount = SolutionUtils.getMiniSolutionEleCount(solutions);
        int solutionEleCount;
        Set<IntList> solutionSet = new HashSet<>();
        for (int i = 0; i < solutions.length; i++) {
            if (miniSolutionSizeArray[i]!=miniSolutionSize) continue;
            if (solutions[i].size()!=miniSolutionEleCount) continue;
            solutionEleCount =  solutions[i].size();
            // filter Source is NONE
            if (solutionEleCount==3 && (solutions[i].getInt(0)==16)) continue;
            solutionSet.add(solutions[i]);
        }
        List<IntList> solutionList = new ArrayList<>(solutionSet);
        if (solutionList.size()==0) return new IntArrayList();
        int[] targetSourceCumLen = SolutionUtils.getTargetSourceCumLen(solutionList);
        int miniTargetSourceCumLen = Integer.MAX_VALUE;
        for (int i = 0; i < targetSourceCumLen.length; i++) {
            miniTargetSourceCumLen = targetSourceCumLen[i] < miniTargetSourceCumLen ? targetSourceCumLen[i] :
                    miniTargetSourceCumLen;
        }
        IntList miniTargetSourceCumLenSolution=null;
        for (int i = 0; i < targetSourceCumLen.length; i++) {
            if (targetSourceCumLen[i] == miniTargetSourceCumLen){
                miniTargetSourceCumLenSolution = solutionList.get(i);
                break;
            }
        }
        IntList solutionRes = new IntArrayList();
        int seqLen = solutionList.get(0).getInt(1);
        for (int i = 0; i < miniTargetSourceCumLenSolution.size(); i=i+3) {
            solutionRes.add(miniTargetSourceCumLenSolution.getInt(i));
            solutionRes.add(seqLen-miniTargetSourceCumLenSolution.getInt(i+1));
            solutionRes.add(seqLen-miniTargetSourceCumLenSolution.getInt(i+2));
        }
        return solutionRes;
    }

    /**
     *
     * @param candidateSolutions
     * @return final solution, every three numbers as a group, representing a tract
     * the first number is source population index, equal WindowSource.Source.index()
     * the second and third number is start(inclusive) position and end(inclusive) position
     */
    public static IntList calculateBreakPoint(EnumMap<Direction, IntList[]> candidateSolutions){
        if (candidateSolutions.size()==0) return new IntArrayList();
        IntList forwardSolution = SolutionUtils.coalescentForward(candidateSolutions.get(Direction.F));
        if (forwardSolution.size()==0) return new IntArrayList();
        IntList singleSourceFeatureList = SourceType.getSingleSourceFeatureList();
        int forwardCount, reverseCount;
        forwardCount = forwardSolution.size()/3;
        int forwardLastSourceFeature = forwardSolution.getInt(forwardSolution.size()-3);
        IntList reverseSolution;
        if (singleSourceFeatureList.contains(forwardLastSourceFeature)){
            reverseSolution = SolutionUtils.coalescentReverse(candidateSolutions.get(Direction.R));
            if (reverseSolution.size()==0) return forwardSolution;
            reverseCount = reverseSolution.size()/3;
            if (forwardCount==1 || reverseCount==1) return forwardSolution;

            // 双向动态规划的末端调整需要满足reverseCount > forwardCount的条件
            if (reverseCount < forwardCount) return reverseSolution;
            if (reverseSolution.getInt(reverseSolution.size()-6)==(forwardLastSourceFeature)){
                if((reverseSolution.getInt(reverseSolution.size()-6)!=(reverseSolution.getInt(reverseSolution.size()-3)))){

                    // make sure final range is right
                    if (reverseSolution.getInt(reverseSolution.size()-4) > forwardSolution.getInt(forwardSolution.size()-2)){
                        forwardSolution.set((forwardCount-1)*3+2, reverseSolution.getInt((reverseCount-2)*3+2));
                        forwardSolution.add(reverseSolution.getInt((reverseCount-1)*3+0));
                        forwardSolution.add(reverseSolution.getInt((reverseCount-1)*3+1));
                        forwardSolution.add(reverseSolution.getInt((reverseCount-1)*3+2));
                        return forwardSolution;
                    }

                }
            }
        }
//        IntList engravedSolution = SolutionUtils.engrave(forwardSolution);
        return forwardSolution;
    }

    public static IntList engrave(IntList solution){
        IntList res = new IntArrayList(solution);
        IntList singleSourceFeatureList = SourceType.getSingleSourceFeatureList();
        int preFeature, currentFeature, nextFeature;
        for (int i = 3; i < res.size()-3; i=i+3) {
            preFeature = res.getInt(i-3);
            currentFeature = res.getInt(i);
            nextFeature = res.getInt(i+3);
            if (preFeature!=nextFeature) continue;
            if (singleSourceFeatureList.contains(preFeature) && (!singleSourceFeatureList.contains(currentFeature))){
                if ((preFeature & currentFeature) == preFeature){
                    res.set(i-1, res.getInt(i+5));
                    res.removeElements(i, i+6);
                    i=i-3;
                }
            }
        }
        return res;
    }


    /**
     *
     * @param srcGenotype the first dim is haplotype, the second dim is SNP position
     * @return 反向序列
     */
    public static double[][] reverseSrcGenotype(double[][] srcGenotype){
        double[][] reverseGenotype = new double[srcGenotype.length][];
        for (int i = 0; i < reverseGenotype.length; i++) {
            reverseGenotype[i] = new double[srcGenotype[0].length];
            Arrays.fill(reverseGenotype[i], -1);
        }

        for (int i = 0; i < srcGenotype.length; i++) {
            for (int j = 0; j < srcGenotype[i].length; j++) {
                reverseGenotype[i][srcGenotype[i].length-1-j]=srcGenotype[i][j];
            }
        }
        return reverseGenotype;
    }

    /**
     *
     * @param genotype
     * @return 反向序列
     */
    public static double[] reverseGenotype(double[] genotype){
        double[] reverseGenotype = new double[genotype.length];
        Arrays.fill(reverseGenotype, -1);
        for (int i = 0; i < genotype.length; i++) {
            reverseGenotype[genotype.length-1-i]=genotype[i];
        }
        return reverseGenotype;
    }


}
