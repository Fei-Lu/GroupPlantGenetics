package daxing.v2.localAncestryInfer;

import it.unimi.dsi.fastutil.ints.*;
import org.omg.CORBA.PUBLIC_MEMBER;

import java.util.*;
import java.util.List;

public class SolutionUtils {

    static int iteration=0;

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

    public static IntList[] getCandidateSolution2(double[][] srcGenotype, double[] queryGenotype, double switchCostScore, List<String> srcIndiList,
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

    /**
     *
     * @param srcGenotype
     * @param queryGenotype
     * @param switchCostScore
     * @param srcIndiList
     * @param taxaSourceMap
     * @return
     */
    public static EnumMap<Solution.Direction, IntList[]> getCandidateSourceSolution2(double[][] srcGenotype,
                                                                                                           double[] queryGenotype,
                                                                                                           double switchCostScore,
                                                                                                           List<String> srcIndiList,
                                                                                                           Map<String, WindowSource.Source> taxaSourceMap,
                                                                                                           int maxSolutionCount){

        IntList[] forwardCandidateSolutionCurrent = SolutionUtils.getCandidateSolution2(srcGenotype, queryGenotype, switchCostScore,
                        srcIndiList, taxaSourceMap);
        IntList[] forwardCandidateSolutionNext =
                SolutionUtils.getCandidateSolution2(srcGenotype, queryGenotype, switchCostScore+1,
                        srcIndiList, taxaSourceMap);

        int totalSolutionSizeCurrent = SolutionUtils.getMiniOptimalSolutionSize(forwardCandidateSolutionCurrent);
        int totalSolutionSizeNext = SolutionUtils.getMiniOptimalSolutionSize(forwardCandidateSolutionNext);

//        StringBuilder log = new StringBuilder();

        //  totalSolutionSizeCurrent < 0 是因为 两个Int相乘的结果大于Int max
        if ((totalSolutionSizeCurrent > maxSolutionCount && totalSolutionSizeNext < totalSolutionSizeCurrent/2) || totalSolutionSizeCurrent <= 0){
//            log.setLength(0);
////            log.append("\n").append("iteration "+iteration).append("\n");
//            log.append("Switch cost score is "+switchCostScore).append("\n");
//            log.append("Total solution size is "+totalSolutionSizeCurrent).append("\n").append("\n");
//            System.out.println(log);
            return SolutionUtils.getCandidateSourceSolution2(srcGenotype, queryGenotype, switchCostScore+1, srcIndiList,
                    taxaSourceMap, maxSolutionCount);
        }

        if (totalSolutionSizeCurrent > maxSolutionCount){
//            log.setLength(0);
////            log.append("\n").append("iteration "+iteration).append("\n");
//            log.append("Switch cost score is "+switchCostScore).append("\n");
//            log.append("Total solution size is "+totalSolutionSizeCurrent).append("\n");
//            log.append("Total solution size greater than maxSolutionCount "+maxSolutionCount).append("\n");
//            log.append("continue iterate").append("\n");
//            System.out.println(log);
//            return SolutionUtils.getCandidateSourceSolution2(srcGenotype, queryGenotype, switchCostScore+1, srcIndiList,
//                    taxaSourceMap, maxSolutionCount);
            return new EnumMap<>(Solution.Direction.class);
        }

//        System.out.println();
////        System.out.println("iteration "+iteration);
//        System.out.println("Switch cost score is "+switchCostScore);
//        System.out.println("Total solution size is "+totalSolutionSizeCurrent);
//        System.out.println();
        IntList[] reverseCandidateSolution =
                SolutionUtils.getCandidateSolution2(SolutionUtils.reverseSrcGenotype(srcGenotype),
                        SolutionUtils.reverseGenotype(queryGenotype),switchCostScore, srcIndiList, taxaSourceMap);
        EnumMap<Solution.Direction, IntList[]> enumMap = new EnumMap<>(Solution.Direction.class);
        enumMap.put(Solution.Direction.F, forwardCandidateSolutionCurrent);
        enumMap.put(Solution.Direction.R, reverseCandidateSolution);
        return enumMap;
    }

    public static IntList coalescentForward(IntList[] solutions){
        int[] miniSolutionSizeArray = SolutionUtils.getOptimalSolutionsSize(solutions);
        int miniSolutionSize = SolutionUtils.getMiniOptimalSolutionSize(solutions);
        int solutionEleCount;
        Set<IntList> solutionSet = new HashSet<>();
        for (int i = 0; i < solutions.length; i++) {
            if (miniSolutionSizeArray[i]!=miniSolutionSize) continue;
            solutionEleCount =  solutions[i].size();
            // filter Source is NONE
            if (solutionEleCount==3 && (solutions[i].getInt(0)==16)) continue;
            solutionSet.add(solutions[i]);
        }
        List<IntList> solutionList = new ArrayList<>(solutionSet);
        if (solutionList.size()==0) return new IntArrayList();
        IntList singleSourceFeatureList = SourceType.getSingleSourceFeatureList();
        int sourceFeature;
        IntList solutionRes = new IntArrayList();
        int singleSourceStart, singleSourceEnd, mutipleSourceStart, mutipleSouceEnd;
        for (int i = solutionList.get(0).size()-1; i > 0 ; i=i-3) {
            singleSourceStart = Integer.MIN_VALUE;
            singleSourceEnd = Integer.MAX_VALUE;
            mutipleSourceStart = Integer.MAX_VALUE;
            mutipleSouceEnd = Integer.MIN_VALUE;
            sourceFeature = solutionList.get(0).getInt(i-2);
            if (singleSourceFeatureList.contains(sourceFeature)){
                for (int j = 0; j < solutionList.size(); j++) {
                    singleSourceStart = solutionList.get(j).getInt(i) > singleSourceStart ? solutionList.get(j).getInt(i): singleSourceStart;
                    singleSourceEnd = solutionList.get(j).getInt(i-1) < singleSourceEnd ? solutionList.get(j).getInt(i-1) : singleSourceEnd;
                }
                solutionRes.add(sourceFeature);
                solutionRes.add(singleSourceStart);
                solutionRes.add(singleSourceEnd);
            }else {
                for (int j = 0; j < solutionList.size(); j++) {
                    mutipleSourceStart = solutionList.get(j).getInt(i) < mutipleSourceStart ? solutionList.get(j).getInt(i): mutipleSourceStart;
                    mutipleSouceEnd = solutionList.get(j).getInt(i-1) > mutipleSouceEnd ? solutionList.get(j).getInt(i-1) : mutipleSouceEnd;
                }
                solutionRes.add(sourceFeature);
                solutionRes.add(mutipleSourceStart);
                solutionRes.add(mutipleSouceEnd);
            }
        }
        return solutionRes;
    }

    public static IntList coalescentReverse(IntList[] solutions){
        int[] miniSolutionSizeArray = SolutionUtils.getOptimalSolutionsSize(solutions);
        int miniSolutionSize = SolutionUtils.getMiniOptimalSolutionSize(solutions);
        int solutionEleCount;
        Set<IntList> solutionSet = new HashSet<>();
        for (int i = 0; i < solutions.length; i++) {
            if (miniSolutionSizeArray[i]!=miniSolutionSize) continue;
            solutionEleCount =  solutions[i].size();
            // filter Source is NONE
            if (solutionEleCount==3 && (solutions[i].getInt(0)==16)) continue;
            solutionSet.add(solutions[i]);
        }
        List<IntList> solutionList = new ArrayList<>(solutionSet);
        if (solutionList.size()==0) return new IntArrayList();
        IntList singleSourceFeatureList = SourceType.getSingleSourceFeatureList();
        int sourceFeature;
        IntList solutionRes = new IntArrayList();
        int singleSourceStart, singleSourceEnd, mutipleSourceStart, mutipleSouceEnd;
        int seqLen = solutionList.get(0).getInt(1);
        for (int i = 0; i < solutionList.get(0).size(); i=i+3) {
            singleSourceStart = Integer.MIN_VALUE;
            singleSourceEnd = Integer.MAX_VALUE;
            mutipleSourceStart = Integer.MAX_VALUE;
            mutipleSouceEnd = Integer.MIN_VALUE;
            sourceFeature = solutionList.get(0).getInt(i);
            if (singleSourceFeatureList.contains(sourceFeature)){
                for (int j = 0; j < solutionList.size(); j++) {
                    singleSourceStart = seqLen-solutionList.get(j).getInt(i+1) > singleSourceStart ?
                            seqLen-solutionList.get(j).getInt(i+1): singleSourceStart;
                    singleSourceEnd = seqLen-solutionList.get(j).getInt(i+2) < singleSourceEnd ?
                            seqLen-solutionList.get(j).getInt(i+2) : singleSourceEnd;
                }
                solutionRes.add(sourceFeature);
                solutionRes.add(singleSourceStart);
                solutionRes.add(singleSourceEnd);
            }else {
                for (int j = 0; j < solutionList.size(); j++) {
                    mutipleSourceStart = seqLen-solutionList.get(j).getInt(i+1) < mutipleSourceStart ?
                            seqLen-solutionList.get(j).getInt(i+1): mutipleSourceStart;
                    mutipleSouceEnd = seqLen-solutionList.get(j).getInt(i+2) > mutipleSouceEnd ?
                            seqLen-solutionList.get(j).getInt(i+2) : mutipleSouceEnd;
                }
                solutionRes.add(sourceFeature);
                solutionRes.add(mutipleSourceStart);
                solutionRes.add(mutipleSouceEnd);
            }
        }
        return solutionRes;
    }

    public static IntList calculateBreakPoint(EnumMap<Solution.Direction, IntList[]> candidateSolutions){
        if (candidateSolutions.size()==0) return new IntArrayList();
        IntList forwardSolution = SolutionUtils.coalescentForward(candidateSolutions.get(Solution.Direction.F));
        if (forwardSolution.size()==0) return new IntArrayList();
        IntList singleSourceFeatureList = SourceType.getSingleSourceFeatureList();
        int forwardCount, reverseCount;
        forwardCount = forwardSolution.size()/3;
        int forwardLastSourceFeature = forwardSolution.getInt(forwardSolution.size()-3);
        IntList reverseSolution;
        if (singleSourceFeatureList.contains(forwardLastSourceFeature)){
            reverseSolution = SolutionUtils.coalescentReverse(candidateSolutions.get(Solution.Direction.R));
            if (reverseSolution.size()==0) return forwardSolution;
            reverseCount = reverseSolution.size()/3;
            if (forwardCount==1 || reverseCount==1) return forwardSolution;
            if (reverseSolution.getInt(reverseSolution.size()-6)==(forwardLastSourceFeature)){
                if((reverseSolution.getInt(reverseSolution.size()-6)!=(reverseSolution.getInt(reverseSolution.size()-3)))){
                    forwardSolution.set((forwardCount-1)*3+2, reverseSolution.getInt((reverseCount-2)*3+2));
                    forwardSolution.add(reverseSolution.getInt((reverseCount-1)*3+0));
                    forwardSolution.add(reverseSolution.getInt((reverseCount-1)*3+1));
                    forwardSolution.add(reverseSolution.getInt((reverseCount-1)*3+2));
                    return forwardSolution;
                }
            }
        }
        IntList engravedSolution = SolutionUtils.engrave(forwardSolution);
        return engravedSolution;
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
            if (singleSourceFeatureList.contains(preFeature)){
                if (currentFeature==preFeature+16){
                    res.set(i-1, res.getInt(i+5));
                    res.removeElements(i, i+6);
                    i=i-3;
                }
            }
        }
        return res;
    }

    public static EnumSet<WindowSource.Source>[][] getCandidateSourceSolutionEnumSet(double[][] srcGenotype, double[] queryGenotype,
                                                                                     double switchCostScore,
                                                                                     List<String> srcIndiList,
                                                                                     Map<String, WindowSource.Source> taxaSourceMap){
        double[][] miniCost = SolutionUtils.getMiniCostScore(srcGenotype, queryGenotype, switchCostScore);
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
     * @param srcGenotype
     * @param queryGenotype
     * @param switchCostScore
     * @param srcIndiList
     * @param taxaSourceMap
     * @return
     */
    public static EnumMap<Solution.Direction, EnumSet<WindowSource.Source>[][]> getCandidateSourceSolution(double[][] srcGenotype,
                                                                                double[] queryGenotype,
                                                                                double switchCostScore,
                                                                                List<String> srcIndiList,
                                                                                Map<String, WindowSource.Source> taxaSourceMap,
                                                                                                           int maxSolutionCount){

        EnumSet<WindowSource.Source>[][] forwardCandidateSolutionCurrent =
                SolutionUtils.getCandidateSourceSolutionEnumSet(srcGenotype, queryGenotype, switchCostScore,
                        srcIndiList, taxaSourceMap);
        EnumSet<WindowSource.Source>[][] forwardCandidateSolutionNext =
                SolutionUtils.getCandidateSourceSolutionEnumSet(srcGenotype, queryGenotype, switchCostScore+1,
                        srcIndiList, taxaSourceMap);

        int totalSolutionSizeCurrent = SolutionUtils.getTotalOptimalSolutionSize(forwardCandidateSolutionCurrent);
        int totalSolutionSizeNext = SolutionUtils.getTotalOptimalSolutionSize(forwardCandidateSolutionNext);

        StringBuilder log = new StringBuilder();
        //  totalSolutionSizeCurrent < 0 是因为 两个Int相乘的结果大于Int max
        if ((totalSolutionSizeCurrent > maxSolutionCount && totalSolutionSizeNext < totalSolutionSizeCurrent/2) || totalSolutionSizeCurrent <= 0){
            log.setLength(0);
//            log.append("\n").append("iteration "+iteration).append("\n");
            log.append("Switch cost score is "+switchCostScore).append("\n");
            log.append("Total solution size is "+totalSolutionSizeCurrent).append("\n").append("\n");
            System.out.println(log);
            return SolutionUtils.getCandidateSourceSolution(srcGenotype, queryGenotype, switchCostScore+1, srcIndiList,
                    taxaSourceMap, maxSolutionCount);
        }

        if (totalSolutionSizeCurrent > maxSolutionCount){
            log.setLength(0);
//            log.append("\n").append("iteration "+iteration).append("\n");
            log.append("Switch cost score is "+switchCostScore).append("\n");
            log.append("Total solution size is "+totalSolutionSizeCurrent).append("\n");
            log.append("Total solution size greater than maxSolutionCount "+maxSolutionCount).append("\n");
            log.append("continue iterate").append("\n");
            System.out.println(log);
            return SolutionUtils.getCandidateSourceSolution(srcGenotype, queryGenotype, switchCostScore+1, srcIndiList,
                    taxaSourceMap, maxSolutionCount);
        }


        System.out.println();
//        System.out.println("iteration "+iteration);
        System.out.println("Switch cost score is "+switchCostScore);
        System.out.println("Total solution size is "+totalSolutionSizeCurrent);
        System.out.println();
//        iteration =0;
        EnumSet<WindowSource.Source>[][] reverseCandidateSolution =
                SolutionUtils.getCandidateSourceSolutionEnumSet(SolutionUtils.reverseSrcGenotype(srcGenotype),
                        SolutionUtils.reverseGenotype(queryGenotype),switchCostScore, srcIndiList, taxaSourceMap);
        EnumMap<Solution.Direction, EnumSet<WindowSource.Source>[][]> enumMap = new EnumMap<>(Solution.Direction.class);
        enumMap.put(Solution.Direction.F, forwardCandidateSolutionCurrent);
        enumMap.put(Solution.Direction.R, reverseCandidateSolution);
        return enumMap;
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

    /**
     *
     * @param candidateSourceSolution the first dim is SNP position
     * @param count the count of mini cost score index, to prevent the sum of these paths from being greater than Integer.MAX_VALUE
     * @return the total number of optimal solutions corresponding to a mini cost score
     */
    public static int getOptimalSolutionsSize(EnumSet<WindowSource.Source>[] candidateSourceSolution, int count){
        int size=candidateSourceSolution[0].size();
        for (int i = 1; i < candidateSourceSolution.length; i++) {
            if (candidateSourceSolution[i].size() == 1) continue;
            if (candidateSourceSolution[i-1].equals(candidateSourceSolution[i])) continue;
            size *=candidateSourceSolution[i].size();
            if (size > ((Integer.MAX_VALUE)/count)){
                return (Integer.MAX_VALUE)/count;
            }
        }
        return size;
    }

    /**
     *
     * @param candidateSourceSolutions the first dim is mini cost score index, the second dim is SNP position
     * @return the total number of optimal candidateSourceSolutions corresponding to all mini cost scores
     */
    public static int[] getOptimalSolutionsSize(EnumSet<WindowSource.Source>[][] candidateSourceSolutions){

        int[] size = new int[candidateSourceSolutions.length];
        for (int i = 0; i < candidateSourceSolutions.length; i++) {
            size[i]= SolutionUtils.getOptimalSolutionsSize(candidateSourceSolutions[i], candidateSourceSolutions.length);
        }
        return size;
    }

    public static int getTotalOptimalSolutionSize(EnumSet<WindowSource.Source>[][] candidateSourceSolutions){
        int[] size = SolutionUtils.getOptimalSolutionsSize(candidateSourceSolutions);
        int sum = 0;
        for (int i = 0; i < size.length; i++) {
            sum+=size[i];
        }
        return sum;
    }


}
