package daxing.v2.localAncestryInfer.laidp;

import it.unimi.dsi.fastutil.ints.*;
import java.util.stream.IntStream;

public class HMM {

    /**
     *
     * @param alts dim1 is different source population (order seen Source), dim2 is different variant
     * @param states_trans_prob dim1 and dim2 are different source population (order seen Source)
     * @param start_prob dim1 is different source population (order seen Source)
     * @return hide state, dim1 is different path, dim2 is different variant, number represent source feature (may
     * contain multiple sources)
     */
    public static int[][] hmm(double[][] alts, double[][] states_trans_prob, double[] start_prob,
                              int[] obs){
        int stateNum = alts.length;
        int variantNum = alts[0].length;

        double[][] dp = new double[stateNum][variantNum];
        IntList[][] backPath = new IntList[stateNum][variantNum];


        for (int i = 0; i < stateNum; i++) {
            dp[i][0] = Math.log(start_prob[i]) +  (obs[0] == 1 ? Math.log(alts[i][0]) : Math.log(1 - alts[i][0]));
            backPath[i][0] = new IntArrayList();
            backPath[i][0].add(i);
        }

        for (int variantIndex = 1; variantIndex < variantNum; variantIndex++) {
            for (int currentState = 0; currentState < stateNum; currentState++) {
                double maxScore = Double.NEGATIVE_INFINITY;
                IntList prevStateList = new IntArrayList();
                for (int fromState = 0; fromState < stateNum; fromState++) {
                    double score = dp[fromState][variantIndex - 1] + Math.log(states_trans_prob[fromState][currentState]) +
                            (obs[variantIndex] == 1 ? Math.log(alts[currentState][variantIndex]) :
                                    Math.log(1 - alts[currentState][variantIndex]));
                    if (score > maxScore){
                        maxScore =  score;
                        prevStateList.clear();
                        prevStateList.add(fromState);
                    }else if (score == maxScore){
                        prevStateList.add(fromState);
                    }
                }

                dp[currentState][variantIndex] = maxScore;
                backPath[currentState][variantIndex] = prevStateList;
            }
        }

        double bestScore = Double.NEGATIVE_INFINITY;
        IntList bestStateList = new IntArrayList();
        for (int stateIndex = 0; stateIndex < stateNum; stateIndex++) {
            if (dp[stateIndex][variantNum-1] > bestScore){
                bestScore = dp[stateIndex][variantNum-1];
                bestStateList.clear();
                bestStateList.add(stateIndex);
            }else if(dp[stateIndex][variantNum-1] == bestScore){
                bestStateList.add(stateIndex);
            }
        }


        int[][] paths = new int[bestStateList.size()][variantNum];

        IntSet currentStateSet;
        IntSet prevStateSet;
        int currentFeature, prevFeature;

        for (int i = 0; i < bestStateList.size(); i++) {
            currentStateSet = new IntArraySet(backPath[bestStateList.getInt(i)][variantNum-1]);
            currentFeature = HMM.transform(currentStateSet);
            int[] path = new int[variantNum];
            path[variantNum-1] = currentFeature;
            for (int j = variantNum -2; j >= 0; j--) {
                prevStateSet = new IntArraySet();
                for (int currentState : currentStateSet){
                    prevStateSet.addAll(backPath[currentState][j]);
                }
                prevFeature = HMM.transform(prevStateSet);
                path[j] = prevFeature;
                currentStateSet = prevStateSet;
            }
            paths[i] = path;
        }
        return paths;
    }

    public static Int2IntMap indices2FeatureMap(){
        int[] indices = IntStream.range(0, 6).toArray();
        int[] future = {1,2,4,8,16,32};
        Int2IntMap int2IntMap = new Int2IntArrayMap();
        for (int i = 0; i < indices.length; i++) {
            int2IntMap.put(indices[i], future[i]);
        }
        return int2IntMap;
    }

    public static int transform(IntSet indexSet){
        IntIterator integerIterator = indexSet.iterator();
        int index = integerIterator.nextInt();
        while (integerIterator.hasNext()){
            index = index | indices2FeatureMap().get(integerIterator.nextInt());
        }
        return index;
    }

    /**
     *
     * @param alts dim1 is different source population (order seen Source), dim2 is different variant
     * @param states_trans_prob dim1 and dim2 are different source population (order seen Source)
     * @param start_prob dim1 is different source population (order seen Source)
     * @return hide state, dim1 is different variant, number represent source index
     */
    public static int[] hmm2(double[][] alts, double[][] states_trans_prob, double[] start_prob, int[] obs){
        int stateNum = alts.length;
        int variantNum = alts[0].length;

        double[][] dp = new double[stateNum][variantNum];
        int[][] backPath = new int[stateNum][variantNum];


        for (int i = 0; i < stateNum; i++) {
            dp[i][0] = Math.log(start_prob[i]) +  (obs[0] == 1 ? Math.log(alts[i][0]) : Math.log(1 - alts[i][0])) ;
        }

        for (int variantIndex = 1; variantIndex < variantNum; variantIndex++) {
            for (int currentState = 0; currentState < stateNum; currentState++) {
                double maxScore = Double.NEGATIVE_INFINITY;
                int prevState = 0;

                for (int fromState = 0; fromState < stateNum; fromState++) {
                    double score = dp[fromState][variantIndex - 1] + Math.log(states_trans_prob[fromState][currentState]) +
                            (obs[variantIndex] == 1 ? Math.log(alts[currentState][variantIndex]) :
                                    Math.log(1 - alts[currentState][variantIndex]));
                    if (score > maxScore){
                        maxScore =  score;
                        prevState = fromState;
                    }
                }

                dp[currentState][variantIndex] = maxScore;
                backPath[currentState][variantIndex-1] = prevState;
            }
        }

        double bestScore = Double.NEGATIVE_INFINITY;
        int bestState = -1;

        for (int stateIndex = 0; stateIndex < stateNum; stateIndex++) {
            if (dp[stateIndex][variantNum-1] > bestScore){
                bestScore = dp[stateIndex][variantNum-1];
                bestState = stateIndex;
            }
        }

        int[] path = new int[variantNum];
        path[variantNum-1] = bestState;

        for (int i = variantNum - 2; i >=0; i--) {
            path[i] = backPath[path[i+1]][i];
        }

        return path;
    }
}


