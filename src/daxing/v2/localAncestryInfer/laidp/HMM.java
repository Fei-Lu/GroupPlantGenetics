package daxing.v2.localAncestryInfer.laidp;

import it.unimi.dsi.fastutil.ints.*;
import java.util.Arrays;
import java.util.BitSet;
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
            dp[i][0] = Math.log(start_prob[i]) +  (obs[0]==1 ? Math.log(alts[i][0]) : Math.log(1 - alts[i][0])) ;
        }

        for (int variantIndex = 1; variantIndex < variantNum; variantIndex++) {
            for (int currentState = 0; currentState < stateNum; currentState++) {
                double maxScore = Double.NEGATIVE_INFINITY;
                int prevState = 0;

                for (int fromState = 0; fromState < stateNum; fromState++) {
                    double score = dp[fromState][variantIndex - 1] + Math.log(states_trans_prob[fromState][currentState]) +
                            (obs[variantIndex] ==1 ? Math.log(alts[currentState][variantIndex]) :
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

        bestState = bestState < 0 ? 0 : bestState;
        int[] path = new int[variantNum];
        path[variantNum-1] = bestState;

        for (int i = variantNum - 2; i >=0; i--) {
            path[i] = backPath[path[i+1]][i];
        }

        return path;
    }

    /**
     *
     * @param alts alternative allele frequency of different source population, order seen Source class
     * @param states_trans_prob dim1 and dim2 are different source population (order seen Source)
     * @param start_prob dim1 is different source population (order seen Source)
     * @param obs genotype of an admixed haplotype, 0 means reference allele, 1 means alternative allele
     * @return forward matrix, dim1 is different source population, dim2 is different variant
     */
    public static double[][] getForward(double[][] alts, double[][] states_trans_prob, double[] start_prob, int[] obs){

        int stateNum = alts.length;
        int variantNum = alts[0].length;

        double[][] forward = new double[stateNum][variantNum];
        for (double[] doubles : forward) {
            Arrays.fill(doubles, -1);
        }

        // Initializing the Forward Matrix
        for (int i = 0; i < stateNum; i++) {
            forward[i][0] = start_prob[i] * (obs[0]==1 ? (alts[i][0]) : (1 - alts[i][0]));
        }

        for (int variantIndex = 1; variantIndex < variantNum; variantIndex++) {
            for (int currentStateIndex = 0; currentStateIndex < stateNum; currentStateIndex++) {
                forward[currentStateIndex][variantIndex] = 0;
                for (int fromStateIndex = 0; fromStateIndex < stateNum; fromStateIndex++) {
                    forward[currentStateIndex][variantIndex] += forward[fromStateIndex][variantIndex-1] * states_trans_prob[fromStateIndex][currentStateIndex];
                }
                forward[currentStateIndex][variantIndex] *= (obs[variantIndex]==1 ?
                        (alts[currentStateIndex][variantIndex]) :
                        (1 - alts[currentStateIndex][variantIndex]));
            }
        }

        return forward;
    }

    /**
     *
     * @param alts alternative allele frequency of different source population, order seen Source class
     * @param states_trans_prob  dim1 and dim2 are different source population (order seen Source)
     * @param obs genotype of an admixed haplotype, 0 means reference allele, 1 means alternative allele
     * @return backward matrix, dim1 is different source population, dim2 is different variant
     */
    public static double[][] getBackward(double[][] alts, double[][] states_trans_prob, int[] obs){
        int stateNum = alts.length;
        int variantNum = alts[0].length;

        double[][] backward = new double[stateNum][variantNum];
        for (double[] doubles : backward) {
            Arrays.fill(doubles, -1);
        }

        // Initializing the Backward Matrix
        for(int i = 0; i < stateNum; i++){
            backward[i][variantNum -1] = 1;
        }

        for (int variantIndex = variantNum - 2; variantIndex >= 0; variantIndex--) {
            for (int currentStateIndex = 0; currentStateIndex < stateNum; currentStateIndex++) {
                backward[currentStateIndex][variantIndex] = 0;
                for (int fromStateIndex = 0; fromStateIndex < stateNum; fromStateIndex++) {
                    backward[currentStateIndex][variantIndex] += backward[fromStateIndex][variantIndex+1] * states_trans_prob[fromStateIndex][currentStateIndex];
                }
                backward[currentStateIndex][variantIndex] *= (obs[variantIndex+1]==1 ?
                        (alts[currentStateIndex][variantIndex+1]) :
                        (1 - alts[currentStateIndex][variantIndex+1]));
            }
        }
        return backward;
    }

    /**
     * @param forwardMatrix forward matrix
     * @param variantNum variant number used in forward matrix
     * @return likelihood probability of the input admixed haplotype
     */
    public static double getForwardLogLikelihood(double[][] forwardMatrix, int variantNum){
        double likelihood = 0;
        for (double[] matrix : forwardMatrix) {
            likelihood += matrix[variantNum - 1];
        }
        return Math.log(likelihood);
    }

    /**
     *
     * @param backwardMatrix backward matrix
     * @param alts alternative allele frequency of different source population, order seen Source class
     * @param start_prob dim1 is different source population (order seen Source)
     * @param obs genotype of an admixed haplotype, 0 means reference allele, 1 means alternative allele
     * @return likelihood probability of the input admixed haplotype
     */
    public static double getBackwardLogLikelihood(double[][] backwardMatrix, double[][] alts,
                                                  double[] start_prob, int[] obs){
        int stateNum = backwardMatrix.length;
        double likelihood = 0;
        for (int i = 0; i < stateNum; i++) {
            likelihood += backwardMatrix[i][0] * (obs[0] == 1 ? (alts[i][0]) : (1 - alts[i][0])) * start_prob[i];
        }
        return Math.log(likelihood);
    }

    /**
     *
     * @param alts alternative allele frequency of different source population, order seen Source class
     * @param states_trans_prob dim1 and dim2 are different source population (order seen Source)
     * @param obs genotype of an admixed haplotype, 0 means reference allele, 1 means alternative allele
     * @param forward forward matrix
     * @param backward backward matrix
     * @return states_trans_prob after one iteration of the EM (expectation-maximization) algorithm,
     * dim1 and dim2 are different source population (order seen Source)
     */
    private static double[][] em_step(double[][] alts, double[][] states_trans_prob, int[] obs,
                                       double[][] forward, double[][] backward){
        int stateNum = forward.length;
        int variantNum = alts[0].length;

        // E-step
        double[] denominator_notQuite = new double[variantNum];
        double[][][] notQuiteXi = new double[variantNum-1][stateNum][stateNum];
        double[][][] xi = new double[variantNum-1][stateNum][stateNum];

        Arrays.fill(denominator_notQuite, -1);

        for (int variantIndex = 0; variantIndex < variantNum; variantIndex++) {
            denominator_notQuite[variantIndex] = 0;
            for (int stateIndex = 0; stateIndex < stateNum; stateIndex++) {
                denominator_notQuite[variantIndex]+= forward[stateIndex][variantIndex] * backward[stateIndex][variantIndex];
            }
        }

        for (int variantIndex = 0; variantIndex < variantNum-1; variantIndex++) {
            for (int fromStateIndex = 0; fromStateIndex < stateNum; fromStateIndex++) {
                for (int toStateIndex = 0; toStateIndex < stateNum; toStateIndex++) {
                    notQuiteXi[variantIndex][fromStateIndex][toStateIndex]= forward[fromStateIndex][variantIndex] *
                            states_trans_prob[fromStateIndex][toStateIndex] *
                            (obs[variantIndex+1]==1 ? (alts[toStateIndex][variantIndex+1]) :
                            (1 - alts[toStateIndex][variantIndex+1])) *
                            backward[toStateIndex][variantIndex+1];
                    xi[variantIndex][fromStateIndex][toStateIndex]= notQuiteXi[variantIndex][fromStateIndex][toStateIndex]/denominator_notQuite[variantIndex];
                }
            }
        }

        // M-step
        double[] denominator_xi = new double[stateNum];
        double[][] trans_prob_em = new double[stateNum][stateNum];

        Arrays.fill(denominator_xi, -1);

        for (int fromStateIndex = 0; fromStateIndex < stateNum; fromStateIndex++) {
            denominator_xi[fromStateIndex] = 0;
            for (int variantIndex = 0; variantIndex < variantNum-1; variantIndex++) {
                for (int toStateIndex = 0; toStateIndex < stateNum; toStateIndex++) {
                    denominator_xi[fromStateIndex]+=xi[variantIndex][fromStateIndex][toStateIndex];
                }
            }
        }


        for (int fromStateIndex = 0; fromStateIndex < stateNum; fromStateIndex++) {
            for (int toStateIndex = 0; toStateIndex <stateNum; toStateIndex++) {
                double trans_prob_em_numerator = 0;
                for (int variantIndex = 0; variantIndex < variantNum - 1; variantIndex++) {
                    trans_prob_em_numerator+=xi[variantIndex][fromStateIndex][toStateIndex];
                }
                trans_prob_em[fromStateIndex][toStateIndex]=trans_prob_em_numerator/denominator_xi[fromStateIndex];
            }
        }

        return trans_prob_em;
    }

    /**
     *
     * @param lastForward last forward matrix
     * @param currentForward current forward matrix
     * @param variantNum variant number
     * @param logThreshold Threshold of log-likelihood
     * @return true if and only if the difference between the current log-likelihood value and the previous
     * one is less than the threshold
     */
    public static boolean ifConvergence(double[][] lastForward, double[][] currentForward, int variantNum,
                                    double logThreshold){
        double lastLikelihood = HMM.getForwardLogLikelihood(lastForward, variantNum);
        double currentLikelihood = HMM.getForwardLogLikelihood(currentForward, variantNum);
        double delta = currentLikelihood - lastLikelihood;
        return Math.abs(delta) < logThreshold;
    }

    /**
     *
     * @param alts alternative allele frequency of different source population, order seen Source class
     * @param states_trans_prob dim1 and dim2 are different source population (order seen Source)
     * @param start_prob dim1 is different source population (order seen Source)
     * @param obs genotype of an admixed haplotype, 0 means reference allele, 1 means alternative allele
     * @param maxEMStep max em-step number
     * @param logThreshold Threshold of log-likelihood
     * @return states_trans_prob estimated by forward-backward algorithm
     * dim1 and dim2 are different source population (order seen Source)
     */
    public static double[][] forwardBackward(double[][] alts, double[][] states_trans_prob,
                                             double[] start_prob, int[] obs, int maxEMStep, double logThreshold){

        double[][] last_trans_prob = states_trans_prob;
        double[][] lastForward = HMM.getForward(alts, states_trans_prob, start_prob, obs);
        double[][] lastBackward= HMM.getBackward(alts, states_trans_prob, obs);

        double[][] current_trans_prob;
        double[][] currentForward;
        for (int i = 0; i < maxEMStep; i++) {
            current_trans_prob = HMM.em_step(alts, last_trans_prob, obs, lastForward, lastBackward);
            currentForward = HMM.getForward(alts, current_trans_prob, start_prob, obs);
            if(HMM.ifConvergence(lastForward, currentForward, alts[0].length, logThreshold)){
                return current_trans_prob;
            }else {
                last_trans_prob = current_trans_prob;
                lastForward = currentForward;
                lastBackward = HMM.getBackward(alts, current_trans_prob, obs);
            }
        }
        return last_trans_prob;

    }

    /**
     *
     * @param alts alternative allele frequency of different source population, order seen Source class
     * @param states_trans_prob dim1 and dim2 are different source population (order seen Source)
     * @param start_prob dim1 is different source population (order seen Source)
     * @param obs genotype of an admixed haplotype, 0 means reference allele, 1 means alternative allele
     * @param emStepNum em-step number
     * @return states_trans_prob estimated by forward-backward algorithm
     * dim1 and dim2 are different source population (order seen Source)
     */
    public static double[][] forwardBackward(double[][] alts, double[][] states_trans_prob,
                                             double[] start_prob, int[] obs, int emStepNum){
        double[][] last_trans_prob = states_trans_prob;
        double[][] lastForward = HMM.getForward(alts, states_trans_prob, start_prob, obs);
        double[][] lastBackward= HMM.getBackward(alts, states_trans_prob, obs);

        double[][] current_trans_prob;
        double[][] currentForward;
        for (int i = 0; i < emStepNum; i++) {
            current_trans_prob = HMM.em_step(alts, last_trans_prob, obs, lastForward, lastBackward);
            currentForward = HMM.getForward(alts, current_trans_prob, start_prob, obs);

            last_trans_prob = current_trans_prob;
            lastForward = currentForward;
            lastBackward = HMM.getBackward(alts, current_trans_prob, obs);
        }
        return last_trans_prob;
    }

    public static int[] viterbi(double[][] alts, double[][] states_trans_prob, BitSet obs){
        int stateNum = alts.length;
        int variantNum = alts[0].length;

        double[][] dp = new double[stateNum][variantNum];
        int[][] backPath = new int[stateNum][variantNum];

        for (int i = 0; i < stateNum; i++) {
            dp[i][0] = obs.get(0) ? alts[i][0] : (1 - alts[i][0]);
        }

        for (int variantIndex = 1; variantIndex < variantNum; variantIndex++) {
            for (int currentState = 0; currentState < stateNum; currentState++) {
                double maxScore = Double.NEGATIVE_INFINITY;
                int prevState = 0;
                for (int fromState = 0; fromState < stateNum; fromState++) {
                    double score = dp[fromState][variantIndex - 1] * states_trans_prob[fromState][currentState] +
                            (obs.get(variantIndex) ? (alts[currentState][variantIndex]) :
                                    (1 - alts[currentState][variantIndex]));
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


