package daxing.v2.localAncestryInfer;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;

import java.util.*;

public class LocalAncestryInference {

    public static List<TIntList> getMiniPath(double[][] srcGenotype, double[] queryGenotype){
        double switchCostScore= 1.5;
//        int[][] srcGenotype = {{0,1,0,1,0,1,0,0,0,0,1,1},
//                            {0,0,0,1,0,1,1,0,0,0,1,1},
//                            {0,0,1,0,1,0,0,0,1,0,1,1},
//                            {0,0,0,0,1,0,1,0,1,1,1,1},
//                            {1,1,0,0,0,0,1,1,1,1,0,0},
//                            {1,0,0,1,0,0,1,1,1,1,0,0}};
//        int[] queryGenotype =       {1,1,0,0,0,1,0,0,1,1,1,1};


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


        // initialize solution of path
        TIntList initializeSolution;
        List<TIntList>[] solutionCurrent = new List[srcGenotype.length];
        List<TIntList>[] solutionPast = new List[srcGenotype.length];
        for (int i = 0; i < solutionCurrent.length; i++) {
            solutionCurrent[i] = new ArrayList<>();
            initializeSolution = new TIntArrayList();
            initializeSolution.add(i);
            solutionCurrent[i].add(initializeSolution);
        }

        for (int i = 0; i < solutionPast.length; i++) {
            solutionPast[i] = new ArrayList<>();
            initializeSolution = new TIntArrayList();
            initializeSolution.add(i);
            solutionPast[i].add(initializeSolution);
        }


        // i is SNP position
        // j is haplotype index of source population
        TIntList switchIndexList;
        TIntList switchHaplotype;
        for (int i = 1; i < distance[0].length; i++) {

            // j-1 SNP位置，单倍型路径发生switch对应的最小Cost
            double miniCostSwitch=Double.MAX_VALUE;
            switchIndexList = new TIntArrayList();
            for (int j = 0; j < distance.length; j++) {
                miniCostSwitch = miniCost[j][i-1] < miniCostSwitch ? miniCost[j][i-1] : miniCostSwitch;
            }
            for (int j = 0; j < distance.length; j++) {
                if (miniCost[j][i-1] == miniCostSwitch){
                    switchIndexList.add(j);
                }
            }

            for (int j = 0; j < distance.length; j++) {
                // 最小cost路径对应当前haplotype
                if (miniCost[j][i-1] < miniCostSwitch+switchCostScore){
                    miniCost[j][i] = miniCost[j][i-1] + distance[j][i];
                    for (int k = 0; k < solutionCurrent[j].size(); k++) {
                        solutionCurrent[j].get(k).add(j);
                    }
                }else {
                    // 最小cost路径对应转换单倍型
                    miniCost[j][i] = miniCostSwitch+switchCostScore+distance[j][i];
                    solutionCurrent[j] = new ArrayList<>();
                    for (int k = 0; k < switchIndexList.size(); k++) {
                        if (switchIndexList.get(k)==j) continue;
                        for (int l = 0; l < solutionPast[switchIndexList.get(k)].size(); l++) {
                            switchHaplotype = new TIntArrayList();
                            switchHaplotype.addAll(solutionPast[switchIndexList.get(k)].get(l));
                            solutionCurrent[j].add(switchHaplotype);
                        }
                    }

                    // 添加当前单倍型所在的索引
                    for (int k = 0; k < solutionCurrent[j].size(); k++) {
                        solutionCurrent[j].get(k).add(j);
                    }

                    // 如果转换单倍型和当前单倍型最小Cost相同，则两种路径同时保留
                    if (miniCost[j][i-1] == miniCostSwitch+switchCostScore){
                        for (int k = 0; k < solutionPast[j].size(); k++) {
                            switchHaplotype = new TIntArrayList();
                            switchHaplotype.addAll(solutionPast[j].get(k));
                            switchHaplotype.add(j);
                            solutionCurrent[j].add(switchHaplotype);
                        }
                    }
                }
            }

            // assign solutionPast to solutionCurrent
            for (int j = 0; j < solutionCurrent.length; j++) {
                solutionPast[j] = new ArrayList<>();
                for (int k = 0; k < solutionCurrent[j].size(); k++) {
                    solutionPast[j].add(new TIntArrayList());
                }
                for (int k = 0; k < solutionCurrent[j].size(); k++) {
                    solutionPast[j].get(k).addAll(solutionCurrent[j].get(k));
                }
            }
        }

        // miniCost index
        double mini = Double.MAX_VALUE;
        for (int i = 0; i < miniCost.length; i++) {
            mini = miniCost[i][miniCost[0].length-1] < mini ? miniCost[i][miniCost[0].length-1]: mini;
        }
        TIntList solutionIndexList = new TIntArrayList();
        for (int i = 0; i < miniCost.length; i++) {
            if (miniCost[i][miniCost[0].length-1]==mini){
                solutionIndexList.add(i);
            }
        }

        // optium solution
        List<TIntList> optium=new ArrayList<>();
        for (int i = 0; i < solutionIndexList.size(); i++) {
            for (int j = 0; j < solutionCurrent[solutionIndexList.get(i)].size(); j++) {
                optium.add(solutionCurrent[solutionIndexList.get(i)].get(j));
            }
        }

        return optium;
    }

    public static List<TIntList> getMiniCostPath(BitSet[] scrGenotype, BitSet queryGenotype, int seqLength, double switchCostScore) {
//        double switchCostScore = 1.5;
        BitSet[] distance = new BitSet[scrGenotype.length];
        for (int i = 0; i < scrGenotype.length; i++) {
            scrGenotype[i].xor(queryGenotype);
            distance[i] = (BitSet) scrGenotype[i].clone();
        }

        // initialize mini cost score
        double[][] miniCost = new double[distance.length][];
        for (int i = 0; i < miniCost.length; i++) {
            miniCost[i] = new double[seqLength];
            Arrays.fill(miniCost[i], -1);
            miniCost[i][0] = distance[i].get(0) == true ? 1 : 0;
        }

        // initialize solution of path
        TIntList initializeSolution;
        List<TIntList>[] solutionCurrent = new List[distance.length];
        List<TIntList>[] solutionPast = new List[distance.length];
        for (int i = 0; i < solutionCurrent.length; i++) {
            solutionCurrent[i] = new ArrayList<>();
            initializeSolution = new TIntArrayList();
            initializeSolution.add(i);
            solutionCurrent[i].add(initializeSolution);
        }

        for (int i = 0; i < solutionPast.length; i++) {
            solutionPast[i] = new ArrayList<>();
            initializeSolution = new TIntArrayList();
            initializeSolution.add(i);
            solutionPast[i].add(initializeSolution);
        }

        // i is SNP position
        // j is haplotype index of source population
        TIntList switchIndexList;
        TIntList switchHaplotype;
        for (int i = 1; i < seqLength; i++) {

            // j-1 SNP位置，单倍型路径发生switch对应的最小Cost
            double miniCostSwitch = Double.MAX_VALUE;
            switchIndexList = new TIntArrayList();
            for (int j = 0; j < distance.length; j++) {
                miniCostSwitch = miniCost[j][i - 1] < miniCostSwitch ? miniCost[j][i - 1] : miniCostSwitch;
            }
            for (int j = 0; j < distance.length; j++) {
                if (miniCost[j][i - 1] == miniCostSwitch) {
                    switchIndexList.add(j);
                }
            }

            for (int j = 0; j < distance.length; j++) {
                // 最小cost路径对应当前haplotype
                if (miniCost[j][i - 1] < miniCostSwitch + switchCostScore) {
                    miniCost[j][i] = miniCost[j][i - 1] + (distance[j].get(i) == true ? 1 : 0);
                    for (int k = 0; k < solutionCurrent[j].size(); k++) {
                        solutionCurrent[j].get(k).add(j);
                    }
                } else {
                    // 最小cost路径对应转换单倍型
                    miniCost[j][i] = miniCostSwitch + switchCostScore + (distance[j].get(i) == true ? 1 : 0);
                    solutionCurrent[j] = new ArrayList<>();
                    for (int k = 0; k < switchIndexList.size(); k++) {
                        if (switchIndexList.get(k) == j) continue;
                        for (int l = 0; l < solutionPast[switchIndexList.get(k)].size(); l++) {
                            switchHaplotype = new TIntArrayList();
                            switchHaplotype.addAll(solutionPast[switchIndexList.get(k)].get(l));
                            solutionCurrent[j].add(switchHaplotype);
                        }
                    }

                    // 添加当前单倍型所在的索引
                    for (int k = 0; k < solutionCurrent[j].size(); k++) {
                        solutionCurrent[j].get(k).add(j);
                    }

                    // 如果转换单倍型和当前单倍型最小Cost相同，则两种路径同时保留
                    if (miniCost[j][i - 1] == miniCostSwitch + switchCostScore) {
                        for (int k = 0; k < solutionPast[j].size(); k++) {
                            switchHaplotype = new TIntArrayList();
                            switchHaplotype.addAll(solutionPast[j].get(k));
                            switchHaplotype.add(j);
                            solutionCurrent[j].add(switchHaplotype);
                        }
                    }
                }
            }

            // assign solutionPast to solutionCurrent
            for (int j = 0; j < solutionCurrent.length; j++) {
                solutionPast[j] = new ArrayList<>();
                for (int k = 0; k < solutionCurrent[j].size(); k++) {
                    solutionPast[j].add(new TIntArrayList());
                }
                for (int k = 0; k < solutionCurrent[j].size(); k++) {
                    solutionPast[j].get(k).addAll(solutionCurrent[j].get(k));
                }
            }

        }

        // miniCost index
        double mini = Double.MAX_VALUE;
        for (int i = 0; i < miniCost.length; i++) {
            mini = miniCost[i][miniCost[0].length-1] < mini ? miniCost[i][miniCost[0].length-1]: mini;
        }
        TIntList solutionIndexList = new TIntArrayList();
        for (int i = 0; i < miniCost.length; i++) {
            if (miniCost[i][miniCost[0].length-1]==mini){
                solutionIndexList.add(i);
            }
        }

        // optium solution
        List<TIntList> optium=new ArrayList<>();
        for (int i = 0; i < solutionIndexList.size(); i++) {
            for (int j = 0; j < solutionCurrent[solutionIndexList.get(i)].size(); j++) {
                optium.add(solutionCurrent[solutionIndexList.get(i)].get(j));
            }
        }

        return optium;
    }
}
