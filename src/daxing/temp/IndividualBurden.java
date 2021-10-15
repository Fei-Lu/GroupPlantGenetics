package daxing.temp;

import com.google.common.base.Joiner;
import com.google.common.primitives.Ints;
import daxing.common.factors.LoadType;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.wheat.RefV1Utils;
import java.util.Arrays;

public class IndividualBurden {

    //0: 0/0 1: 1/1 2: 0/1 0为ancestral, 1为derived, 2为杂合
    String taxon;
    int chrID;
    TIntArrayList[][][] genotypeCount; // dim1: fdBins, dim2: p3, dim3: loadGroup
    int[][][] burdenCount;         // dim1: fdBins, dim2: p3, dim3: burden 9 column

    public IndividualBurden(String taxon, int chrID, int fdBinsNum){
        TIntArrayList[][][] genotypeCount=new TIntArrayList[fdBinsNum][][];
        int[][][] burdenCount= new int[fdBinsNum][][];
        for (int j = 0; j < genotypeCount.length; j++) {
            genotypeCount[j] = new TIntArrayList[IntrogressionRegion.P3.values().length][];
            burdenCount[j] = new int[IntrogressionRegion.P3.values().length][];
            for (int k = 0; k < genotypeCount[j].length; k++) {
                genotypeCount[j][k]=new TIntArrayList[LoadType.values().length];
                burdenCount[j][k] = new int[9];
                for (int l = 0; l < genotypeCount[j][k].length; l++) {
                    genotypeCount[j][k][l]=new TIntArrayList();
                }
                Arrays.fill(burdenCount[j][k], -1);
            }
        }
        this.chrID= chrID;
        this.genotypeCount=genotypeCount;
        this.taxon=taxon;
        this.burdenCount=burdenCount;
    }

    public void addLoadCount(double fdBin, IntrogressionRegion.P3 p3, byte[] indexGenotype){
        int fdBinIndex= (int)(fdBin*10-1);
        this.genotypeCount[fdBinIndex][p3.index][indexGenotype[0]].add(indexGenotype[1]);
    }

    private void calculate(){
        int[][] synNonDelBurdenCount3; // dim1: loadGroup, dim2: num,numDerived, numHeter
        TIntArrayList burdenCountList; //        int numSyn, numDerivedInSyn, numHeterInSyn, numNonsyn, numDerivedInNonsyn, numHeterInNonsyn, numHGDeleterious
//                , numDerivedInHGDeleterious, numHeterInHGDeleterious;
        for (int i = 0; i < genotypeCount.length; i++) {
            for (int k = 0; k < genotypeCount[i].length; k++) {
                synNonDelBurdenCount3 = new int[LoadType.values().length][];
                for (int l = 0; l < synNonDelBurdenCount3.length; l++) {
                    synNonDelBurdenCount3[l]=new int[3];
                }
                for (int l = 0; l < genotypeCount[i][k].length; l++) {
                    synNonDelBurdenCount3[l]=calculateNumNumDerivedNumHeter(genotypeCount[i][k][l]);
                }
                burdenCountList = new TIntArrayList();
                for (int l = 0; l < synNonDelBurdenCount3.length; l++) {
                    for (int m = 0; m < synNonDelBurdenCount3[l].length; m++) {
                        burdenCountList.add(synNonDelBurdenCount3[l][m]);
                    }
                }

                this.burdenCount[i][k]=burdenCountList.toArray();
            }
        }
    }

    //genotypeList 0: 0/0 1: 1/1 2: 0/1 0为ancestral, 1为derived, 2为杂合
    private int[] calculateNumNumDerivedNumHeter(TIntArrayList genotypeList){
        int[] genotypeCount= new int[3]; //num,numDerived, numHeter
        genotypeCount[0]=genotypeList.size();
        for (int i = 0; i < genotypeList.size(); i++) {
            if (genotypeList.get(i)==1){
                genotypeCount[1]++;
            }else if (genotypeList.get(i)==2){
                genotypeCount[2]++;
            }
        }
        return genotypeCount;
    }

    public String toString(){
        this.calculate();
        StringBuilder sb = new StringBuilder();
        // dim1: fdBins, dim2: p3, dim3: burden 9 column
        for (int i = 0; i < this.burdenCount.length; i++) {
            for (int k = 0; k < burdenCount[i].length; k++) {
                if ((!RefV1Utils.getChromosome(chrID,0).contains("D") && k==3)) continue;
                if ((RefV1Utils.getChromosome(chrID,0).contains("D") && k!=3)) continue;
                sb.append(taxon).append("\t");
                sb.append(chrID).append("\t");
                sb.append((i+1)/(double)10).append("\t");
                sb.append(IntrogressionRegion.P3.newInstanceFrom(k).name()).append("\t");
                sb.append(Joiner.on("\t").join(Ints.asList(burdenCount[i][k])));
                sb.append("\n");
            }
        }
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }
}
