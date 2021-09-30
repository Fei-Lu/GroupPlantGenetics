package daxing.temp;

import com.google.common.base.Joiner;
import daxing.common.LoadType;
import daxing.common.WheatLineage;
import gnu.trove.list.array.TIntArrayList;
import java.util.Arrays;

public class IndividualBurden {

    //0: 0/0 1: 1/1 2: 0/1 0为ancestral, 1为derived, 2为杂合
    String taxon;
    TIntArrayList[][][][] genotypeCount; // dim1: Sub, dim2: fdBins, dim3: p3, dim4: loadGroup
    int[][][][] burdenCount; // dim1: Sub, dim2: fdBins, dim3: p3, dim4: burden 9 column

    public IndividualBurden(String taxon, int fdBinsNum){
        TIntArrayList[][][][] genotypeCount=new TIntArrayList[WheatLineage.values().length][][][];
        for (int i = 0; i < genotypeCount.length; i++) {
            genotypeCount[i]=new TIntArrayList[fdBinsNum][][];
            for (int j = 0; j < genotypeCount[i].length; j++) {
                genotypeCount[i][j] = new TIntArrayList[IntrogressionRegion.P3.values().length][];
                for (int k = 0; k < genotypeCount[i][j].length; k++) {
                    genotypeCount[i][j][k]=new TIntArrayList[LoadType.values().length];
                    for (int l = 0; l < genotypeCount[i][j][k].length; l++) {
                        genotypeCount[i][j][k][l]=new TIntArrayList();
                    }
                }
            }
        }
        this.genotypeCount=genotypeCount;
        this.taxon=taxon;
    }

    public void addLoadCount(WheatLineage wheatLineage, double fdBin, IntrogressionRegion.P3 p3, byte[] indexGenotype){
        int fdBinIndex=((int)fdBin)*10-1;
        this.genotypeCount[wheatLineage.getIndex()][fdBinIndex][p3.index][indexGenotype[0]].add(indexGenotype[1]);
    }

    private void calculate(){
        int[][] synNonDelBurdenCount3; // dim1: loadGroup, dim2: num,numDerived, numHeter
        TIntArrayList burdenCountList; //        int numSyn, numDerivedInSyn, numHeterInSyn, numNonsyn, numDerivedInNonsyn, numHeterInNonsyn, numHGDeleterious
//                , numDerivedInHGDeleterious, numHeterInHGDeleterious;
        for (int i = 0; i < genotypeCount.length; i++) {
            for (int j = 0; j < genotypeCount[i].length; j++) {
                for (int k = 0; k < genotypeCount[i][j].length; k++) {
                    synNonDelBurdenCount3 = new int[LoadType.values().length][];
                    for (int l = 0; l < synNonDelBurdenCount3.length; l++) {
                        synNonDelBurdenCount3[l]=new int[3];
                    }
                    for (int l = 0; l < genotypeCount[i][j][k].length; l++) {
                        synNonDelBurdenCount3[l]=calculateNumNumDerivedNumHeter(genotypeCount[i][j][k][l]);
                    }
                    burdenCountList = new TIntArrayList();
                    for (int l = 0; l < synNonDelBurdenCount3.length; l++) {
                        for (int m = 0; m < synNonDelBurdenCount3[l].length; m++) {
                            burdenCountList.add(synNonDelBurdenCount3[l][m]);
                        }
                    }
                    this.burdenCount[i][j][k]=burdenCountList.toArray(new int[0]);
                }
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
        StringBuilder sb = new StringBuilder();
        // dim1: Sub, dim2: fdBins, dim3: p3, dim4: burden 9 column
        for (int i = 0; i < this.burdenCount.length; i++) {
            for (int j = 0; j < burdenCount[i].length; j++) {
                for (int k = 0; k < burdenCount[i][j].length; k++) {
                    sb.append(WheatLineage.newInstanceFrom(i).name()).append("\t");
                    sb.append((j+1)/10).append("\t");
                    sb.append(IntrogressionRegion.P3.newInstanceFrom(k).name()).append("\t");
                    sb.append(Joiner.on("\t").join(Arrays.asList(burdenCount[i][j][k])));
                    sb.append("\n");
                }
            }
        }
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }
}
