package daxing.individualIntrogression.individual;

import com.google.common.base.Joiner;
import com.google.common.primitives.Ints;
import daxing.common.factors.LoadType;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.wheat.RefV1Utils;
import java.util.Arrays;

public class IndividualBurden_individualFd {

    //0: 0/0 1: 1/1 2: 0/1 0为ancestral, 1为derived, 2为杂合
    String taxon;
    int chrID;
    TIntArrayList[][] genotypeCount; // dim1: Donor, dim2: loadGroup
    int[][] burdenCount;         // dim1: Donor, dim2: burden 9 column

    public IndividualBurden_individualFd(String taxon, int chrID){
        TIntArrayList[][] genotypeCount=new TIntArrayList[Donor.values().length][];
        int[][] burdenCount= new int[Donor.values().length][];
        for (int j = 0; j < genotypeCount.length; j++) {
            genotypeCount[j] = new TIntArrayList[LoadType.values().length];
            burdenCount[j] = new int[9];
            for (int i = 0; i < genotypeCount[j].length; i++) {
                genotypeCount[j][i]=new TIntArrayList();
            }
            Arrays.fill(burdenCount[j], -1);
        }
        this.chrID= chrID;
        this.genotypeCount=genotypeCount;
        this.taxon=taxon;
        this.burdenCount=burdenCount;
    }

    public void addLoadCount(Donor donor, byte[] indexGenotype){
        this.genotypeCount[donor.getIndex()][indexGenotype[0]].add(indexGenotype[1]);
    }

    private void calculate(){
        // genotypeCount: dim1: Donor, dim2: loadGroup
        int[][] synNonDelBurdenCount3;
        TIntArrayList burdenCountList; //        int numSyn, numDerivedInSyn, numHeterInSyn, numNonsyn, numDerivedInNonsyn, numHeterInNonsyn, numHGDeleterious
//                , numDerivedInHGDeleterious, numHeterInHGDeleterious;
        for (int i = 0; i < genotypeCount.length; i++) {
            synNonDelBurdenCount3 = new int[LoadType.values().length][];
            for (int l = 0; l < synNonDelBurdenCount3.length; l++) {
                synNonDelBurdenCount3[l]=new int[3];
            }
            for (int l = 0; l < genotypeCount[i].length; l++) {
                synNonDelBurdenCount3[l]=calculateNumNumDerivedNumHeter(genotypeCount[i][l]);
            }
            burdenCountList = new TIntArrayList();
            for (int l = 0; l < synNonDelBurdenCount3.length; l++) {
                for (int m = 0; m < synNonDelBurdenCount3[l].length; m++) {
                    burdenCountList.add(synNonDelBurdenCount3[l][m]);
                }
            }

            this.burdenCount[i]=burdenCountList.toArray();
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
        // dim1: Donor, dim2: burden 9 column
        for (int i = 0; i < this.burdenCount.length; i++) {
            if (((!RefV1Utils.getChromosome(chrID,0).contains("D")) && i==3)) continue;
            if ((RefV1Utils.getChromosome(chrID,0).contains("D") && i==0)) continue;
            if ((RefV1Utils.getChromosome(chrID,0).contains("D") && i==1)) continue;
            if ((RefV1Utils.getChromosome(chrID,0).contains("D") && i==2)) continue;
            sb.append(taxon).append("\t");
            sb.append(chrID).append("\t");
            sb.append(Donor.values()[i]).append("\t");
            sb.append(Joiner.on("\t").join(Ints.asList(burdenCount[i])));
            sb.append("\n");
        }
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }
}
