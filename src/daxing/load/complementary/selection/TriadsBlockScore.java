package daxing.load.complementary.selection;

import daxing.common.ChrRange;
import daxing.load.complementary.TriadsBlock;
import gnu.trove.list.array.TDoubleArrayList;
import java.util.List;

/**
 * selection of triads A B D using TriadsBlock class
 */
public class TriadsBlockScore extends TriadsBlock {

    /**
     * subgenome xpclr score list;
     */
    TDoubleArrayList[] xpclrScoreListArray;

    public TriadsBlockScore(String triadsID, String[] geneName, int[] triadsCDSLen, List<String>[] blockGeneName,
                            ChrRange[] chrRange) {
        super(triadsID, geneName, triadsCDSLen, blockGeneName, chrRange);
        this.xpclrScoreListArray=new TDoubleArrayList[3];
        for (int i = 0; i < this.xpclrScoreListArray.length; i++) {
            this.xpclrScoreListArray[i]=new TDoubleArrayList();
        }
    }

    public TriadsBlockScore(ChrRange chrRange){
        super(chrRange);
    }

    public void addXPCLRScore(DominanceAdditiveSelection.Sub sub, double score){
        this.xpclrScoreListArray[sub.getIndex()].add(score);
    }

    public double[] calculateXPCLRScore(){
        double[] abdScore=new double[DominanceAdditiveSelection.Sub.values().length];
        for (int i = 0; i < this.xpclrScoreListArray.length; i++) {
            abdScore[i]=(this.xpclrScoreListArray[i].sum())/(this.xpclrScoreListArray[i].size());
        }
        return abdScore;
    }

}
