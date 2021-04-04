package daxing;

import daxing.load.complementary.TriadsBlockUtils;

public class Start {

    public static void main(String[] args) {
//        String inputDir=args[0];
//        String outFile=args[1];
//        String file="/Users/xudaxing/Desktop/Fig5_som/triadsLoadHexaploidPseudo.txt.gz";
//        String outFile="/Users/xudaxing/Desktop/res.txt";
//        WilcoxonRankUtils wilcoxonRankUtils=new WilcoxonRankUtils(file);
//        wilcoxonRankUtils.writeWilcoxonSignedRank(outFile);
//        System.out.println();
        String triadsBlockInputFile="/Users/xudaxing/Desktop/triadsBlock.txt.gz";
        String pgfFile="/Users/xudaxing/Data/wheatReference/v1.0/annotation/wheat_v1.1_Lulab.pgf";
        String triadsBlockOutFile="/Users/xudaxing/Desktop/triadsBlockChrRange.txt.gz";
        TriadsBlockUtils.writeTriadsBlock(triadsBlockInputFile, pgfFile, triadsBlockOutFile);
    }
}