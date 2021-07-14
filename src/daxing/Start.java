package daxing;

import daxing.common.MD5;
import daxing.common.Matrix;

public class Start {

    public static void main(String[] args) {

//        SNPEff.extractEffectAndImpact(args[0],args[1]);
//        VEP.extractEffectAndImpact(args[0],args[1]);
//        ComplementaryGo.start();
//        MD5.checkMD5("/Volumes/LuLab3T-41/Wheat1700/wheat1700_remain.md5.txt");
//        String inputFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/005_vmap2_1000/001_GermplasmDetermination/d_ibs.txt";
//        String outFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/005_vmap2_1000" +
//                "/001_GermplasmDetermination/d_ibs.pairwise.txt";
//        IBS(inputFile, outFile);
        MD5.getMD5FromDir(args[0]);
    }

    public static void IBS(String ibsFile, String outFile){
        Matrix m = new Matrix(ibsFile, false);
        m.writePairWise(outFile);

    }


}