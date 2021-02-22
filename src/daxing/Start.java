package daxing;

import daxing.applets.Ne;
import daxing.common.Group;
import daxing.common.SubgenomeCombination;

public class Start {

    public static void main(String[] args) {
        String vcfDir=args[0];
        String taxaInfo_depth=args[1];
        String pseudoDiploidDir=args[2];
        Group group=Group.valueOf(args[3]);
        int pseudoDiploidNum=Integer.parseInt(args[4]);
        String ancestralDir=args[5];
        String ancestralVCFDir=args[6];

        Ne.randomSyntheticPseudoDiploid(vcfDir, taxaInfo_depth, pseudoDiploidDir, group,
                SubgenomeCombination.D, pseudoDiploidNum);
        Ne.randomSyntheticPseudoDiploid(vcfDir, taxaInfo_depth, pseudoDiploidDir, group,
                SubgenomeCombination.AB, pseudoDiploidNum);
        Ne.transformRefToAncestral(pseudoDiploidDir, ancestralDir, ancestralVCFDir);
    }
}