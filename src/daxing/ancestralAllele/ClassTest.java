package daxing.ancestralAllele;

import daxing.common.ChrConvertionRule;
import format.position.ChrPos;

import java.nio.file.Paths;
import java.util.Map;

public class ClassTest {

    public static void main(String[] args) {
        int b=2;
        ChrConvertionRule chrConvertionRule=new ChrConvertionRule(Paths.get("/Users/xudaxing/Desktop/chrConvertionRule.txt"));
        MAF maf=new MAF(1, chrConvertionRule, Paths.get("/Users/xudaxing/Desktop/test"));
        AllelesInfor allelesInfor=new AllelesInfor(Paths.get("/Users/xudaxing/Desktop/abd/chr001.ABDgenome.filterMiss_subset.vcf"));
        Map<ChrPos, String[]> map=maf.getAllele(allelesInfor);
        b=3;
    }
}
