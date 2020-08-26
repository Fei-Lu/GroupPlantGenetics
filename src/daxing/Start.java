package daxing;

import daxing.load.ancestralSite.LoadGO;

public class Start {

    public static void main(String[] args) {
//        String inputDir="/Users/xudaxing/Documents/deleteriousMutation/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/003_retainTriad/009_retainTriadHexaploid_expected_pos_slidingWindow/002_slidWindow";
//        String outDir="/Users/xudaxing/Documents/deleteriousMutation/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/003_retainTriad/009_retainTriadHexaploid_expected_pos_slidingWindow/004_mergeSynNonDel";
//        WindowLoadABD.mergeIndividual(inputDir, outDir);
//        IOTool.viewHeader("/Users/xudaxing/Documents/deleteriousMutation/003_vmap2.1_20200628/002_deleteriousPosLib/002_exonAnnotationByDerivedSift/001_exonSNPAnnotationByChrID/chr001_SNP_anno.txt.gz");
        LoadGO.start();
    }
}