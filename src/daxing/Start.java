package daxing;

import daxing.common.PGF;
import pgl.infra.range.Range;
import pgl.infra.utils.wheat.RefV1Utils;

public class Start {

    public static void main(String[] args) {
//        int chrSize= RefV1Utils.getChromosomeLength("1A");
//        String inputFile="/Users/xudaxing/Documents/deleteriousMutation/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/002_derivedSift/005_test/003_hexaploid_test_pos_forSlidingWindow";
//        String outFile="/Users/xudaxing/Documents/deleteriousMutation/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/002_derivedSift/005_test/004_hexaploid_test_pos_slidingWindow";
//        int windowSize=10_000_000;
//        int stepSize=1_000_000;
//        Go.slidingWindow(inputFile, "A", windowSize, stepSize, outFile);
//        String inputDir="/Users/xudaxing/Documents/deleteriousMutation/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/002_derivedSift/005_test/004_hexaploid_test_pos_slidingWindow";
//        String outFile="/Users/xudaxing/Documents/deleteriousMutation/003_vmap2" +
//                ".1_20200628/004_deleterious/001_triadsSelection/002_derivedSift/005_test" +
//                "/005_hexaploid_test_pos_slidingWindow_mergeTaxon/triadA.10MbWindow.1MbStep.txt.gz";
//        Go.mergeTaxonDel(inputDir,outFile);
//        SlidingWindowForLoadComplement window=new SlidingWindowForLoadComplement(23, 5, 2);
//        System.out.println();
//        Go.start();
        String pgfFile="/Users/xudaxing/Data/wheatReference/v1.0/annotation/wheat_v1.1_Lulab.pgf";
//        String geneName="TraesCS5A02G391700";
        PGF pgf=new PGF(pgfFile);
        pgf.sortGeneByGeneRange();

//        String refChr="2D";
//        int refPos=33953530;
//        int chrID= RefV1Utils.getChrID(refChr, refPos);
//        int posOnChrID=RefV1Utils.getPosOnChrID(refChr, refPos);
//        int geneIndex=pgf.getGeneIndex(chrID, posOnChrID);
//        System.out.println(geneIndex);
//        Range geneRange=pgf.getGene(geneIndex).getGeneRange();
//        System.out.println(geneIndex);
//        System.out.println(geneRange.chr+" "+geneRange.getRangeStart()+" "+geneRange.end);

        int geneIndex=28919;
        Range geneRange=pgf.getGene(geneIndex).getGeneRange();
        String refChr= RefV1Utils.getChromosome(geneRange.chr, geneRange.start);
        int startPosRef=RefV1Utils.getPosOnChromosome(geneRange.chr, geneRange.start);
        int endPosRef=RefV1Utils.getPosOnChromosome(geneRange.chr, geneRange.end);
        System.out.println(refChr+" "+startPosRef+" "+endPosRef);

//        pgf.sortGeneByName();
//        int geneIndex=pgf.getGeneIndex(geneName);
//        int start=pgf.getGene(geneIndex).getGeneRange().start;
//        int end=pgf.getGene(geneIndex).getGeneRange().end;
//        int chrID=pgf.getGene(geneIndex).getGeneRange().getRangeChromosome();
//        String chr= RefV1Utils.getChromosome(chrID, start);
//        int startRef=RefV1Utils.getPosOnChromosome(chrID, start);
//        int endRef=RefV1Utils.getPosOnChromosome(chrID, end);
//        System.out.println(chr+"\t"+startRef+" "+endRef);
    }
}