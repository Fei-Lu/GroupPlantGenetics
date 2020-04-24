package daxing.load.neutralSiteLoad;

import daxing.common.PGF;
import pgl.infra.range.Range;
import pgl.infra.utils.wheat.RefV1Utils;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GenesDB {

    GeneRange[][] geneRanges;

    public GenesDB(List<String> geneNames, String pgfFile){
        PGF pgf=new PGF(pgfFile);
        pgf.sortGeneByName();
        int chr, geneIndex, geneStart, geneEnd, geneMiddlePos, start, end;
        String geneName;
        GeneRange geneRange;
        List<GeneRange>[] geneRangeListArray=new List[42];
        for (int i = 0; i < geneRangeListArray.length; i++) {
            geneRangeListArray[i]=new ArrayList<>();
        }
        for (int i = 0; i < geneNames.size(); i++) {
            geneName=geneNames.get(i);
            geneIndex=pgf.getGeneIndex(geneName);
            chr=pgf.getGeneChromosome(geneIndex);
            geneStart=pgf.getGeneStart(geneIndex);
            geneEnd=pgf.getGeneEnd(geneIndex);
            geneMiddlePos=(geneStart+geneEnd)/2;
            start=geneMiddlePos-50000 < 0 ? 0 : geneMiddlePos-50000;
            end=geneMiddlePos+50000 > RefV1Utils.getChrIDLength(chr) ? RefV1Utils.getChrIDLength(chr) : geneMiddlePos+50000;
            geneRange=new GeneRange(geneName, chr, start, end);
            geneRangeListArray[chr-1].add(geneRange);
        }
        GeneRange[][] geneRanges=new GeneRange[42][];
        for (int i = 0; i < geneRanges.length; i++) {
            geneRanges[i]=new GeneRange[geneRangeListArray[i].size()];
            for (int j = 0; j < geneRanges[i].length; j++) {
                geneRanges[i][j]=geneRangeListArray[i].get(j);
            }
        }
        this.geneRanges=geneRanges;
        this.sortByGeneRange();
    }

    private void sortByGeneRange(){
        for (int i = 0; i < this.geneRanges.length; i++) {
            Arrays.sort(geneRanges[i]);
        }
    }

    public GeneRange[] getChrGeneRange(int chr){
        return this.geneRanges[chr-1];
    }

    public class GeneRange extends Range{

        String geneName;

        GeneRange(String geneName, int chr, int start, int end){
            super(chr, start, end);
            this.geneName=geneName;
        }

        public boolean contain(int chr, int pos){
            return isContain(chr, pos);
        }

    }
}
