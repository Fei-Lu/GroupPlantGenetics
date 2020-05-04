package daxing.load.slidingWindow;

import daxing.common.PGF;
import pgl.infra.range.Range;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class GeneWindowDB {

    List<GeneWindow>[] geneWindows;

    public GeneWindowDB(List<String> geneNameList, String pgfFile){
        PGF pgfSortByGeneName=new PGF(pgfFile);
        PGF pgfSortByGeneRange=new PGF(pgfFile);
        pgfSortByGeneName.sortGeneByName();
        pgfSortByGeneRange.sortGeneByGeneRange();
        List<GeneWindow>[] geneWindows=new List[42];
        for (int i = 0; i < geneWindows.length; i++) {
            geneWindows[i]=new ArrayList<>();
        }
        int chr, pos,geneIndexByName,transcriptIndexByName,geneIndexByRange;
        int chrStartIndex, chrEndIndex, startIndex, endIndex;
        int geneIndex=0, count=0;
        String geneName;
        List<PGF.Gene> genes;
        GeneWindow geneWindow;
        List<Range> ranges;
        for (int i = 0; i < geneNameList.size(); i++) {
            geneName=geneNameList.get(i);
            geneIndexByName=pgfSortByGeneName.getGeneIndex(geneName);
            transcriptIndexByName=pgfSortByGeneName.getLongestTranscriptIndex(geneIndexByName);
            chr=pgfSortByGeneName.getGeneChromosome(geneIndexByName);
            ranges=pgfSortByGeneName.getCDSList(geneIndexByName, transcriptIndexByName);
            pos=ranges.get(0).start;
            geneIndexByRange=pgfSortByGeneRange.getGeneIndex(chr, pos);
            chrStartIndex=pgfSortByGeneRange.getStartIndexOfChromosome(chr);
            chrEndIndex=pgfSortByGeneRange.getEndIndexOfChromosome(chr);
            startIndex= geneIndexByRange-5 < chrStartIndex ? chrStartIndex : geneIndexByRange-5;
            endIndex= geneIndexByRange+6 > chrEndIndex ? chrEndIndex : geneIndexByRange+6;
            genes=new ArrayList<>();
            for (int j = startIndex; j < endIndex; j++) {
                genes.add(pgfSortByGeneRange.getGene(j));
                if (j==geneIndexByRange){
                    geneIndex=count;
                }else {
                    count++;
                }
            }
            count=0;
            geneWindow=new GeneWindow(genes, geneIndex);
            geneWindows[chr-1].add(geneWindow);
        }
        this.geneWindows=geneWindows;
        this.sortByGeneRange();
    }

    public List<GeneWindow> getGeneWindow(int chr){
        return this.geneWindows[chr-1];
    }

    private void sortByGeneRange(){
        for (int i = 0; i < this.geneWindows.length; i++) {
            Collections.sort(this.geneWindows[i]);
        }
    }

    public class GeneWindow implements Comparable<GeneWindow> {

        private List<PGF.Gene> genes;
        private int geneIndex;

        public GeneWindow(List<PGF.Gene> genes, int geneIndex){
            this.genes=genes;
            this.geneIndex=geneIndex;
        }

        public List<Range> getWindowCDSRanges(){
            List<Range> cdsRangeListAllGenes=new ArrayList<>();
            List<Range> cdsRangeList;
            int longestTranscriptIndex;
            for (int i = 0; i < this.genes.size(); i++) {
                longestTranscriptIndex=this.genes.get(i).getLongestTranscriptIndex();
                cdsRangeList=this.genes.get(i).getTs().get(longestTranscriptIndex).getCdsList();
                cdsRangeListAllGenes.addAll(cdsRangeList);
            }
            Collections.sort(cdsRangeListAllGenes);
            return cdsRangeListAllGenes;
        }

        public int getGeneIndex() {
            return geneIndex;
        }

        public String getGeneName(){
            return this.genes.get(geneIndex).getGeneName();
        }

        public int getStartOfWindowCDSRange(){
            return this.getWindowCDSRanges().get(0).start;
        }

        public int getEndOfWindowCDSRange(){
            int size=this.genes.size();
            return this.getWindowCDSRanges().get(size-1).end;
        }

        @Override
        public int compareTo(GeneWindow o) {
            int thisGeneIndex=this.getGeneIndex();
            int oGeneIndex=o.getGeneIndex();
            return this.genes.get(thisGeneIndex).getGeneRange().compareTo(o.genes.get(oGeneIndex).getGeneRange());
        }

        public boolean containCDSPos(int chr, int pos){
            List<Range> genes=this.getWindowCDSRanges();
            for (int i = 0; i < genes.size(); i++) {
                if (genes.get(i).isContain(chr, pos)) return true;
            }
            return false;
        }
    }

}
