package daxing.load;


import daxing.common.IOTool;
import daxing.common.PGF;
import daxing.common.RowTableTool;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class TranscriptDB {

    private String[][] transcriptNames;
    private int[][] cdsLen;
    private List<SNPAnnotation>[][] transcriptSNPAnno;

    public TranscriptDB(String pgfFile, String exonSNPAnnoDir){
        List<File> exonAnnoFiles=IOUtils.getVisibleFileListInDir(exonSNPAnnoDir);
        String[][] transcriptNames=this.getTranscriptNames(exonAnnoFiles);
        PGF pgf=new PGF(pgfFile);
        pgf.sortGeneByName();
        int[][] cdsLen=this.getCdsLen(transcriptNames, pgf);
        List<SNPAnnotation>[][] transcriptSNPAnno=new List[cdsLen.length][];
        for (int i = 0; i < exonAnnoFiles.size(); i++) {
            transcriptSNPAnno[i]=this.getTranscriptSNPAnno(exonAnnoFiles.get(i), transcriptNames[i]);
        }
        this.transcriptNames=transcriptNames;
        this.cdsLen=cdsLen;
        this.transcriptSNPAnno=transcriptSNPAnno;
    }

    private String[][] getTranscriptNames(List<File> exonAnnoFiles){
        String[][] transcriptNames=new String[42][];
        Set<String> geneNameSet;
        List<String> geneNamesList;
        for (int i = 0; i < exonAnnoFiles.size(); i++) {
            geneNameSet=RowTableTool.getColumnSet(exonAnnoFiles.get(i).getAbsolutePath(), 10);
            geneNamesList=new ArrayList<>(geneNameSet);
            Collections.sort(geneNamesList);
            transcriptNames[i]=geneNamesList.stream().toArray(String[]::new);
        }
        return transcriptNames;
    }

    private int[][] getCdsLen(String[][] transcriptNames, PGF pgf){
        int[][] cdsLen=new int[transcriptNames.length][];
        for (int i = 0; i < cdsLen.length; i++) {
            cdsLen[i]=new int[transcriptNames[i].length];
        }
        int geneIndex, longestTranscriptIndex;
        for (int i = 0; i < transcriptNames.length; i++) {
            for (int j = 0; j < transcriptNames[i].length; j++) {
                geneIndex=pgf.getGeneIndex(transcriptNames[i][j].substring(0,18));
                longestTranscriptIndex=pgf.getLongestTranscriptIndex(geneIndex);
                cdsLen[i][j]=pgf.getCDSLen(geneIndex, longestTranscriptIndex);
            }
        }
        return cdsLen;
    }

    private List<SNPAnnotation>[] getTranscriptSNPAnno(File exonSNPAnnoFile, String[] transcriptNames){
        List<SNPAnnotation>[] geneAnno=new List[transcriptNames.length];
        for (int i = 0; i < geneAnno.length; i++) {
            geneAnno[i]=new ArrayList<>();
        }
        try (BufferedReader br = IOTool.getReader(exonSNPAnnoFile)) {
            br.readLine();
            String line;
            List<String> temp;
            int index;
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                index=Arrays.binarySearch(transcriptNames, temp.get(10));
                geneAnno[index].add(getSNPAnnotation(line));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return geneAnno;
    }

    private SNPAnnotation getSNPAnnotation(String line){
        List<String> temp= PStringUtils.fastSplit(line);
        short chr=Short.parseShort(temp.get(1));
        int pos=Integer.parseInt(temp.get(2));
        char refBase=temp.get(3).charAt(0);
        char altBase=temp.get(4).charAt(0);
        char majorBase=temp.get(5).charAt(0);
        double[] aaf=new double[2];
        aaf[0]=Double.parseDouble(temp.get(8));
        aaf[1]=Double.parseDouble(temp.get(9));
        String transcriptName=temp.get(10);
        SNPAnnotation.Region region=SNPAnnotation.Region.valueOf(temp.get(11));
        SNPAnnotation.Variant_type variant_type=SNPAnnotation.Variant_type.valueOf(temp.get(12));
        String alt_SIFT=temp.get(13);
        String gerp=temp.get(20);
        String daf=temp.get(32);
        double maf=Double.parseDouble(temp.get(7));
        String[] dafs=new String[2];
        dafs[0]=temp.get(33);
        dafs[1]=temp.get(34);
        String ancestral=temp.get(31);
        float recombinationRate=Float.parseFloat(temp.get(21));
        return new SNPAnnotation(chr, pos, refBase, altBase, transcriptName, majorBase, ancestral, maf
                , aaf, daf, dafs, region, variant_type, alt_SIFT, gerp, recombinationRate);
    }

    public List<SNPAnnotation>[][] getTranscriptSNPAnno() {
        return transcriptSNPAnno;
    }

    public String[] getTranscriptNames(int chr){
        return this.transcriptNames[chr-1];
    }

    public String[] getGeneName(int chr){
        String[] transcriptNames=this.transcriptNames[chr-1];
        return Arrays.stream(transcriptNames).map(s->s.substring(0,18)).toArray(String[]::new);
    }

    public String getGeneName(int chr, int pos){
        String[] geneNames=this.getGeneName(chr);
        List<SNPAnnotation>[] genes=this.transcriptSNPAnno[chr-1];
        int index=-1;
        for (int i = 0; i < genes.length; i++) {
            for (int j = 0; j < genes[i].size(); j++) {
                if (genes[i].get(j).check((short) chr, pos)){
                    index=i;
                }
            }
        }
        return geneNames[index];
    }

    /**
     *
     * @param chr
     * @param pos
     * @return -1  表示不包含
     */
    public int[] getGeneIndexSNPIndex(int chr, int pos){
        int[] index=new int[2];
        for (int i = 0; i < index.length; i++) {
            index[i]=-1;
        }
        List<SNPAnnotation>[] genes=this.transcriptSNPAnno[chr-1];
        for (int i = 0; i < genes.length; i++) {
            for (int j = 0; j < genes[i].size(); j++) {
                if (genes[i].get(j).check((short) chr, pos)){
                    index[0]=i;
                    index[1]=j;
                }
            }
        }
        return index;
    }

    public SNPAnnotation getSNP(int chr, int pos){
        List<SNPAnnotation>[] snpAnnoList=this.getTranscriptSNPAnno()[chr-1];
        int[] index=getGeneIndexSNPIndex(chr, pos);
        return snpAnnoList[index[0]].get(index[1]);
    }

    public int getGeneNum(int chr){
        return this.getGeneName(chr).length;
    }

    public boolean contain(int chr, int pos){
        List<SNPAnnotation>[] genes=this.transcriptSNPAnno[chr-1];
        for (int i = 0; i < genes.length; i++) {
            for (int j = 0; j < genes[i].size(); j++) {
                if (genes[i].get(j).check((short) chr, pos)) return true;
            }
        }
        return false;
    }

    public boolean hasAncestral(int chr, int pos){
        return this.getSNP(chr, pos).hasAncestral();
    }

    public boolean isRefAlleleAncestral(int chr, int pos){
        return this.getSNP(chr, pos).isRefAlleleAncestral();
    }

    public boolean isSyn(int chr, int pos){
        return this.getSNP(chr, pos).isSyn();
    }

    public boolean isNonsyn(int chr, int pos){
        return this.getSNP(chr, pos).isNonSyn();
    }

    public boolean isDeleterious(int chr, int pos){
        return this.getSNP(chr, pos).isDeleterious();
    }

}
