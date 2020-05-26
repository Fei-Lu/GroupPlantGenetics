/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.wheat_epi;

import pgl.infra.dna.FastaByte;
import pgl.infra.pos.ChrPos;
import pgl.infra.range.Range;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import pgl.infra.utils.PStringUtils;



/**
 * Class holding protein-coding gene feature for fast access. Annotation information from GFF format.
 * @author Fei Lu
 */
public class WheatGeneFeature {
    Gene[] genes;
    //0 sort by position, 1 by sort by name
    int sortType = 0;
    public WheatGeneFeature () {}
    
    /**
     * Constructs a object from reading pgf (key gene feature) format 
     * @param infileS 
     */
     public WheatGeneFeature (String infileS , FastaByte outfileS1, String outfileS2) {
        this.readFromWheatNumGFF(infileS);
        //this.writeFile(outfileS);
        this.writeGeneSequence(outfileS1, outfileS2);
    }
    
     public WheatGeneFeature (String infileS ,String outfileS) {
        this.readFromWheatNumGFF(infileS);
        this.writeFile(outfileS);
    }
    
    
    /**
     * Write genomic sequence of gene, from start to the end, no flip if the gene is in the minus direction
     * @param genomef
     * @param outfileS 
     */
    public void writeGeneSequence (FastaByte genomef, String outfileS) {
        genomef.sortByName();
        try {
            BufferedWriter bw = XueboIOUtils.getTextWriter(outfileS);
            for (int i = 0; i < this.getGeneNumber(); i++) {
                String title = String.valueOf(this.getGeneChromosome(i))+"_"+String.valueOf(this.getGeneStart(i)+"_"+String.valueOf(this.getGeneEnd(i))+"_"+String.valueOf(this.getGeneName(i)));   
                int index = this.getGeneChromosome(i);
                int chrIndex = genomef.getIndexByName(String.valueOf(index));
                String chrseq = genomef.getSeq(chrIndex);
                String geneSeq = chrseq.substring(this.getGeneStart(i)-1, this.getGeneEnd(i)-1);
                String[] geneSeqs = PStringUtils.getMultilineString(60, geneSeq);
                bw.write(">"+title);
                bw.newLine();
                for (int j = 0; j < geneSeqs.length; j++) {
                    bw.write(geneSeqs[j]);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Write pgf file of gene annotation
     * @param outfileS 
     */
    public void writeFile (String outfileS) {
        try {
            BufferedWriter bw = XueboIOUtils.getTextWriter(outfileS);
            bw.write("GeneNumber\t"+String.valueOf(this.getGeneNumber()));
            bw.newLine();
            for (int i = 0; i < this.getGeneNumber(); i++) {
                StringBuilder sb = new StringBuilder();
                //sb.append("Gene\t").append(this.getGeneName(i)).append("\t").append(this.getGeneChromosome(i)).append("\t").append(this.getGeneStart(i)).append("\t").append(this.getGeneEnd(i)).append("\t").append(this.getGeneStrand(i));
                sb.append("Gene\t").append(this.getGeneName(i)).append("\t").append(this.getGeneChromosome(i)).append("\t").append(this.getGeneStart(i)).append("\t").append(this.getGeneEnd(i)).append("\t").append(this.getGeneStrand(i));
                
                //sb.append("\t").append(this.getGeneBiotype(i)).append("\t").append(this.getGeneDescription(i));
                System.out.println(this.getGeneName(i));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Return transcription start site (TSS) of a gene, null if 5'UTR does not exist.
     * @param geneIndex
     * @return 
     */
    public ChrPos getTSSOfGene (int geneIndex) {
        int index = this.getLongestTranscriptIndex(geneIndex);
        List<Range> utr5 = this.get5UTRList(geneIndex, index);
        if (utr5.isEmpty()) return null;
        if (this.getTranscriptStrand(geneIndex, index) == 1) {
            return new ChrPos((short)this.getGeneChromosome(geneIndex),utr5.get(0).getRangeStart());
        }
        else {
            return new ChrPos((short)this.getGeneChromosome(geneIndex),utr5.get(utr5.size()-1).getRangeEnd()-1);
        }
    }
    
    public int getCDSIndex (int geneIndex, int transcriptIndex, int chr, int pos) {
        return this.genes[geneIndex].ts.get(transcriptIndex).getCDSIndex(chr, pos);
    }
    
    /**
     * Return gene index from of a position, make sure the genes are sorted by position
     * @param chr
     * @param pos
     * @return negative value if the position is not in the range of any gene
     */
    public int getGeneIndex (int chr, int pos) {
        if (this.sortType != 0) {
            System.out.println("Genes needs to be sorted by position, program quits");
            System.exit(0);
        }
        Gene query = new Gene(chr, pos, pos+1);
        int hit = Arrays.binarySearch(this.genes, query);
        int index = hit;
        if (index < -1) {
            index = -index-2;
            if (this.isWithinThisGene(index, chr, pos)) return index;
        }
        return hit;
    }
    
    /**
     * Return if a position belong to a gene model
     * @param geneIndex
     * @param chr
     * @param pos
     * @return 
     */
    public boolean isWithinThisGene (int geneIndex, int chr, int pos) {
        if (chr != this.getGeneChromosome(geneIndex)) return false;
        if (pos < this.getGeneStart(geneIndex)) return false;
        if (pos >= this.getGeneEnd(geneIndex)) return false;
        return true;
    }
    
    /**
     * Return index of a gene, make sure the genes are sorted by name first
     * @param geneName
     * @return negative value if the gene name is not found
     */
    public int getGeneIndex (String geneName) {
        if (this.sortType != 1) {
            System.out.println("Genes needs to be sorted by name, program quits");
            System.exit(0);
        }
        return Arrays.binarySearch(genes, new Gene(geneName));
    }
    
    /**
     * Return the starting index of the gene on a chromosome, inclusive
     * @param chromosome
     * @return 
     */
    public int getStartIndexOfChromosome (int chromosome) {
        if (this.sortType != 0) {
            System.out.println("Genes needs to be sorted by position, program quits");
            System.exit(0);
        }
        Gene query  = new Gene ("Query", chromosome, Integer.MIN_VALUE, Integer.MIN_VALUE, Byte.MIN_VALUE);
        int hit  = Arrays.binarySearch(genes, query);
        int index = -hit-1;
        if (this.getGeneChromosome(index) == chromosome) return index;
        return hit;
    }
    
    /**
     * Return the ending index of the gene on a chromosome, exclusive
     * @param chromosome
     * @return 
     */
    public int getEndIndexOfChromosome (int chromosome) {
        if (this.sortType != 0) {
            System.out.println("Genes needs to be sorted by position, program quits");
            System.exit(0);
        }
        Gene query  = new Gene ("Query", chromosome+1, Integer.MIN_VALUE, Integer.MIN_VALUE, Byte.MIN_VALUE);
        int hit  = Arrays.binarySearch(genes, query);
        int index = -hit-1;
        if (this.getGeneChromosome(index-1) == chromosome) return index;
        return hit;
    }
    
    /**
     * Return the name of the gene
     * @param index
     * @return 
     */
    public String getGeneName(int index) {
        return genes[index].geneName;
    }
    
    /**
     * Return the chromosome number of the gene
     * @param index
     * @return 
     */
    public int getGeneChromosome (int index) {
        return genes[index].geneRange.chr;
    }
    
    /**
     * Return the starting position of the gene
     * @param index
     * @return 
     */
    public int getGeneStart (int index) {
        return genes[index].geneRange.start;
    }
    
    /**
     * Return the end position of the gene
     * @param index
     * @return 
     */
    public int getGeneEnd (int index) {
        return genes[index].geneRange.end;
    }
    
    /**
     * Return the strand of the gene, 1 is plus, 0 is minus
     * @param index
     * @return 
     */
    public byte getGeneStrand (int index) {
        return genes[index].strand;
    }
    
    /**
     * Return the biotype of the gene
     * @param index
     * @return 
     */
//    public String getGeneBiotype (int index) {
//        return genes[index].biotype;
//    }
    
    /**
     * Return the description of the gene
     * @param index
     * @return 
     */
//    public String getGeneDescription (int index) {
//        return genes[index].description;
//    }
    
    /**
     * Return the name of transcript
     * @param i
     * @param j
     * @return 
     */
    public String getTranscriptName (int i, int j) {
        return genes[i].ts.get(j).transcriptName;
    }
    
    /**
     * Return the name of transcript
     * @param i
     * @param j
     * @return 
     */
    public int getTranscriptChromosome (int i, int j) {
        return genes[i].ts.get(j).transcriptRange.chr;
    }
    
    /**
     * Return the start position of the transcript
     * @param i
     * @param j
     * @return 
     */
    public int getTranscriptStart (int i, int j) {
        return genes[i].ts.get(j).transcriptRange.start;
    }
    
    /**
     * Return the end position of transcript
     * @param i
     * @param j
     * @return 
     */
    public int getTranscriptEnd (int i, int j) {
        return genes[i].ts.get(j).transcriptRange.end;
    }
    
    /**
     * Return the strand of the transcript, 1 is plus, 0 is minus
     * @param i
     * @param j
     * @return 
     */
    public byte getTranscriptStrand (int i, int j) {
        return genes[i].ts.get(j).strand;
    }
    
    /**
     * Return a range list of CDS, return an empty list if there is no CDS
     * @param i
     * @param j
     * @return 
     */
    public List<Range> getCDSList (int i, int j) {
        return genes[i].ts.get(j).cdsList;
    }
    
    /**
     * Return a range list of 5UTR, return an empty list if there is no 5UTR
     * @param i
     * @param j
     * @return 
     */
    public List<Range> get5UTRList (int i, int j) {
        return genes[i].ts.get(j).utr5List;
    }
    
    private String getRangePositionString (List<Range> rList) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < rList.size(); i++) {
            sb.append(rList.get(i).getRangeStart()).append(":").append(rList.get(i).getRangeEnd()).append(";");
        }
        if (rList.size()!=0)sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }
    /**
     * Return a range list of 3UTR, return an empty list if there is no 3UTR
     * @param i
     * @param j
     * @return 
     */
    public List<Range> get3UTRList (int i, int j) {
        return genes[i].ts.get(j).utr3List;
    }
    
    /**
     * Return intron Range list, return an empty list if there is no intron
     * @param i
     * @param j
     * @return 
     */
    public List<Range> getIntronList (int i, int j) {
        return genes[i].ts.get(j).intronList;
    }
    
    /**
     * Return a string of CDS positions for output, return empty when there is no CDS
     * @param i
     * @param j
     * @return 
     */
    public String getCDSPositionString (int i, int j) {
        return getRangePositionString(this.getCDSList(i, j));
    }
    
    /**
     * Return a string of intron positions, return empty when there is no intron
     * @param i
     * @param j
     * @return 
     */
    public String getIntronPositionString (int i, int j) {
        return getRangePositionString(this.getIntronList(i, j));
    }
    
    /**
     * Return a string of 5UTR positions, return empty when there is no 5UTR
     * @param i
     * @param j
     * @return 
     */
    public String get5UTRPositionString (int i, int j) {
        return getRangePositionString(this.get5UTRList(i, j));
    }
    
    /**
     * Return a string of 3UTR positions, return empty when there is no 3UTR
     * @param i
     * @param j
     * @return 
     */
    public String get3UTRPositionString (int i, int j) {
        return getRangePositionString(this.get3UTRList(i, j));
    }
    
    /**
     * Return if there is 5UTR
     * @param i
     * @param j
     * @return 
     */
    public boolean isThere5UTR (int i, int j) {
        return genes[i].ts.get(j).utr5List.isEmpty() ? false:true;
    }
    
    /**
     * Return if there is 3UTR
     * @param i
     * @param j
     * @return 
     */
    public boolean isThere3UTR (int i, int j) {
        return genes[i].ts.get(j).utr3List.isEmpty() ? false:true;
    }
    
    /**
     * Read from Maize GFF file (ftp://ftp.ensemblgenomes.org/pub/plants/release-38/gff3/zea_mays)
     * @param infileS 
     */
    public void readFromWheatNumGFF (String infileS) {
        try {
            BufferedReader br;
            if (infileS.endsWith("gz")) br = XueboIOUtils.getTextGzipReader(infileS);
            else br = XueboIOUtils.getTextReader(infileS);
            //String temp = br.readLine();
            String temp = null;
            ArrayList<String> infoList = new ArrayList();
            ArrayList<String> geneList = new ArrayList();
            String[] tem = null;
            while ((temp = br.readLine()) != null) {
                //char s = temp.charAt(0);
                //if ((int)s < 48 || (int)s > 89) continue;
                List<String> tList= PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                //if (tem[2].startsWith("exon")) continue;
                //if (tem[2].startsWith("chromosome")) continue;
                //System.out.println(tem[2]);
                if (tem[2].startsWith("gene")) {
                    String[] te = tem[8].split(";");
                    geneList.add(te[0].split("=")[1]);
                }
                infoList.add(temp);
            }
            String[] geneNames = geneList.toArray(new String[geneList.size()]);
            Arrays.sort(geneNames);
            genes = new Gene[geneNames.length];
            String[] info = infoList.toArray(new String[infoList.size()]);
            for (int i = 0; i < info.length; i++) {
                tem = info[i].split("\t");
                if (tem[2].startsWith("gene")) {
                    String[] te = tem[8].split(";");
                    String query = te[0].split("=")[1];
                    int index = Arrays.binarySearch(geneNames, query);
                    genes[index] = new Gene (query, Integer.valueOf(tem[0]), Integer.valueOf(tem[3]), Integer.valueOf(tem[4])+1, (byte)(tem[6].equals("+")? 1:0));
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.sortGeneByStartPosition();
    }
    
    /**
     * Return number of genes
     * @return 
     */
    public int getGeneNumber () {
        return genes.length;
    }
    
    /**
     * Return transcript number of a gene
     * @param index
     * @return 
     */
    public int getTranscriptNumber (int index) {
        return genes[index].getTranscriptNumber();
    }
    
    /**
     * Return the index of the longest transcript of a gene
     * @param index
     * @return 
     */
    public int getLongestTranscriptIndex (int index) {
        return genes[index].longestTranscriptIndex;
    }
    
    /**
     * Return transcript number of all genes
     * @return 
     */
    public int getTotalTranscriptNumber () {
        int n = 0;
        for (int i = 0; i < this.getGeneNumber(); i++) n+=this.getTranscriptNumber(i);
        return n;
    }
    
    public void sortGeneByName () {
        this.sortType = 1;
        Arrays.sort(genes);
    }
    
    /**
     * Sort genes by their chromosome positions
     */
    public void sortGeneByStartPosition () {
        this.sortType = 0;
        Arrays.sort(genes);
    }
    
    class Gene implements Comparable<Gene> {
        String geneName = null;
        Range geneRange = null;
        byte strand = Byte.MIN_VALUE;
        //String biotype = null;
        //String description = null;
        ArrayList<Transcript> ts = new ArrayList();
        int longestTranscriptIndex = -1;
        
        public Gene (String geneName, int chr, int start, int end, byte strand) {
            this.geneName = geneName;
            this.strand = strand;
            geneRange = new Range(chr, start, end);
            //this.biotype = biotype;
            //this.description = discription;
        }
        
        public Gene (String geneName) {
            this.geneName = geneName;
        }
        
        public Gene (int chr, int start, int end) {
            geneRange = new Range(chr, start, end);
        }
        
        public void addTranscript (Transcript t) {
            ts.add(t);
        }
        
        public int getTranscriptNumber () {
            return ts.size();
        }
        
        public int getLongestTranscriptIndex () {
            return this.longestTranscriptIndex;
        }
        
        public int getTranscriptIndex (String transcriptName) {
            return Collections.binarySearch(ts, new Transcript(transcriptName));
        }
        
        public void calculateLongestTranscriptIndex () {
            int index = -1;
            int len = -1;
            for (int i = 0; i < this.getTranscriptNumber(); i++) {
                if (ts.get(i).transcriptRange.getRangeSize()>len) {
                    len = ts.get(i).transcriptRange.getRangeSize();
                    index = i;
                }
            }
            this.longestTranscriptIndex = index;
        }
        
        public void setLongestTranscriptIndex (int index) {
            this.longestTranscriptIndex = index;
        }
        
        public void sortTranscriptsByName () {
            Collections.sort(ts);
        }
        
        @Override
        public int compareTo(Gene t) {
            if (sortType == 0) {
                return geneRange.compareTo(t.geneRange);
            }
            else if (sortType == 1) {
                return this.geneName.compareTo(t.geneName);
            }
            return 0;
        }
    }
    
    class Transcript implements Comparable<Transcript> {
        String transcriptName = null;
        Range transcriptRange = null;
        byte strand = Byte.MIN_VALUE;
        List<Range> cdsList = new ArrayList();
        List<Range> intronList = new ArrayList();
        List<Range> utr5List = new ArrayList();
        List<Range> utr3List = new ArrayList();
        
        public Transcript (String transcriptName) {
            this.transcriptName = transcriptName;
        }
        
        public Transcript (String transcriptName, int chr, int start, int end, byte strand) {
            this.transcriptName = transcriptName;
            this.strand = strand;
            transcriptRange = new Range(chr, start, end);
        }
        
        public void addCDS (int chr, int start, int end) {
            cdsList.add(new Range(chr, start, end));
        }
        
        public void add5UTR (int chr, int start, int end) {
            utr5List.add(new Range(chr, start, end));
        }
        
        public void add3UTR (int chr, int start, int end) {
            utr3List.add(new Range(chr, start, end));
        }
        
        public int getChromosome () {
            return this.transcriptRange.chr;
        }
        
        public void sortCDSByPosition () {
            Collections.sort(cdsList);
        }
        
        public void sort5UTRByPosition () {
            if (utr5List.isEmpty()) return;
            Collections.sort(utr5List);
        }
        
        public void sort3UTRByPosition () {
            if (utr3List.isEmpty()) return;
            Collections.sort(utr5List);
        }
        
        public int getCDSIndex (int chr, int pos) {
            Range query = new Range(chr, pos, pos+1);
            int hit = Collections.binarySearch(cdsList, query);
            int index = hit;
            if (index < -1) {
                index = -index-2;
                if (this.isWithinThisCDS(index, chr, pos)) return index;
            }
            return hit;
        }
        
        public boolean isWithinThisCDS (int cdsIndex, int chr, int pos) {
            if (cdsList.get(cdsIndex).chr != chr) return false;
            if (pos < cdsList.get(cdsIndex).start) return false;
            if (pos >= cdsList.get(cdsIndex).end) return false;
            return true;
        }
        
        public void calculateIntron () {
            if (cdsList.size() < 2) return;
            for (int i = 0; i < cdsList.size()-1; i++) {
                Range r = cdsList.get(i);
                Range nr = cdsList.get(i+1);
                intronList.add(new Range(r.chr, r.end, nr.start ));
            }
        }

        @Override
        public int compareTo(Transcript t) {
            return transcriptName.compareTo(t.transcriptName);
        }
    }
}

