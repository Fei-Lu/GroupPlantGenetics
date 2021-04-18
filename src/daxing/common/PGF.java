/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package daxing.common;

import edu.umd.cs.findbugs.annotations.SuppressFBWarnings;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import org.apache.commons.lang3.ArrayUtils;
import pgl.infra.dna.FastaByte;
import pgl.infra.dna.SequenceByte;
import pgl.infra.pos.ChrPos;
import pgl.infra.range.Range;
import pgl.infra.range.RangeInterface;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.IntStream;

/**
 *  modified from GeneFuture class.
 * @author Daxing Xu
 */
public class PGF {
    Gene[] genes;
    //0 sort by position, 1 by sort by name
    int sortType = 0;

    /**
     * Constructs a object from reading pgf (key gene feature) format
     * @param infileS
     */
    public PGF (String infileS) {
        this.readFile(infileS);
    }

    public PGF(Gene[] genes){
        this.genes=genes;
    }

    public Gene[] getGenes() {
        return genes;
    }

    public Gene getGene(int geneIndex){
        return this.getGenes()[geneIndex];
    }

    /**
     * Read from pgf file of gene annotation
     * @param infileS
     */
    @SuppressFBWarnings("DM_BOXED_PRIMITIVE_FOR_PARSING")
    private void readFile (String infileS) {
        try {
            BufferedReader br = IOTool.getReader(infileS);
            int geneNumber = Integer.parseInt(br.readLine().split("\t")[1]);
            genes = new Gene[geneNumber];
            String temp;
            for (int i = 0; i < geneNumber; i++) {
                temp = br.readLine();
                String[] tem = temp.split("\t");
                genes[i] = new Gene(tem[1], Integer.parseInt(tem[2]), Integer.parseInt(tem[3]), Integer.parseInt(tem[4]), Byte.parseByte(tem[5]), tem[6], tem[7]);
                tem = br.readLine().split("\t");
                int transcriptNumber = Integer.parseInt(tem[1]);
                genes[i].setLongestTranscriptIndex(Integer.parseInt(tem[2]));
                for (int j = 0; j < transcriptNumber; j++) {
                    temp = br.readLine();
                    tem = temp.split("\t");
                    int chr = Integer.parseInt(tem[2]);
                    Transcript t = new Transcript(tem[1], chr, Integer.parseInt(tem[3]), Integer.parseInt(tem[4]), Byte.parseByte(tem[5]));
                    tem = br.readLine().split("\t");
                    if (!tem[1].startsWith("NA")) {
                        tem = tem[1].split(";");
                        for (int k = 0; k < tem.length; k++) {
                            String[] te = tem[k].split(":");
                            t.add5UTR(chr, Integer.parseInt(te[0]), Integer.parseInt(te[1]));
                        }
                    }
                    tem = br.readLine().split("\t");
                    tem = tem[1].split(";");
                    for (int k = 0; k < tem.length; k++) {
                        String[] te = tem[k].split(":");
                        t.addCDS(chr, Integer.parseInt(te[0]), Integer.parseInt(te[1]));
                    }
                    br.readLine();
                    tem = br.readLine().split("\t");
                    if (!tem[1].startsWith("NA")) {
                        tem = tem[1].split(";");
                        for (int k = 0; k < tem.length; k++) {
                            String[] te = tem[k].split(":");
                            t.add3UTR(chr, Integer.parseInt(te[0]), Integer.parseInt(te[1]));
                        }
                    }
                    t.calculateIntron();
                    genes[i].addTranscript(t);
                }
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.sortGeneByGeneRange();
    }

    /**
     * Write all CDS of the longest transcripts to a file
     * @param genomef
     * @param outfileS
     */
    public void writeCDSSequence (FastaByte genomef, String outfileS) {
        genomef.sortByName();
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String title, chrseq, cdsSeq;
            StringBuilder sb;
            List<Range> cdsList;
            SequenceByte s;
            for (int i = 0; i < this.getGeneNumber(); i++) {
                sb=new StringBuilder(50);
                sb.append(this.getGeneChromosome(i)).append("_").append(this.getGeneStart(i)).append("_");
                sb.append(this.getGeneEnd(i)).append("_").append(this.getGeneName(i));
                title=sb.toString();
                int chrIndex = genomef.getIndexByName(String.valueOf(this.getGeneChromosome(i)));
                chrseq = genomef.getSeq(chrIndex);
                sb = new StringBuilder(2100);
                int longestTranscriptIndex = this.getLongestTranscriptIndex(i);
                cdsList = this.getCDSList(i, longestTranscriptIndex);
                for (int j = 0; j < cdsList.size(); j++) {
                    sb.append(chrseq.subSequence(cdsList.get(j).getRangeStart()-1, cdsList.get(j).getRangeEnd()-1));
                }
                cdsSeq = sb.toString();
                if (this.getTranscriptStrand(i, longestTranscriptIndex) == 0) {
                    s = new SequenceByte(cdsSeq);
                    cdsSeq = s.getReverseComplementarySeq();
                }
                sb=new StringBuilder(51);
                sb.append(">").append(title);
                bw.write(sb.toString());
                bw.newLine();
                bw.write(PStringUtils.getMultiplelineString(60, cdsSeq));
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
     * Write CDS sequences to different chromosome files
     * @param genomeFa_Dir
     * @param outDir
     */
    public void writeCDSSequencePerChr(String genomeFa_Dir, String outDir){
        long start=System.nanoTime();
        File[] genomeFa=IOUtils.listRecursiveFiles(new File(genomeFa_Dir));
        Predicate<File> p=File::isHidden;
        TIntArrayList p_1_42 = new TIntArrayList(IntStream.range(1, 43).toArray());
        Predicate<File> p1=f->p_1_42.contains(StringTool.getNumFromString(f.getName()));
        File[] fa=Arrays.stream(genomeFa).filter(p.negate()).filter(p1).toArray(File[]::new);
        String[] outNames=Arrays.stream(fa).map(File::getName).map(str->str.substring(0, 6)+".genes.fa").toArray(String[]::new);
        PGF[] chrPGF=this.getGeneOnAllChr();
        TIntArrayList chrs=this.getChrs();
//        if (fa.length!=chrs.size()){
//            System.out.println("error, check "+genomeFa_Dir+" and "+outDir);
//            System.exit(1);
//        }
        IntStream.range(0, fa.length).parallel().forEach(e->
                writeCDSSequencePerChr(chrs.get(e), chrPGF[e], fa[e].getAbsolutePath(), new File(outDir, outNames[e]).getAbsolutePath()));
        System.out.println(" pgf cds sequences per chromosome were write to "+outDir+" in "+Benchmark.getTimeSpanMinutes(start)+" minutes");
    }

    private void writeCDSSequencePerChr(int chr, PGF chrPGF, String chrFaFile, String outFile){
        long start=System.nanoTime();
        try {
            FastaByte chrFa=new FastaByte(chrFaFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String title, chrseq, cdsSeq;
            StringBuilder sb;
            List<Range> cdsList;
            SequenceByte s;
            for (int i = 0; i < chrPGF.getGeneNumber(); i++) {
                sb=new StringBuilder(100);
                int longestTranscriptIndex = chrPGF.getLongestTranscriptIndex(i);
                sb.append(chrPGF.getTranscriptName(i, longestTranscriptIndex));
                title=sb.toString();
                int chrIndex = chrFa.getIndexByName(String.valueOf(chr));
                chrseq = chrFa.getSeq(chrIndex);
                cdsList = chrPGF.getCDSList(i, longestTranscriptIndex);
                sb = new StringBuilder(2100);
                for (int j = 0; j < cdsList.size(); j++) {
                    sb.append(chrseq.subSequence(cdsList.get(j).getRangeStart()-1, cdsList.get(j).getRangeEnd()-1));
                }
                cdsSeq = sb.toString();
                if (chrPGF.getTranscriptStrand(i, longestTranscriptIndex)== 0) {
                    s = new SequenceByte(cdsSeq);
                    cdsSeq = s.getReverseComplementarySeq();
                }
                sb=new StringBuilder(101);
                sb.append(">").append(title);
                bw.write(sb.toString());
                bw.newLine();
                bw.write(PStringUtils.getMultiplelineString(60, cdsSeq));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println(chr+" completed in "+ Benchmark.getTimeSpanMinutes(start)+" minutes");
    }

    /**
     * Write genomic sequence of gene, from start to the end, no flip if the gene is in the minus direction
     * @param genomef
     * @param outfileS
     */
    public void writeGeneSequence (FastaByte genomef, String outfileS) {
        genomef.sortByName();
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < this.getGeneNumber(); i++) {
                String title = this.getGeneChromosome(i) +"_"+ this.getGeneStart(i) + "_" + this.getGeneEnd(i) + "_" + this.getGeneName(i);
                int chrIndex = genomef.getIndexByName(String.valueOf(this.getGeneChromosome(i)));
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
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("GeneNumber\t"+ this.getGeneNumber());
            bw.newLine();
            for (int i = 0; i < this.getGeneNumber(); i++) {
                StringBuilder sb = new StringBuilder();
                sb.append("Gene\t").append(this.getGeneName(i)).append("\t").append(this.getGeneChromosome(i)).append("\t").append(this.getGeneStart(i)).append("\t").append(this.getGeneEnd(i)).append("\t").append(this.getGeneStrand(i));
                sb.append("\t").append(this.getGeneBiotype(i)).append("\t").append(this.getGeneDescription(i));
                bw.write(sb.toString());
                bw.newLine();
                sb = new StringBuilder("TranscriptNumber\t");
                sb.append(this.getTranscriptNumber(i)).append("\t").append(genes[i].longestTranscriptIndex);
                bw.write(sb.toString());
                bw.newLine();
                for (int j = 0; j < this.getTranscriptNumber(i); j++) {
                    sb = new StringBuilder();
                    sb.append("Transcript\t").append(this.getTranscriptName(i, j)).append("\t").append(this.getTranscriptChromosome(i, j)).append("\t").append(this.getTranscriptStart(i, j)).append("\t").append(this.getTranscriptEnd(i, j)).append("\t").append(this.getTranscriptStrand(i,j));
                    bw.write(sb.toString());
                    bw.newLine();
                    sb = new StringBuilder();
                    sb.append("5'UTR\t");
                    if (this.isThere5UTR(i, j)) {
                        sb.append(this.get5UTRPositionString(i, j));
                    }
                    else {
                        sb.append("NA");
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                    bw.write("CDS\t"+this.getCDSPositionString(i, j));
                    bw.newLine();
                    if (this.getIntronPositionString(i, j).equals("")) {
                        bw.write("Intron\t"+"NA");
                    }
                    else {
                        bw.write("Intron\t"+this.getIntronPositionString(i, j));
                    }
                    bw.newLine();
                    sb = new StringBuilder();
                    sb.append("3'UTR\t");
                    if (this.isThere3UTR(i, j)) {
                        sb.append(this.get3UTRPositionString(i, j));
                    }
                    else {
                        sb.append("NA");
                    }
                    bw.write(sb.toString());
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
     * Removes the genes that satisfy the given predicate
     * @param predicate 针对基因进行过滤的函数式接口
     * @return true, if any gene was removed
     */
    public boolean removeIf(Predicate<Gene> predicate){
        List<Gene> genes=CollectionTool.changeToList(this.genes);
        boolean res=genes.removeIf(predicate);
        this.genes=genes.toArray(new Gene[genes.size()]);
        return res;
    }

    public TIntArrayList getCDSPosOnGene(int geneIndex, int transcriptIndex){
        List<Range> cdsList=this.getCDSList(geneIndex, transcriptIndex);
        TIntArrayList cdsPosListOnGene=new TIntArrayList();
        for (int i = 0; i < cdsList.size(); i++) {
            int start=cdsList.get(i).start;
            int end=cdsList.get(i).end;
            cdsPosListOnGene.add(IntStream.range(start, end).toArray());
        }
        cdsPosListOnGene.sort();
        return cdsPosListOnGene;
    }

    public int getCDSLen(int geneIndex, int transcriptIndex){
        List<Range> cdsList=this.getCDSList(geneIndex, transcriptIndex);
        int len=0;
        for (int i = 0; i < cdsList.size(); i++) {
            len+=cdsList.get(i).getRangeSize();
        }
        return len;
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

    public Map<Integer, TIntArrayList> getOverlappedGeneIndex(){
        RangeInterface interSection;
        Map<Integer, TIntArrayList> geneIndexOverlappedGeneindexMap=new HashMap<>();
        List<RangeInterface>  interSectionList;
        TIntArrayList overlappedGeneIndex;
        for (int i = 0; i < genes.length-1; i++) {
            interSectionList=new ArrayList<>();
            overlappedGeneIndex=new TIntArrayList();
            for (int j = i+1; j < genes.length; j++) {
                interSection= genes[i].geneRange.getIntersection(genes[j].geneRange);
                if (interSection==null) break;
                interSectionList.add(interSection);
                overlappedGeneIndex.add(j);
            }
            if (interSectionList.size()==0) continue;
//            System.out.println("overlapped gene "+genes[i].geneName+" overlapped count "+interSectionList.size());
            geneIndexOverlappedGeneindexMap.put(i, overlappedGeneIndex);
        }
        return geneIndexOverlappedGeneindexMap;
    }

    public void writeOverlappedGene(String outFile){
        Map<Integer, TIntArrayList> geneIndexOverlappedGeneindexMap = this.getOverlappedGeneIndex();
        List<Integer> geneIndexList=new ArrayList<>(geneIndexOverlappedGeneindexMap.keySet());
        Collections.sort(geneIndexList);
        try (BufferedWriter bw = IOUtils.getTextWriter(outFile)) {
            bw.write("GeneNameIndex\tOverlappedGeneNameIndex");
            bw.newLine();
            StringBuilder sb;
            int geneIndex;
            TIntArrayList values;
            for (int i = 0; i < geneIndexList.size(); i++) {
                sb=new StringBuilder();
                geneIndex=geneIndexList.get(i);
                values=geneIndexOverlappedGeneindexMap.get(geneIndex);
                sb.append(genes[geneIndex].geneName).append(",").append(geneIndex).append("\t");
                for (int j = 0; j < values.size(); j++) {
                    sb.append(genes[values.get(j)].geneName).append(",").append(values.get(j)).append(";");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public int getCDSIndex (int geneIndex, int transcriptIndex, int chr, int pos) {
        return this.genes[geneIndex].ts.get(transcriptIndex).getCDSIndex(chr, pos);
    }

    public int get5UTRIndex (int geneIndex, int transcriptIndex, int chr, int pos) {
        return this.genes[geneIndex].ts.get(transcriptIndex).get5UTRIndex(chr, pos);
    }

    public int get3UTRIndex (int geneIndex, int transcriptIndex, int chr, int pos) {
        return this.genes[geneIndex].ts.get(transcriptIndex).get3UTRIndex(chr, pos);
    }

    /**
     *
     * @param geneIndex
     * @param chr
     * @param pos
     * @return negative value if the position is not in the range of any gene's longest transcript
     */
    public boolean isWithinThisGeneExon(int geneIndex, int chr, int pos){
        int longestTranscriptIndex=this.getLongestTranscriptIndex(geneIndex);
        int cdsIndex=this.getCDSIndex(geneIndex, longestTranscriptIndex, chr, pos);
        int utr5Index=this.get5UTRIndex(geneIndex, longestTranscriptIndex, chr, pos);
        int utr3Index=this.get3UTRIndex(geneIndex, longestTranscriptIndex, chr, pos);
        if (cdsIndex < 0 && utr5Index < 0 && utr3Index < 0){
            return false;
        }else {
            return true;
        }
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
        Gene query  = new Gene ("Query", chromosome, Integer.MIN_VALUE, Integer.MIN_VALUE, Byte.MIN_VALUE, "", "");
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
        Gene query  = new Gene ("Query", chromosome+1, Integer.MIN_VALUE, Integer.MIN_VALUE, Byte.MIN_VALUE, "", "");
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
    public String getGeneBiotype (int index) {
        return genes[index].biotype;
    }

    /**
     * Return the description of the gene
     * @param index
     * @return
     */
    public String getGeneDescription (int index) {
        return genes[index].description;
    }

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
     * 获取指定转录本的exon总长度
     * @param i
     * @param j
     * @return
     */
    public int getTranscriptExonLen(int i, int j){
        return genes[i].ts.get(j).getExonLen();
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
     *
     * @param infileS
     */
    public void readFromWheatGFF (String infileS) {
        try {
            BufferedReader br;
            if (infileS.endsWith("gz")) br = IOUtils.getTextGzipReader(infileS);
            else br = IOUtils.getTextReader(infileS);
            String temp = br.readLine();
            ArrayList<String> infoList = new ArrayList();
            ArrayList<String> geneList = new ArrayList();
            String[] tem = null;
            while ((temp = br.readLine()) != null) {
                char s = temp.charAt(0);
                if ((int)s < 48 || (int)s > 57) continue;
                List<String> tList= PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                if (tem[2].startsWith("exon")) continue;
                if (tem[2].startsWith("chromosome")) continue;
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
                    String biotype = "NA";
                    String description = "NA";
                    for (int j = 1; j < te.length; j++) {
                        if (te[j].startsWith("biotype")) {
                            biotype = te[j].replaceFirst("biotype=", "");
                        }
                        else if (te[j].startsWith("description")) {
                            description = te[j].replaceFirst("description=", "");
                        }
                    }
                    genes[index] = new Gene (query, Integer.parseInt(tem[0]), Integer.parseInt(tem[3]), Integer.parseInt(tem[4])+1, (byte)(tem[6].equals("+")? 1:0), biotype, description);
                }
            }
            for (int i = 0; i < info.length; i++) {
                tem = info[i].split("\t");
                if (tem[2].startsWith("mRNA")) {
                    String[] te = tem[8].split(";");
                    String geneName = te[1].split("=")[1];
                    int geneIndex = Arrays.binarySearch(geneNames, geneName);
                    Transcript t = new Transcript (te[0].split("=")[1], Integer.parseInt(tem[0]), Integer.parseInt(tem[3]), Integer.parseInt(tem[4])+1, (byte)(tem[6].equals("+")? 1:0));
                    genes[geneIndex].addTranscript(t);
                }
            }
            for (int i = 0; i < genes.length; i++) genes[i].sortTranscriptsByName();
            for (int i = 0; i < info.length; i++) {
                tem = info[i].split("\t");
                if (tem[2].startsWith("CDS")) {
                    String[] te = tem[8].split(";");
                    String transcriptName = te[1].split("=")[1];
                    String geneName;
                    geneName = transcriptName.split("\\.")[0];
                    int geneIndex = Arrays.binarySearch(geneNames, geneName);
                    int transcriptIndex = genes[geneIndex].getTranscriptIndex(transcriptName);
                    genes[geneIndex].ts.get(transcriptIndex).addCDS(Integer.parseInt(tem[0]), Integer.parseInt(tem[3]), Integer.parseInt(tem[4])+1);
                }
                else if (tem[2].startsWith("five_prime_UTR")) {
                    String[] te = tem[8].split(";");
                    String transcriptName = te[1].split("=")[1];
                    String geneName;
                    geneName = transcriptName.split("\\.")[0];
                    int geneIndex = Arrays.binarySearch(geneNames, geneName);
                    int transcriptIndex = genes[geneIndex].getTranscriptIndex(transcriptName);
                    genes[geneIndex].ts.get(transcriptIndex).add5UTR(Integer.parseInt(tem[0]), Integer.parseInt(tem[3]), Integer.parseInt(tem[4])+1);
                }
                else if (tem[2].startsWith("three_prime_UTR")) {
                    String[] te = tem[8].split(";");
                    String transcriptName = te[1].split("=")[1];
                    String geneName;
                    geneName = transcriptName.split("\\.")[0];
                    int geneIndex = Arrays.binarySearch(geneNames, geneName);
                    int transcriptIndex = genes[geneIndex].getTranscriptIndex(transcriptName);
                    genes[geneIndex].ts.get(transcriptIndex).add3UTR(Integer.parseInt(tem[0]), Integer.parseInt(tem[3]), Integer.parseInt(tem[4])+1);
                }
            }
            for (int i = 0; i < this.genes.length; i++) {
                for (int j = 0; j < genes[i].ts.size(); j++) {
                    genes[i].ts.get(j).sort5UTRByPosition();
                    genes[i].ts.get(j).sortCDSByPosition();
                    genes[i].ts.get(j).sort3UTRByPosition();
                    genes[i].ts.get(j).calculateIntron();
                }
                genes[i].calculateLongestTranscriptIndex();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.sortGeneByGeneRange();
    }

    /**
     * Read from Maize GFF file (ftp://ftp.ensemblgenomes.org/pub/plants/release-38/gff3/zea_mays)
     * @param infileS
     */
    public void readFromMaizeGFF (String infileS) {
        try {
            BufferedReader br;
            if (infileS.endsWith("gz")) br = IOUtils.getTextGzipReader(infileS);
            else br = IOUtils.getTextReader(infileS);
            String temp = br.readLine();
            ArrayList<String> infoList = new ArrayList();
            ArrayList<String> geneList = new ArrayList();
            String[] tem = null;
            while ((temp = br.readLine()) != null) {
                char s = temp.charAt(0);
                if ((int)s < 48 || (int)s > 57) continue;
                List<String> tList= PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                if (tem[2].startsWith("exon")) continue;
                if (tem[2].startsWith("chromosome")) continue;
                if (tem[2].startsWith("gene")) {
                    String[] te = tem[8].split(";");
                    geneList.add(te[0].split(":")[1]);
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
                    String query = te[0].split(":")[1];
                    int index = Arrays.binarySearch(geneNames, query);
                    String biotype = "NA";
                    String description = "NA";
                    for (int j = 1; j < te.length; j++) {
                        if (te[j].startsWith("biotype")) {
                            biotype = te[j].replaceFirst("biotype=", "");
                        }
                        else if (te[j].startsWith("description")) {
                            description = te[j].replaceFirst("description=", "");
                        }
                    }
                    genes[index] = new Gene (query, Integer.parseInt(tem[0]), Integer.parseInt(tem[3]), Integer.parseInt(tem[4])+1, (byte)(tem[6].equals("+")? 1:0), biotype, description);
                }
            }
            for (int i = 0; i < info.length; i++) {
                tem = info[i].split("\t");
                if (tem[2].startsWith("mRNA")) {
                    String[] te = tem[8].split(";");
                    String geneName = te[1].split(":")[1];
                    int geneIndex = Arrays.binarySearch(geneNames, geneName);
                    Transcript t = new Transcript (te[0].split(":")[1], Integer.parseInt(tem[0]), Integer.parseInt(tem[3]), Integer.parseInt(tem[4])+1, (byte)(tem[6].equals("+")? 1:0));
                    genes[geneIndex].addTranscript(t);
                }
            }
            for (int i = 0; i < genes.length; i++) genes[i].sortTranscriptsByName();
            for (int i = 0; i < info.length; i++) {
                tem = info[i].split("\t");
                if (tem[2].startsWith("CDS")) {
                    String[] te = tem[8].split(";");
                    String transcriptName = te[1].split(":")[1];
                    String geneName;
                    geneName = transcriptName.split("_")[0];
                    int geneIndex = Arrays.binarySearch(geneNames, geneName);
                    int transcriptIndex = genes[geneIndex].getTranscriptIndex(transcriptName);
                    genes[geneIndex].ts.get(transcriptIndex).addCDS(Integer.parseInt(tem[0]), Integer.parseInt(tem[3]), Integer.parseInt(tem[4])+1);
                }
                else if (tem[2].startsWith("five_prime_UTR")) {
                    String[] te = tem[8].split(":");
                    String transcriptName = te[1];
                    String geneName;
                    geneName = transcriptName.split("_")[0];
                    int geneIndex = Arrays.binarySearch(geneNames, geneName);
                    int transcriptIndex = genes[geneIndex].getTranscriptIndex(transcriptName);
                    genes[geneIndex].ts.get(transcriptIndex).add5UTR(Integer.parseInt(tem[0]), Integer.parseInt(tem[3]), Integer.parseInt(tem[4])+1);
                }
                else if (tem[2].startsWith("three_prime_UTR")) {
                    String[] te = tem[8].split(":");
                    String transcriptName = te[1];
                    String geneName;
                    geneName = transcriptName.split("_")[0];
                    int geneIndex = Arrays.binarySearch(geneNames, geneName);
                    int transcriptIndex = genes[geneIndex].getTranscriptIndex(transcriptName);
                    genes[geneIndex].ts.get(transcriptIndex).add3UTR(Integer.parseInt(tem[0]), Integer.parseInt(tem[3]), Integer.parseInt(tem[4])+1);
                }
            }
            for (int i = 0; i < this.genes.length; i++) {
                for (int j = 0; j < genes[i].ts.size(); j++) {
                    genes[i].ts.get(j).sort5UTRByPosition();
                    genes[i].ts.get(j).sortCDSByPosition();
                    genes[i].ts.get(j).sort3UTRByPosition();
                    genes[i].ts.get(j).calculateIntron();
                }
                genes[i].calculateLongestTranscriptIndex();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.sortGeneByGeneRange();
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
     * Sort genes by their gene range
     */
    public void sortGeneByGeneRange() {
        this.sortType = 0;
        Arrays.sort(genes);
    }

    /**
     *
     * @return all chromosomes in this object
     */
    public TIntArrayList getChrs(){
        TIntHashSet chrSet=new TIntHashSet();
        for (int i = 0; i < this.getGeneNumber(); i++) {
            chrSet.add(genes[i].geneRange.chr);
        }
        TIntArrayList chrList=new TIntArrayList(chrSet);
        chrList.sort();
        return chrList;
    }

    /**
     *
     * @param chr
     * @return all genes on target chromosome
     */
    public PGF getGeneOnChr(int chr){
        this.sortGeneByGeneRange();
        int chrIndexA=this.getStartIndexOfChromosome(chr);
        int chrIndexB=this.getEndIndexOfChromosome(chr);
        Gene[] genes=ArrayUtils.subarray(this.genes, chrIndexA, chrIndexB);
        return new PGF(genes);
    }

    /**
     *
     * @param chr
     * @return all genes on all target chromosome
     */
    public PGF getGeneOnChr(int[] chr){
        this.sortGeneByGeneRange();
        List<Gene> list=new ArrayList<>();
        Gene[] genes;
        for (int i = 0; i < chr.length; i++) {
            int chrIndexA=this.getStartIndexOfChromosome(chr[i]);
            int chrIndexB=this.getEndIndexOfChromosome(chr[i]);
            genes=ArrayUtils.subarray(this.genes, chrIndexA, chrIndexB);
            for (int j = 0; j < genes.length; j++) {
                list.add(genes[j]);
            }
        }
        Gene[] resGene=new Gene[list.size()];
        for (int i = 0; i < list.size(); i++) {
            resGene[i]=list.get(i);
        }
        return new PGF(resGene);
    }

    /**
     *
     * @param chrs
     * @return all genesName on all target chromosome sorted by geneRange if chrs is sorted
     */
    public List<String> getGeneNameOnChr(int[] chrs){
        this.sortGeneByGeneRange();
        List<String> list=new ArrayList<>();
        Gene[] genes;
        for (int i = 0; i < chrs.length; i++) {
            int chrIndexA=this.getStartIndexOfChromosome(chrs[i]);
            int chrIndexB=this.getEndIndexOfChromosome(chrs[i]);
            genes=ArrayUtils.subarray(this.genes, chrIndexA, chrIndexB);
            for (int j = 0; j < genes.length; j++) {
                list.add(genes[j].getGeneName());
            }
        }
        return list;
    }

    public List<Gene> getGeneList(int[] chr){
        this.sortGeneByGeneRange();
        List<Gene> list=new ArrayList<>();
        Gene[] genes;
        for (int i = 0; i < chr.length; i++) {
            int chrIndexA=this.getStartIndexOfChromosome(chr[i]);
            int chrIndexB=this.getEndIndexOfChromosome(chr[i]);
            genes=ArrayUtils.subarray(this.genes, chrIndexA, chrIndexB);
            for (int j = 0; j < genes.length; j++) {
                list.add(genes[j]);
            }
        }
        Collections.sort(list);
        return list;
    }

    /**
     *
     * @return genes on each chromosome
     */
    public PGF[] getGeneOnAllChr(){
        TIntArrayList chrList=this.getChrs();
        PGF[] genes=new PGF[chrList.size()];
        for (int i = 0; i < chrList.size(); i++) {
            genes[i]=this.getGeneOnChr(chrList.get(i));
        }
        return genes;
    }

    /**
     * 返回指定范围内的基因总数
     * @param chromosome
     * @param posStart inclusive
     * @param posEnd exclusive
     * @return
     */
    public int getGeneNum(String chromosome, int posStart, int posEnd){
        int chrStart= RefV1Utils.getChrID(chromosome, posStart);
        int chrEnd=RefV1Utils.getChrID(chromosome, posEnd);
        int posOnChrIDStart=RefV1Utils.getPosOnChrID(chromosome, posStart);
        int posOnChrIDEnd=RefV1Utils.getPosOnChrID(chromosome, posEnd);
        this.sortGeneByGeneRange();
        Gene queryStart=new Gene(chrStart, posOnChrIDStart, posOnChrIDStart+1);
        Gene queryEnd=new Gene(chrEnd, posOnChrIDEnd, posOnChrIDEnd+1);
        int startHit = Arrays.binarySearch(this.genes, queryStart);
        int endHit = Arrays.binarySearch(this.genes, queryEnd);
        int indexStart, indexEnd;
        if (startHit < -1){
            indexStart = -startHit-2;
        }else if (startHit > -1){
            indexStart = startHit;
        }else {
            indexStart = 0;
        }
        if (endHit < -1){
            indexEnd = -endHit-2;
        }else if (endHit > -1){
            indexEnd = endHit;
        }else {
            indexEnd = 0;
        }
        return indexEnd-indexStart;
    }

    public class Gene implements Comparable<Gene> {
        String geneName = null;
        Range geneRange = null;
        byte strand = Byte.MIN_VALUE;
        String biotype = null;
        String description = null;
        ArrayList<Transcript> ts = new ArrayList();
        int longestTranscriptIndex = -1;

        public Gene (String geneName, int chr, int start, int end, byte strand, String biotype, String discription) {
            this.geneName = geneName;
            this.strand = strand;
            geneRange = new Range(chr, start, end);
            this.biotype = biotype;
            this.description = discription;
        }

        public Gene (String geneName) {
            this.geneName = geneName;
        }

        public Gene (int chr, int start, int end) {
            geneRange = new Range(chr, start, end);
        }

        public String getGeneName() {
            return geneName;
        }

        public Range getGeneRange() {
            return geneRange;
        }

        public byte getStrand() {
            return strand;
        }

        public String getBiotype() {
            return biotype;
        }

        public String getDescription() {
            return description;
        }

        public ArrayList<Transcript> getTs() {
            return ts;
        }

        public int getLongestTranscriptIndex () {
            return this.longestTranscriptIndex;
        }

        public void addTranscript (Transcript t) {
            ts.add(t);
        }

        public int getTranscriptNumber () {
            return ts.size();
        }

        public int getTranscriptIndex (String transcriptName) {
            return Collections.binarySearch(ts, new Transcript(transcriptName));
        }

        public int getLongestTranscriptCDSLen(){
            int longestTranscriptIndex=this.getLongestTranscriptIndex();
            return this.getTs().get(longestTranscriptIndex).getCDSLen();
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

    public class Transcript implements Comparable<Transcript> {
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

        public byte getStrand() {
            return strand;
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

        public int get5UTRIndex (int chr, int pos) {
            Range query = new Range(chr, pos, pos+1);
            int hit = Collections.binarySearch(utr5List, query);
            int index = hit;
            if (index < -1) {
                index = -index-2;
                if (this.isWithinThis5UTR(index, chr, pos)) return index;
            }
            return hit;
        }

        public boolean isWithinThis5UTR (int cdsIndex, int chr, int pos) {
            if (utr5List.get(cdsIndex).chr != chr) return false;
            if (pos < utr5List.get(cdsIndex).start) return false;
            if (pos >= utr5List.get(cdsIndex).end) return false;
            return true;
        }

        public int get3UTRIndex (int chr, int pos) {
            Range query = new Range(chr, pos, pos+1);
            int hit = Collections.binarySearch(utr3List, query);
            int index = hit;
            if (index < -1) {
                index = -index-2;
                if (this.isWithinThis3UTR(index, chr, pos)) return index;
            }
            return hit;
        }

        public boolean isWithinThis3UTR (int cdsIndex, int chr, int pos) {
            if (utr3List.get(cdsIndex).chr != chr) return false;
            if (pos < utr3List.get(cdsIndex).start) return false;
            if (pos >= utr3List.get(cdsIndex).end) return false;
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

        public int getExonLen(){
            int len=0;
            for (int i = 0; i < utr5List.size(); i++) {
                len+=utr5List.get(i).getRangeSize();
            }
            for (int i = 0; i < utr3List.size(); i++) {
                len+=utr3List.get(i).getRangeSize();
            }
            for (int i = 0; i < cdsList.size(); i++) {
                len+=cdsList.get(i).getRangeSize();
            }
            return len;
        }

        public int getCDSLen(){
            int len=0;
            for (int i = 0; i < cdsList.size(); i++) {
                len+=cdsList.get(i).getRangeSize();
            }
            return len;
        }

        public List<Range> getCdsList(){
            return this.cdsList;
        }

        @Override
        public int compareTo(Transcript t) {
            return transcriptName.compareTo(t.transcriptName);
        }
    }
}

