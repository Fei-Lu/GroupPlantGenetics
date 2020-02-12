/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.annotation;

import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TFloatArrayList;
import java.awt.PageAttributes;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
//import utils.FStringUtils;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author xuebozhao
 */
public class GeneFeatureReadGff3 {
    Gene[] genes;
    int sortType = 0;
    public GeneFeatureReadGff3 () {}
    
    public GeneFeatureReadGff3 (String infileS,String outfileS) {
//        this.readFile(infileS);
      this.readFromMaizeGFF(infileS);     
      this.writeFile(outfileS);
    }

    /**
     * Read from txt file of gene annotation
     * @param infileS 
     */
    private void readFile (String infileS) {
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            int geneNumber = Integer.valueOf(br.readLine().split("\t")[1]);
            genes = new Gene[geneNumber];
            String temp;
            for (int i = 0; i < geneNumber; i++) {
                temp = br.readLine();
                String[] tem = temp.split("\t");
                genes[i] = new Gene(tem[1], Integer.valueOf(tem[2]), Integer.valueOf(tem[3]), Integer.valueOf(tem[4]), Byte.valueOf(tem[5]));
                int transcriptNumber = Integer.valueOf(br.readLine().split("\t")[1]);
                for (int j = 0; j < transcriptNumber; j++) {
                    temp = br.readLine();
                    tem = temp.split("\t");
                    int chr = Integer.valueOf(tem[2]);
                    Transcript t = new Transcript(tem[1], chr, Integer.valueOf(tem[3]), Integer.valueOf(tem[4]), Byte.valueOf(tem[5]));
                    tem = br.readLine().split("\t");
                    if (!tem[1].startsWith("null")) {
                        tem = tem[1].split(";");
                        for (int k = 0; k < tem.length; k++) {
                            String[] te = tem[k].split(":");
                            t.add5UTR(chr, Integer.valueOf(te[0]), Integer.valueOf(te[1]));
                        }
                    }
                    tem = br.readLine().split("\t");
                    tem = tem[1].split(";");
                    for (int k = 0; k < tem.length; k++) {
                        String[] te = tem[k].split(":");
                        t.addCDS(chr, Integer.valueOf(te[0]), Integer.valueOf(te[1]));
                    }
                    br.readLine();
                    tem = br.readLine().split("\t");
                    if (!tem[1].startsWith("null")) {
                        tem = tem[1].split(";");
                        for (int k = 0; k < tem.length; k++) {
                            String[] te = tem[k].split(":");
                            t.add3UTR(chr, Integer.valueOf(te[0]), Integer.valueOf(te[1]));
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
        this.sortGeneByStartPosition();
    }
    
    public void writeCDSSequence (Fasta genomef, String outfileS) {
        genomef.sortRecordByName();
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < this.getGeneNumber(); i++) {
                String title = String.valueOf(this.getGeneChromosome(i))+"_"+String.valueOf(this.getGeneStart(i)+"_"+String.valueOf(this.getGeneEnd(i))+"_"+String.valueOf(this.getGeneName(i)));
                int chrIndex = genomef.getIndex(String.valueOf(this.getGeneChromosome(i)));
                String chrseq = genomef.getSeq(chrIndex);
                StringBuilder sb = new StringBuilder();
                ArrayList<Range> cdsList = this.getCDSList(i, 0);
                for (int j = 0; j < cdsList.size(); j++) {
                    sb.append(chrseq.subSequence(cdsList.get(j).getRangeStart()-1, cdsList.get(j).getRangeEnd()-1));
                }
                String cdsSeq = sb.toString();
                if (this.getTranscriptStrand(i, 0) == 0) {
                    Sequence s = new Sequence(cdsSeq);
                    cdsSeq = s.getReverseComplementarySeq();
                }
                bw.write(">"+title);
                bw.newLine();
                bw.write(FStringUtils.getMultiplelineString(60, cdsSeq));
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
     * Write genomic sequence of gene, from start to the end, no flip if the gene is in the minus direction
     * @param genomef
     * @param outfileS 
     */
    public void writeGeneSequence (Fasta genomef, String outfileS) {
        genomef.sortRecordByName();
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < this.getGeneNumber(); i++) {
                String title = String.valueOf(this.getGeneChromosome(i))+"_"+String.valueOf(this.getGeneStart(i)+"_"+String.valueOf(this.getGeneEnd(i))+"_"+String.valueOf(this.getGeneName(i)));
                int chrIndex = genomef.getIndex(String.valueOf(this.getGeneChromosome(i)));
                String chrseq = genomef.getSeq(chrIndex);
                String geneSeq = chrseq.substring(this.getGeneStart(i)-1, this.getGeneEnd(i)-1);
                String[] geneSeqs = FStringUtils.getMultilineString(60, geneSeq);
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
     * Write txt file of gene annotation
     * @param outfileS 
     */
    public void writeFile (String outfileS) {
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("GeneNumber\t"+String.valueOf(this.getGeneNumber()));
            bw.newLine();
            for (int i = 0; i < this.getGeneNumber(); i++) {
                StringBuilder sb = new StringBuilder();
                sb.append("Gene\t").append(this.getGeneName(i)).append("\t").append(this.getGeneChromosome(i)).append("\t").append(this.getGeneStart(i)).append("\t").append(this.getGeneEnd(i)).append("\t").append(this.getGeneStrand(i));
                bw.write(sb.toString());
                bw.newLine();
                bw.write("TranscriptNumber\t"+String.valueOf(this.getTranscriptNumber(i)));
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
                        sb.append("null");
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                    bw.write("CDS\t"+this.getCDSPositionString(i, j));
                    bw.newLine();
                    if (this.getIntronPositionString(i, j).equals("")) {
                        bw.write("Intron\t"+"null");
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
                        sb.append("null");
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
    
    public int getCDSIndex (int geneIndex, int transcriptIndex, int chr, int pos) {
        return this.genes[geneIndex].ts.get(transcriptIndex).getCDSIndex(chr, pos);
    }
    
    /**
     * Return gene index from of a position
     * @param chr
     * @param pos
     * @return 
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
     * Return index of a gene
     * @param geneName
     * @return 
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
    public ArrayList<Range> getCDSList (int i, int j) {
        return genes[i].ts.get(j).cdsList;
    }
    
    /**
     * Return a range list of 5UTR, return an empty list if there is no 5UTR
     * @param i
     * @param j
     * @return 
     */
    public ArrayList<Range> get5UTRList (int i, int j) {
        return genes[i].ts.get(j).utr5List;
    }
    
    /**
     * Return a range list of 3UTR, return an empty list if there is no 3UTR
     * @param i
     * @param j
     * @return 
     */
    public ArrayList<Range> get3UTRList (int i, int j) {
        return genes[i].ts.get(j).utr3List;
    }
    
    /**
     * Return intron Range list, return an empty list if there is no intron
     * @param i
     * @param j
     * @return 
     */
    public ArrayList<Range> getIntronList (int i, int j) {
        return genes[i].ts.get(j).intronList;
    }
    
    /**
     * Return a string of CDS positions for output, return empty when there is no CDS
     * @param i
     * @param j
     * @return 
     */
    public String getCDSPositionString (int i, int j) {
        return Ranges.getRangePositionString(this.getCDSList(i, j));
    }
    
    /**
     * Return a string of intron positions, return empty when there is no intron
     * @param i
     * @param j
     * @return 
     */
    public String getIntronPositionString (int i, int j) {
        return Ranges.getRangePositionString(this.getIntronList(i, j));
    }
    
    /**
     * Return a string of 5UTR positions, return empty when there is no 5UTR
     * @param i
     * @param j
     * @return 
     */
    public String get5UTRPositionString (int i, int j) {
        return Ranges.getRangePositionString(this.get5UTRList(i, j));
    }
    
    /**
     * Return a string of 3UTR positions, return empty when there is no 3UTR
     * @param i
     * @param j
     * @return 
     */
    public String get3UTRPositionString (int i, int j) {
        return Ranges.getRangePositionString(this.get3UTRList(i, j));
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
     * Read from cassava GFF
     * @param infileS 
     */
    public void readFromCassavaGFF (String infileS) {
        try {
            BufferedReader br;
            if (infileS.endsWith("gz")) br = IOUtils.getTextGzipReader(infileS);
            else br = IOUtils.getTextReader(infileS);
            String temp = br.readLine();
            ArrayList<String> infoList = new ArrayList();
            ArrayList<String> geneList = new ArrayList();
            while ((temp = br.readLine()) != null) {
                char s = temp.charAt(0);
                if ((int)s < 48 || (int)s > 57) continue;
                String[] tem = temp.split("\t");
                if (tem[1].equals("JGI_gene_alt")) continue;
                if (tem[1].equals("JGI_gene")) continue;
                if (tem[2].startsWith("repeat")) continue;
                if (tem[2].startsWith("exon")) continue;
                if (tem[2].startsWith("gene")) {
                    String[] te = tem[8].split("Name=");
                    geneList.add(te[1]);
                }
                infoList.add(temp);
            }
            String[] geneNames = geneList.toArray(new String[geneList.size()]);
            Arrays.sort(geneNames);
            genes = new Gene[geneNames.length];
            String[] info = infoList.toArray(new String[infoList.size()]);
            for (int i = 0; i < info.length; i++) {
                String[] tem = info[i].split("\t");
                if (tem[2].startsWith("gene")) {
                    String[] te = tem[8].split("Name=");
                    String query = te[1];
                    int index = Arrays.binarySearch(geneNames, query);
                    genes[index] = new Gene (query, Integer.valueOf(tem[0]), Integer.valueOf(tem[3]), Integer.valueOf(tem[4])+1, (byte)(tem[6].equals("+")? 1:0));
                }
            }
            HashMap<String, Integer> transcriptGeneMap = new HashMap();
            for (int i = 0; i < info.length; i++) {
                String[] tem = info[i].split("\t");
                if (tem[2].startsWith("mRNA")) {
                    String[] te = tem[8].split(";");
                    String geneName = te[3].replaceFirst("Parent=", "");
                    int geneIndex = Arrays.binarySearch(geneNames, geneName);
                    String transcpritName = te[0].replaceFirst("ID=", "");
                    transcriptGeneMap.put(transcpritName, geneIndex);
                    Transcript t = new Transcript (transcpritName, Integer.valueOf(tem[0]), Integer.valueOf(tem[3]), Integer.valueOf(tem[4])+1, (byte)(tem[6].equals("+")? 1:0));
                    genes[geneIndex].addTranscript(t);
                }
            }
            for (int i = 0; i < genes.length; i++) genes[i].sortTranscriptsByName();
            for (int i = 0; i < info.length; i++) {
                String[] tem = info[i].split("\t");
                if (tem[2].startsWith("CDS")) {
                    String[] te = tem[8].split(";");
                    String transcriptName = te[1].replaceFirst("Parent=", "");
                    int geneIndex = transcriptGeneMap.get(transcriptName);
                    int transcriptIndex = genes[geneIndex].getTranscriptIndex(transcriptName);
                    genes[geneIndex].ts.get(transcriptIndex).addCDS(Integer.valueOf(tem[0]), Integer.valueOf(tem[3]), Integer.valueOf(tem[4])+1);
                }
                else if (tem[2].startsWith("five_prime_UTR")) {
                    String[] te = tem[8].split(";");
                    String transcriptName = te[1].replaceFirst("Parent=", "");
                    
                    int geneIndex = transcriptGeneMap.get(transcriptName);
                    int transcriptIndex = genes[geneIndex].getTranscriptIndex(transcriptName);
                    genes[geneIndex].ts.get(transcriptIndex).add5UTR(Integer.valueOf(tem[0]), Integer.valueOf(tem[3]), Integer.valueOf(tem[4])+1);
                }
                else if (tem[2].startsWith("three_prime_UTR")) {
                    String[] te = tem[8].split(";");
                    String transcriptName = te[1].replaceFirst("Parent=", "");
                    int geneIndex = transcriptGeneMap.get(transcriptName);
                    int transcriptIndex = genes[geneIndex].getTranscriptIndex(transcriptName);
                    genes[geneIndex].ts.get(transcriptIndex).add3UTR(Integer.valueOf(tem[0]), Integer.valueOf(tem[3]), Integer.valueOf(tem[4])+1);
                }
            }
            for (int i = 0; i < this.genes.length; i++) {
                for (int j = 0; j < genes[i].ts.size(); j++) {
                    genes[i].ts.get(j).sort5UTRByPosition();
                    genes[i].ts.get(j).sortCDSByPosition();
                    genes[i].ts.get(j).sort3UTRByPosition();
                    genes[i].ts.get(j).calculateIntron();
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.sortGeneByStartPosition();
    }
    
    /**
     * Read from Maize GFF
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
            while ((temp = br.readLine()) != null) {
                char s = temp.charAt(0);
                if ((int)s < 48 || (int)s > 57) continue;
                String[] tem = temp.split("\t");
                if (tem[2].startsWith("chromosome")) continue;
                if (tem[2].startsWith("repeat")) continue;
                if (tem[2].startsWith("exon")) continue;
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
                String[] tem = info[i].split("\t");
                if (tem[2].startsWith("gene")) {
                    String[] te = tem[8].split(";");
                    String query = te[0].split(":")[1];
                    int index = Arrays.binarySearch(geneNames, query);
                    genes[index] = new Gene (query, Integer.valueOf(tem[0]), Integer.valueOf(tem[3]), Integer.valueOf(tem[4])+1, (byte)(tem[6].equals("+")? 1:0));
//                    System.out.println(genes[index]);
                }
            }
            for (int i = 0; i < info.length; i++) {
                String[] tem = info[i].split("\t");
                if (tem[2].startsWith("mRNA")) {
                    String[] te = tem[8].split(";");
                    String geneName = te[1].split(":")[1];
                    int geneIndex = Arrays.binarySearch(geneNames, geneName);
                    Transcript t = new Transcript (te[0].split(":")[1], Integer.valueOf(tem[0]), Integer.valueOf(tem[3]), Integer.valueOf(tem[4])+1, (byte)(tem[6].equals("+")? 1:0));
                    genes[geneIndex].addTranscript(t);
//                    System.out.println(t);
                }
            }
            for (int i = 0; i < genes.length; i++) genes[i].sortTranscriptsByName();
            for (int i = 0; i < info.length; i++) {
                String[] tem = info[i].split("\t");
                if (tem[2].startsWith("CDS")) {
                    String[] te = tem[8].split(";");
                    String transcriptName = te[1].split(":")[1];
                    String geneName;
                    if (transcriptName.startsWith("Zm")) {
                        geneName = transcriptName.split("_")[0];
                    }
                    else {
                        geneName = transcriptName.replaceFirst("_FGT", "_FG");
                    }
                    int geneIndex = Arrays.binarySearch(geneNames, geneName);
                    int transcriptIndex = genes[geneIndex].getTranscriptIndex(transcriptName);
                    genes[geneIndex].ts.get(transcriptIndex).addCDS(Integer.valueOf(tem[0]), Integer.valueOf(tem[3]), Integer.valueOf(tem[4])+1);
                }
                else if (tem[2].startsWith("five_prime_UTR")) {
                    String[] te = tem[8].split(";");
                    String transcriptName = te[0].split("=")[1];
                    transcriptName = transcriptName.split(":")[1];
                    String geneName;
                    if (transcriptName.startsWith("Zm")) {
                        geneName = transcriptName.split("_")[0];
                    }
                    else {
                        geneName = transcriptName.replaceFirst("_FGT", "_FG");
                    }
                    int geneIndex = Arrays.binarySearch(geneNames, geneName);
                    int transcriptIndex = genes[geneIndex].getTranscriptIndex(transcriptName);
                     genes[geneIndex].ts.get(transcriptIndex).add5UTR(Integer.valueOf(tem[0]), Integer.valueOf(tem[3]), Integer.valueOf(tem[4])+1);
                }
                else if (tem[2].startsWith("three_prime_UTR")) {
                    String[] te = tem[8].split(";");
                    String transcriptName = te[0].split("=")[1];
                    transcriptName = transcriptName.split(":")[1];
                    String geneName;
                    if (transcriptName.startsWith("Zm")) {
                        geneName = transcriptName.split("_")[0];
                    }
                    else {
                        geneName = transcriptName.replaceFirst("_FGT", "_FG");
                    }
                    int geneIndex = Arrays.binarySearch(geneNames, geneName);
                    int transcriptIndex = genes[geneIndex].getTranscriptIndex(transcriptName);
                    genes[geneIndex].ts.get(transcriptIndex).add3UTR(Integer.valueOf(tem[0]), Integer.valueOf(tem[3]), Integer.valueOf(tem[4])+1);
                }
            }
            for (int i = 0; i < this.genes.length; i++) {
                for (int j = 0; j < genes[i].ts.size(); j++) {
                    genes[i].ts.get(j).sort5UTRByPosition();
                    genes[i].ts.get(j).sortCDSByPosition();
                    genes[i].ts.get(j).sort3UTRByPosition();
                    genes[i].ts.get(j).calculateIntron();
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
    
    /**
     * Return a RangeAttribute object for from all 5UTR
     * @return 
     */
    public RangeAttribute getAllGene5UTRRange () {
        ArrayList<Range> rs = new ArrayList();
        TByteArrayList byteList = new TByteArrayList();
        TFloatArrayList floatList = new TFloatArrayList();
        for (int i = 0; i < this.getGeneNumber(); i++) {
            if (!this.isThere5UTR(i, 0)) continue;
            ArrayList<Range> l = this.get5UTRList(i, 0);
            for (int j = 0; j < l.size(); j++) rs.add(l.get(j));
            byteList.add(this.genes[i].ts.get(0).strand);
            floatList.add(Float.NaN);
        }
        return new RangeAttribute(rs.toArray(new Range[rs.size()]), "5'UTR", byteList.toArray(new byte[byteList.size()]), floatList.toArray(new float[floatList.size()]));
    }
    
    /**
     * Return a RangeAttribute object for from all CDS
     * @return 
     */
    public RangeAttribute getAllGeneCDSRange () {
        ArrayList<Range> rs = new ArrayList();
        TByteArrayList byteList = new TByteArrayList();
        TFloatArrayList floatList = new TFloatArrayList();
        for (int i = 0; i < this.getGeneNumber(); i++) {
            ArrayList<Range> l = this.getCDSList(i, 0);
            for (int j = 0; j < l.size(); j++) {
                rs.add(l.get(j));
                byteList.add(this.genes[i].ts.get(0).strand);
                floatList.add(Float.NaN);
            }
        }
        return new RangeAttribute(rs.toArray(new Range[rs.size()]), "CDS", byteList.toArray(new byte[byteList.size()]), floatList.toArray(new float[floatList.size()]));
    }
    
    /**
     * Return a RangeAttribute object for from all Intron
     * @return 
     */
    public RangeAttribute getAllIntronRange () {
        ArrayList<Range> rs = new ArrayList();
        TByteArrayList byteList = new TByteArrayList();
        TFloatArrayList floatList = new TFloatArrayList();
        for (int i = 0; i < this.getGeneNumber(); i++) {
            ArrayList<Range> l = this.getIntronList(i, 0);
            if (l.isEmpty()) continue;
            for (int j = 0; j < l.size(); j++) {
                rs.add(l.get(j));
                byteList.add(this.genes[i].ts.get(0).strand);
                floatList.add(Float.NaN);
            }
        }
        return new RangeAttribute(rs.toArray(new Range[rs.size()]), "Intron", byteList.toArray(new byte[byteList.size()]), floatList.toArray(new float[floatList.size()]));
    }
    
    /**
     * Return a RangeAttribute object for from all 3UTR
     * @return 
     */
    public RangeAttribute getAllGene3UTRRange () {
        ArrayList<Range> rs = new ArrayList();
        TByteArrayList byteList = new TByteArrayList();
        TFloatArrayList floatList = new TFloatArrayList();
        for (int i = 0; i < this.getGeneNumber(); i++) {
            if (!this.isThere3UTR(i, 0)) continue;
            ArrayList<Range> l = this.get3UTRList(i, 0);
            for (int j = 0; j < l.size(); j++) rs.add(l.get(j));
            byteList.add(this.genes[i].ts.get(0).strand);
            floatList.add(Float.NaN);
        }
        return new RangeAttribute(rs.toArray(new Range[rs.size()]), "3'UTR", byteList.toArray(new byte[byteList.size()]), floatList.toArray(new float[floatList.size()]));
    }
    
    /**
     * Return a RangeAttribute object for from all intergenic region
     * @return 
     */
    public RangeAttribute getIntergeneicRange () {
        ArrayList<Range> rs = new ArrayList();
        TByteArrayList byteList = new TByteArrayList();
        TFloatArrayList floatList = new TFloatArrayList();
        for (int i = 0; i < this.getGeneNumber()-1; i++) {
            int chr = this.genes[i].ts.get(0).transcriptRange.chr;
            int start = this.genes[i].ts.get(0).transcriptRange.end;
            if (chr != this.genes[i+1].ts.get(0).transcriptRange.chr) continue;
            int end = this.genes[i+1].ts.get(0).transcriptRange.start;
            rs.add(new Range(chr, start, end));
            byteList.add(Byte.MIN_VALUE);
            floatList.add(Float.NaN);
        }
        return new RangeAttribute(rs.toArray(new Range[rs.size()]), "Intergenic", byteList.toArray(new byte[byteList.size()]), floatList.toArray(new float[floatList.size()]));
    }
    
    /**
     * Return a RangeAttribute object for from all first transcripts
     * @return 
     */
    public RangeAttribute getAllGeneTranscriptRange () {
        ArrayList<Range> rs = new ArrayList();
        TByteArrayList byteList = new TByteArrayList();
        TFloatArrayList floatList = new TFloatArrayList();
        for (int i = 0; i < this.getGeneNumber(); i++) {
            rs.add(this.genes[i].ts.get(0).transcriptRange);
            byteList.add(this.genes[i].ts.get(0).strand);
            floatList.add(Float.NaN);
        }
        return new RangeAttribute(rs.toArray(new Range[rs.size()]), "mRNA", byteList.toArray(new byte[byteList.size()]), floatList.toArray(new float[floatList.size()]));
    }
    
    /**
     * Return a RangeAttribute object for from all nBp upstream of transcripts
     * @return 
     */
    public RangeAttribute getAllGeneUpstreamRange (int bps) {
        ArrayList<Range> rs = new ArrayList();
        TByteArrayList byteList = new TByteArrayList();
        TFloatArrayList floatList = new TFloatArrayList();
        for (int i = 0; i < this.getGeneNumber(); i++) {
            Range r = null;
            if (this.genes[i].ts.get(0).strand == 1) {
                int start = this.getTranscriptStart(i, 0)-bps;
                if (start < 0) start = 1;
                r = new Range(this.getTranscriptChromosome(i, 0), start, this.getTranscriptStart(i, 0));
            }
            else {
                int start = this.getTranscriptEnd(i, 0)+bps;
                r = new Range(this.getTranscriptChromosome(i, 0), this.getTranscriptEnd(i, 0), start);
            }
            rs.add(r);
            byteList.add(this.genes[i].ts.get(0).strand);
            floatList.add(Float.NaN);
        }
        return new RangeAttribute(rs.toArray(new Range[rs.size()]), "Upstream"+String.valueOf(bps), byteList.toArray(new byte[byteList.size()]), floatList.toArray(new float[floatList.size()]));
    }
    
    /**
     * Return a RangeAttribute object for from all nBp downstream of transcripts
     * @return 
     */
    public RangeAttribute getAllGeneDownstreamRange (int bps) {
        ArrayList<Range> rs = new ArrayList();
        TByteArrayList byteList = new TByteArrayList();
        TFloatArrayList floatList = new TFloatArrayList();
        for (int i = 0; i < this.getGeneNumber(); i++) {
            Range r = null;
            if (this.genes[i].ts.get(0).strand == 1) {
                int end = this.getTranscriptEnd(i, 0)+bps;
                r = new Range(this.getTranscriptChromosome(i, 0), this.getTranscriptEnd(i, 0), end);
            }
            else {
                int end = this.getTranscriptStart(i, 0)-bps;
                if (end < 1) end = 1;
                r = new Range(this.getTranscriptChromosome(i, 0), end, this.getTranscriptStart(i, 0));
            }
            rs.add(r);
            byteList.add(this.genes[i].ts.get(0).strand);
            floatList.add(Float.NaN);
        }
        return new RangeAttribute(rs.toArray(new Range[rs.size()]), "Downstream"+String.valueOf(bps), byteList.toArray(new byte[byteList.size()]), floatList.toArray(new float[floatList.size()]));
    }
    
    class Gene implements Comparable<Gene> {
        String geneName = null;
        Range geneRange = null;
        byte strand = Byte.MIN_VALUE;
        ArrayList<Transcript> ts = new ArrayList();
        
        public Gene (String geneName, int chr, int start, int end, byte strand) {
            this.geneName = geneName;
            this.strand = strand;
            geneRange = new Range(chr, start, end);
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
        
        public int getTranscriptIndex (String transcriptName) {
            return Collections.binarySearch(ts, new Transcript(transcriptName));
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
        ArrayList<Range> cdsList = new ArrayList();
        ArrayList<Range> intronList = new ArrayList();
        ArrayList<Range> utr5List = new ArrayList();
        ArrayList<Range> utr3List = new ArrayList();
        
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

