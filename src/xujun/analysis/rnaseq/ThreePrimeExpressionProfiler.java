/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xujun.analysis.rnaseq;

import format.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import utils.IOUtils;

/**
 *
 * @author feilu
 */
public class ThreePrimeExpressionProfiler {
    //The path of reference genome with .fai file in the same folder
    String referenceGenomeS = null;
    //The SampleInformation file (with header), the format is Taxa\tBarcode\tPlateName\tFastqPath
    String sampleInformationFileS = null;
    //The gene annotation file (GTF format)
    String geneAnnotationFileS = null;
    //The path of STAR alignment software
    String starPath = null;
    //The directory of output
    String outputDirS = null;
    
    String[] subDirS = {"subFastqs", "sams", "geneCount"};
    List<String> fqFileSList = null;
    
    int[] barcodeLengths = null;
    
    List<String>[] barcodeLists = null;
    
    HashMap<String, String>[] barcodeTaxaMaps = null;
    
    List<String>[] taxaLists = null;
    
    public ThreePrimeExpressionProfiler (String parameterFileS) {
        this.parseParameters(parameterFileS);
        this.parseFq();
    }
    
    private void parseFq () {
        fqFileSList.parallelStream().forEach(f -> {
            int fqIndex = Collections.binarySearch(this.fqFileSList, f);
            String subFqDirS = new File (this.outputDirS, subDirS[0]).getAbsolutePath();
            List<String> barcodeList = barcodeLists[fqIndex];
            String[] subFqFileS = new String[barcodeList.size()];
            HashMap<String, String> btMap = barcodeTaxaMaps[fqIndex];
            Set<String> barcodeSet = btMap.keySet();
            BufferedWriter[] bws = new BufferedWriter[subFqFileS.length];
            HashMap<String, BufferedWriter> barcodeWriterMap = new HashMap<>();
            for (int i = 0; i < subFqFileS.length; i++) {
                String taxon = btMap.get(barcodeList.get(i));
                subFqFileS[i] = new File(subFqDirS, taxon+".fq.gz").getAbsolutePath();
                bws[i] = IOUtils.getTextGzipWriter(subFqFileS[i]);
                barcodeWriterMap.put(barcodeList.get(i), bws[i]);
            }
            int barcodeLength = this.barcodeLengths[fqIndex];
            try {
                BufferedReader br = null;
                if (f.endsWith(".gz")) {
                    br = IOUtils.getTextGzipReader(f);
                }
                else {
                    br = IOUtils.getTextReader(f);
                }
                String temp = null;
                String seq = null;
                String currentBarcode = null;
                BufferedWriter tw = null;
                int cnt = 0;
                while((temp = br.readLine())!=null){
                    seq = br.readLine();
                    currentBarcode = seq.substring(0, barcodeLength);
                    int cutIndex = 0;
                    if (barcodeSet.contains(currentBarcode)) {
                        tw = barcodeWriterMap.get(currentBarcode);
                        for(int i=barcodeLength; i<seq.length(); i+=4){
                            if (i+4 > seq.length()) {
                                br.readLine();br.readLine();
                                continue;
                            }
                            byte[] seqB = seq.substring(i, i+4).getBytes();
                            int value = this.getsum(seqB);
                            //317 = TTAA
                            if (value < 317) {
                                cutIndex = i;
                                tw.write(temp);
                                tw.newLine();
                                tw.write(seq.substring(cutIndex));
                                tw.newLine();
                                tw.write(br.readLine());
                                tw.newLine();
                                tw.write(br.readLine().substring(cutIndex));
                                cnt++;
                                if (cnt%100000 == 0) System.out.println("Wrote " + String.valueOf(cnt) + " sequences. " + f);
                                break;
                            }              
                        }
                    }
                    else {
                        br.readLine();br.readLine();
                        continue;
                    }
                }
                for (int i = 0; i < subFqFileS.length; i++) {
                    bws[i].flush();
                    bws[i].close();
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }
    
    private int getsum(byte[] seqB){
        int sum=0;
        for(int i=0;i<seqB.length;i++){           
            sum+=seqB[i];
        }
        return sum;
    }
    
    private void parseParameters(String infileS) {
        List<String> pLineList = new ArrayList<>();
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = null;
            boolean ifOut = false;
            if (!(temp = br.readLine()).equals("Three' Expression Profiler (TEP)")) ifOut = true;
            if (!(temp = br.readLine()).equals("Author: Jun Xu, Fei Lu")) ifOut = true;
            if (!(temp = br.readLine()).equals("Email: liuzhongxujun@163.com; flu@genetics.ac.cn")) ifOut = true;
            if (!(temp = br.readLine()).equals("Homepage: http://plantgeneticslab.weebly.com/")) ifOut = true;
            if (ifOut) {
                System.out.println("Thanks for using TEP.");
                System.out.println("Please keep the authorship in the parameter file. Program stops.");
                System.exit(0);
            }
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) continue;
                if (temp.isEmpty()) continue;
                pLineList.add(temp);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        this.referenceGenomeS = pLineList.get(0);
        this.sampleInformationFileS = pLineList.get(1);
        this.geneAnnotationFileS = pLineList.get(2);
        this.starPath = pLineList.get(3);
        this.outputDirS = pLineList.get(4);
        this.processTaxaAndBarcode();
        for (int i = 0; i < this.subDirS.length; i++) {
            new File(this.outputDirS, subDirS[i]).mkdir();
        }
    }
    
    private void processTaxaAndBarcode () {
        RowTable<String> t = new RowTable<>(this.sampleInformationFileS);
        Set<String> fqSet = new HashSet<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            fqSet.add(t.getCell(i, 3));
        }
        fqFileSList = new ArrayList<>(fqSet);
        Collections.sort(fqFileSList);
        barcodeLengths = new int[fqFileSList.size()];
        barcodeLists = new ArrayList[fqFileSList.size()];
        taxaLists = new ArrayList[fqFileSList.size()];
        barcodeTaxaMaps = new HashMap[fqFileSList.size()];
        int[] cnts = new int[fqFileSList.size()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = Collections.binarySearch(fqFileSList, t.getCell(i, 3));
            cnts[index]++;
        }
        for (int i = 0; i < cnts.length; i++) {
            barcodeLists[i] = new ArrayList<>();
            taxaLists[i] = new ArrayList<>();
            barcodeTaxaMaps[i] = new HashMap<>();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = Collections.binarySearch(fqFileSList, t.getCell(i, 3));
            cnts[index]++;
            String taxon = t.getCell(i, 0) + "_"+ t.getCell(i, 2);
            taxaLists[index].add(taxon);
            barcodeLists[index].add(t.getCell(i, 1));
            barcodeTaxaMaps[index].put(t.getCell(i, 1), taxon);
            barcodeLengths[index] = t.getCell(i, 1).length();
        }
    }
    
    public static void main(String arg[]) {
        new ThreePrimeExpressionProfiler(arg[0]);
    }
    
}
