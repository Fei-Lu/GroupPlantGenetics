package daxing.v2.ancestryHmm;

import daxing.common.factors.HexaploidBySubcontinent;
import daxing.common.table.RowTableTool;
import daxing.common.utiles.IOTool;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.*;

public class Panel {

    /**
     * allele count file for ancestry_HMM software
     */
    public static void preparePanelForAncestryHMM(File genotypeFile, File sample2PopFile, File recombinationMapFile,
                                                  File prunedInSNPFile,
                                                  File outFilePanel){
        long start = System.nanoTime();
        int minRefPopCount = 10;
        double minRefPopPairAlleleFreDifference = 0.1;
        System.out.println("min reference population number is be setting "+minRefPopCount);
        System.out.println("min reference population pair allele frequency difference is "+minRefPopPairAlleleFreDifference);
        Map<String, String> sample2PopMap = RowTableTool.getMap(sample2PopFile.getAbsolutePath(), 0 ,1, false);
        try (BufferedReader br = IOTool.getReader(genotypeFile);
             BufferedReader brSNPIn = IOTool.getReader(prunedInSNPFile);
             BufferedWriter bw = IOTool.getWriter(outFilePanel)) {
            RecombinationMaps recombinationMaps = new RecombinationMaps(recombinationMapFile.getAbsolutePath());
            List<String> snpInList = new ArrayList<>(42*1000*1000);
            String line;
            List<String> temp;
            while ((line=brSNPIn.readLine())!=null){
                snpInList.add(line);
            }
            System.out.println(" reading LD pruned SNP take "+ Benchmark.getTimeSpanSeconds(start)+ " s");
            Comparator<String> c1 = Comparator.comparing(l->Integer.parseInt(PStringUtils.fastSplit(l,"-").get(0)));
            Comparator<String> c2 = c1.thenComparing(l->Integer.parseInt(PStringUtils.fastSplit(l, "-").get(1)));
            Collections.sort(snpInList, c2);

            NumberFormat numberFormat = NumberFormat.getNumberInstance();
            numberFormat.setGroupingUsed(true);
            System.out.println("Total "+ numberFormat.format(snpInList.size())+" SNPs in "+prunedInSNPFile.getName());
            while ((line=br.readLine()).startsWith("##")) continue;
            temp = PStringUtils.fastSplit(line);
            List<String> headerList = temp.subList(9, temp.size());
            String[] pops = sample2PopMap.values().stream().distinct().sorted().toArray(String[]::new);
            String[] refPops = Arrays.stream(pops).filter(s -> !s.equals("admixed")).sorted().toArray(String[]::new);
            Map<String, TIntArrayList> pop2IndexesMap = new HashMap<>();
            for (int i = 0; i < pops.length; i++) {
                pop2IndexesMap.put(pops[i], new TIntArrayList());
            }
            String pop;
            for (int i = 0; i < headerList.size(); i++) {
                if (!sample2PopMap.containsKey(headerList.get(i))) continue;
                pop = sample2PopMap.get(headerList.get(i));
                pop2IndexesMap.get(pop).add(i);
            }
            TIntArrayList indexes;
            StringBuilder sb = new StringBuilder();
            int refAlleleCount, altAlleleCount, refSampleCount;
            double[] refFrequency;
            int lastPosition=-1, currentPosition=-1, distance=-1;
            double lastGeneticPosition=0, currentGeneticPosition, geneticsDistance;
            String currentChr="NA";
            int snpIndex;
            int retainSNPCount=0, totalSNPCountInVCF=0;
            int removedCount_sample=0, removedCount_freDiff=0, removedCount_prunedSNP=0, removed_distance=0;

            while ((line = br.readLine())!=null){
                totalSNPCountInVCF++;
                if (totalSNPCountInVCF % 5000000 ==0){
                    System.out.println("reading "+totalSNPCountInVCF+" variants from "+genotypeFile.getName());
                }
                temp = PStringUtils.fastSplit(line);

                // pruned SNP
                snpIndex = Collections.binarySearch(snpInList, temp.get(2), c2);
                if (snpIndex < 0){
                    removedCount_prunedSNP++;
                    continue;
                }

                // filter SNPs less than 1000bp
                if (currentChr.equals("NA")){
                    currentChr=temp.get(0);
                    lastPosition = 2*1000*1000*1000;
                }
                if (currentChr.equals(temp.get(0))){
                    currentPosition = Integer.parseInt(temp.get(1));
                    distance = lastPosition-currentPosition;
                    lastPosition=currentPosition;
                }else {
                    currentPosition = Integer.parseInt(temp.get(1));
                    lastPosition = 2*1000*1000*1000;
                    distance = lastPosition-currentPosition;
                    lastPosition=currentPosition;
                }
                if (distance < 1000) {
                    removed_distance++;
                    continue;
                }

                sb.setLength(0);
                sb.append(String.join("\t", temp.subList(0, 2))).append("\t");
                refFrequency = new double[refPops.length];
                Arrays.fill(refFrequency, -1);
                int refSampleMinCount=10000;
                for (int i = 0; i < refPops.length; i++) {
                    indexes = pop2IndexesMap.get(refPops[i]);
                    refAlleleCount = 0;
                    altAlleleCount = 0;
                    refSampleCount = 0;
                    for (int j = 0; j < indexes.size(); j++) {
                        if (temp.get(indexes.get(j)).startsWith("0/0")){
                            refAlleleCount += 2;
                            refSampleCount +=1;
                        }else if (temp.get(indexes.get(j)).startsWith("1/1")){
                            altAlleleCount += 2;
                            refSampleCount +=1;
                        }else if (temp.get(indexes.get(j)).startsWith("0/1")){
                            refAlleleCount += 1;
                            altAlleleCount += 1;
                            refSampleCount +=1;
                        }else if (temp.get(indexes.get(j)).startsWith("0")){
                            refAlleleCount += 1;
                            refSampleCount +=1;
                        }else if (temp.get(indexes.get(j)).startsWith("1")){
                            altAlleleCount += 1;
                            refSampleCount +=1;
                        }
                    }
                    if (refSampleCount < refSampleMinCount) {
                        refSampleMinCount = refSampleCount;
                    }
                    refFrequency[i] = (double) altAlleleCount/(altAlleleCount+refAlleleCount);
                    sb.append(refAlleleCount).append("\t").append(altAlleleCount).append("\t");
                }
                // make sure we have enough samples at this position
                if (refSampleMinCount < minRefPopCount){
                    removedCount_sample++;
                    continue;
                }

                // make sure the largest population-pair frequency difference greater than 0.1
                double maxFreDiff = 0;
                for (int i = 0; i < refFrequency.length-1; i++) {
                    for (int j = i+1; j < refFrequency.length; j++) {
                        if (refFrequency[j]-refFrequency[i] > maxFreDiff){
                            maxFreDiff = refFrequency[j]-refFrequency[i];
                        }
                    }
                }
                if (maxFreDiff < minRefPopPairAlleleFreDifference){
                    removedCount_freDiff++;
                    continue;
                }

                retainSNPCount++;
                if (retainSNPCount % 5000000 ==0){
                    System.out.println("writing "+retainSNPCount+" variants to "+outFilePanel.getName());
                }

                // add genetic map distance
                currentGeneticPosition=recombinationMaps.getGeneticPositions(temp.get(0),
                        Integer.parseInt(temp.get(1)));
                if (currentChr.equals(temp.get(0))){
                    geneticsDistance = currentGeneticPosition - lastGeneticPosition;
                    lastGeneticPosition = currentGeneticPosition;
                }else {
                    geneticsDistance = currentGeneticPosition - 0;
                    lastGeneticPosition = currentGeneticPosition;
                    currentChr=temp.get(0);
                }
                sb.append(geneticsDistance).append("\t");


                indexes = pop2IndexesMap.get("admixed");
                for (int i = 0; i < indexes.size(); i++) {
                    refAlleleCount = 0;
                    altAlleleCount = 0;
                    if (temp.get(indexes.get(i)).startsWith("0/0")){
                        refAlleleCount += 2;
                    }else if (temp.get(indexes.get(i)).startsWith("1/1")){
                        altAlleleCount += 2;
                    }else if (temp.get(indexes.get(i)).startsWith("0/1")){
                        refAlleleCount += 1;
                        altAlleleCount += 1;
                    }else if (temp.get(indexes.get(i)).startsWith("0")){
                        refAlleleCount +=1;
                    }else if (temp.get(indexes.get(i)).startsWith("1")){
                        altAlleleCount +=1;
                    }
                    sb.append(refAlleleCount).append("\t").append(altAlleleCount).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();

            System.out.println("Total "+numberFormat.format(totalSNPCountInVCF)+" SNPs in "+genotypeFile.getName());
            System.out.println("SNP distance less than 1000 bp "+ removed_distance);
            System.out.println("Sample count filter (less than "+minRefPopCount+" ):"+ removedCount_sample);
            System.out.println("frequency difference filter: "+ removedCount_freDiff);
            System.out.println("pruned snp filter: "+ removedCount_prunedSNP);
            System.out.println(outFilePanel.getName()+" retain "+numberFormat.format(retainSNPCount)+" SNPs after filtering");

        } catch (IOException ioException) {
            ioException.printStackTrace();
        }


        System.out.println("Prepare panel for ancestryHMM take "+ Benchmark.getTimeSpanMinutes(start)+ " minutes");
    }

    /**
     *
     * @param genotypeDir contain two file, chrAB.vcf and chrD.vcf
     * @param sample2PopInfoDir
     * @param recombinationMapFile
     * @param prunedInSNPFile_AB
     * @param prunedInSNPFile_D
     * @param outDirPanel
     */
    public static void preparePanelForAncestryHMM(String genotypeDir, String sample2PopInfoDir,
                                                  String recombinationMapFile, String prunedInSNPFile_AB,
                                                  String prunedInSNPFile_D, String outDirPanel){

        List<File> genotypeFiles = IOTool.getFileListInDirEndsWith(genotypeDir, ".vcf.gz");
        String[] outNamesPrefix = genotypeFiles.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz", "")).toArray(String[]::new);
        for (HexaploidBySubcontinent hexaploidBySubcontinent: HexaploidBySubcontinent.values()){
            preparePanelForAncestryHMM(genotypeFiles.get(0), new File(sample2PopInfoDir,
                    "sampleInfoAB_"+hexaploidBySubcontinent+".txt"), new File(recombinationMapFile),
                    new File(prunedInSNPFile_AB), new File(outDirPanel,
                            outNamesPrefix[0]+"_"+hexaploidBySubcontinent.name()+
                            ".panel"));
        }

        for (HexaploidBySubcontinent hexaploidBySubcontinent: HexaploidBySubcontinent.values()){
            preparePanelForAncestryHMM(genotypeFiles.get(1), new File(sample2PopInfoDir,
                            "sampleInfoD_"+hexaploidBySubcontinent+".txt"), new File(recombinationMapFile),
                    new File(prunedInSNPFile_D), new File(outDirPanel, outNamesPrefix[1]+"_"+hexaploidBySubcontinent+
                            ".panel"));
        }
    }
}
