package daxing.v2.ancestryHmm;

import com.google.common.collect.Comparators;
import com.ibm.icu.text.NumberFormat;
import daxing.common.factors.WheatLineage;
import daxing.common.table.RowTableTool;
import daxing.common.utiles.IOTool;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Panel {

    /**
     * allele count file for ancestry_HMM software
     */
    public static void getAlleleCount(List<File> files, String sample2PopFile, String recombinationMapFile,
                                      String prunedInSNPFile,
                                      String alleleCountDir){

        int minRefPopCount = 10;
        double minRefPopPairAlleleFreDifference = 0.1;

        System.out.println("min reference population number is be setting "+minRefPopCount);
        System.out.println("min reference population pair allele frequency difference is "+minRefPopPairAlleleFreDifference);

        Map<String, String> sample2PopMap = RowTableTool.getMap(sample2PopFile, 0 ,1, false);
        String[] outNames =files.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz","_pruned.panel")).toArray(String[]::new);
        IntStream.range(0, files.size()).forEach(e->{
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedReader brSNPIn = IOTool.getReader(prunedInSNPFile);
                 BufferedWriter bw = IOTool.getWriter(new File(alleleCountDir, outNames[e]))) {
                RecombinationMaps recombinationMaps = new RecombinationMaps(recombinationMapFile);
                List<String> snpInList = new ArrayList<>(42*1000*1000);
                String line;
                List<String> temp;
                long start = System.nanoTime();
                while ((line=brSNPIn.readLine())!=null){
                    snpInList.add(line);
                }
                System.out.println(" reading LD pruned SNP take "+ Benchmark.getTimeSpanSeconds(start)+ " s");
                Comparator<String> c1 = Comparator.comparing(l->Integer.parseInt(PStringUtils.fastSplit(l,"-").get(0)));
                Comparator<String> c2 = c1.thenComparing(l->Integer.parseInt(PStringUtils.fastSplit(l, "-").get(1)));
                Collections.sort(snpInList, c2);

                NumberFormat numberFormat = NumberFormat.getNumberInstance();
                numberFormat.setGroupingUsed(true);
                System.out.println("Total "+ numberFormat.format(snpInList.size())+" SNPs in "+new File(prunedInSNPFile).getName());
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
                boolean ifContinue;
                double lastGeneticPosition=0, currentGeneticPosition, distance;
                String currentChr="NA";
                int snpIndex;
                int retainSNPCount=0, totalSNPCountInVCF=0;
                int removedCount_sample=0, removedCount_freDiff=0, removedCount_prunedSNP=0;

                while ((line = br.readLine())!=null){
                    totalSNPCountInVCF++;
//                    if (totalSNPCountInVCF % 5000000 ==0){
//                        System.out.println("reading "+totalSNPCountInVCF+" variants from "+files.get(e).getName());
//                    }
                    temp = PStringUtils.fastSplit(line);

                    // pruned SNP
                    snpIndex = Collections.binarySearch(snpInList, temp.get(2), c2);
                    if (snpIndex < 0){
                        removedCount_prunedSNP++;
                        continue;
                    }

                    sb.setLength(0);
                    sb.append(String.join("\t", temp.subList(0, 2))).append("\t");
                    refFrequency = new double[refPops.length];
                    Arrays.fill(refFrequency, -1);
                    ifContinue= false;
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
                        if (refSampleCount < minRefPopCount) {
                            ifContinue = true;
                        }
                        refFrequency[i] = (double) altAlleleCount/(altAlleleCount+refAlleleCount);
                        sb.append(refAlleleCount).append("\t").append(altAlleleCount).append("\t");
                    }
                    // make sure we have enough samples at this position
                    if (ifContinue){
                        removedCount_sample++;
                        continue;
                    }

                    // make sure the largest population-pair frequency difference greater than 0.1
                    ifContinue=true;
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
                        System.out.println("writing "+retainSNPCount+" variants to "+new File(alleleCountDir,outNames[e]).getName());
                    }

                    // add genetic map distance
                    if (currentChr.equals("NA")){
                        currentChr = temp.get(0);
                    }
                    currentGeneticPosition=recombinationMaps.getGeneticPositions(temp.get(0),
                            Integer.parseInt(temp.get(1)));
                    if (currentChr.equals(temp.get(0))){
                        distance = currentGeneticPosition - lastGeneticPosition;
                        lastGeneticPosition = currentGeneticPosition;
                    }else {
                        distance = currentGeneticPosition - 0;
                        lastGeneticPosition = currentGeneticPosition;
                        currentChr = temp.get(0);
                    }
                    sb.append(distance).append("\t");


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

                System.out.println("Total "+numberFormat.format(totalSNPCountInVCF)+" SNPs in "+files.get(e).getName());
                System.out.println("Sample count filter (less than "+minRefPopCount+" ):"+ removedCount_sample);
                System.out.println("frequency difference filter: "+ removedCount_freDiff);
                System.out.println("pruned snp filter: "+ removedCount_prunedSNP);
                System.out.println(outNames[e]+" retain "+numberFormat.format(retainSNPCount)+" SNPs after filtering");

            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });

    }

    /**
     *
     * @param genotypeDir contain two file, chrAB.vcf and chrD.vcf
     * @param sample2PopFile_AB
     * @param sample2PopFile_D
     * @param alleleCountDir
     */
    public static void getAlleleCount(String genotypeDir, String sample2PopFile_AB, String sample2PopFile_D,
                                      String recombinationMapFile, String prunedInSNPFile_AB,
                                      String prunedInSNPFile_D, String alleleCountDir){

        List<File> files = IOTool.getFileListInDirEndsWith(genotypeDir, ".vcf.gz");
        Predicate<File> d = file -> WheatLineage.valueOf(file.getName().substring(3,4)).name().equals("D");
        List<File> abFiles = files.stream().filter(d.negate()).collect(Collectors.toList());
        List<File> dFiles = files.stream().filter(d).collect(Collectors.toList());

        getAlleleCount(abFiles, sample2PopFile_AB, recombinationMapFile, prunedInSNPFile_AB, alleleCountDir);
        getAlleleCount(dFiles, sample2PopFile_D, recombinationMapFile, prunedInSNPFile_D, alleleCountDir);
    }
}
