package daxing.v2.ancestryHmm;

import com.ibm.icu.text.NumberFormat;
import daxing.common.factors.WheatLineage;
import daxing.common.table.RowTableTool;
import daxing.common.utiles.IOTool;
import gnu.trove.list.array.TIntArrayList;
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

        Map<String, String> sample2PopMap = RowTableTool.getMap(sample2PopFile, 0 ,1, false);
        String[] outNames = files.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz","_alleleCount.txt")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->{
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedReader brSNPIn = IOTool.getReader(prunedInSNPFile);
                 BufferedWriter bw = IOTool.getWriter(new File(alleleCountDir, outNames[e]))) {
                RecombinationMaps recombinationMaps = new RecombinationMaps(recombinationMapFile);
                List<String> snpInList = new ArrayList<>(42*1000*1000);
                String line;
                List<String> temp;
                while ((line=brSNPIn.readLine())!=null){
                    temp = PStringUtils.fastSplit(line);
                    snpInList.add(temp.get(0));
                }
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

                while ((line = br.readLine())!=null){
                    totalSNPCountInVCF++;
                    temp = PStringUtils.fastSplit(line);

                    // pruned SNP
                    snpIndex = Collections.binarySearch(snpInList, temp.get(2));
                    if (snpIndex < 0) continue;

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
                        sb.append(refAlleleCount).append(",").append(altAlleleCount).append("\t");
                    }
                    // make sure we have enough samples at this position
                    if (ifContinue) continue;

                    // make sure the largest population-pair frequency difference greater than 0.1
                    ifContinue=true;
                    for (int i = 0; i < refFrequency.length-1; i++) {
                        for (int j = i+1; j < refFrequency.length; j++) {
                            if (refFrequency[j]-refFrequency[i] > minRefPopPairAlleleFreDifference){
                                ifContinue = false;
                            }
                        }
                    }
                    if (ifContinue) continue;


                    retainSNPCount++;

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
                        sb.append(refAlleleCount).append(",").append(altAlleleCount).append("\t");
                    }
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                System.out.println("Total "+numberFormat.format(totalSNPCountInVCF)+" SNPs in "+files.get(e).getName());
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
