package daxing.v2.mosaic;

import daxing.common.table.RowTableTool;
import daxing.common.utiles.IOTool;
import daxing.v2.localAncestryInfer.Simulation;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import pgl.infra.utils.PStringUtils;
import java.io.*;
import java.util.*;
import java.util.stream.IntStream;

public class Mosaic_utils {



    public static void transformVCF_to_Mosaic_inputFormat(File genotypeFile, File taxaInfoFile,
                                                          File recombinationMap,
                                                          File outDir_Mosaic){

        String[] subDir = {"mosaicData"};
        File[] subDirFile = new File[subDir.length];
        for (int i = 0; i < subDir.length; i++) {
            subDirFile[i] = new File(outDir_Mosaic, subDir[i]);
            subDirFile[i].mkdir();
        }
        Map<String, String> taxonPopMap = RowTableTool.getMap(taxaInfoFile.getAbsolutePath(), 0, 2);
        IntList posList = new IntArrayList();
        DoubleList geneticsMapList = new DoubleArrayList();
        try (BufferedReader brGenotype = IOTool.getReader(genotypeFile);
             BufferedReader brRecombinationMap = IOTool.getReader(recombinationMap)) {
            String line, popName;
            List<String> temp, headerList;
            StringBuilder sb = new StringBuilder();

            brRecombinationMap.readLine();
            while ((line=brRecombinationMap.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                posList.add(Integer.parseInt(temp.get(0)));
                geneticsMapList.add(Double.parseDouble(temp.get(2)));
            }

            // rates.1
            File rate1 = new File(subDirFile[0], "rates.1");
            BufferedWriter bwRates1 = IOTool.getWriter(rate1);
            sb.setLength(0);
            sb.append(":sites:").append(posList.size());
            bwRates1.write(sb.toString());
            bwRates1.newLine();
            sb.setLength(0);
            for (int i = 0; i < posList.size(); i++) {
                sb.append(posList.getInt(i)).append(" ");
            }
            sb.deleteCharAt(sb.length()-1);
            bwRates1.write(sb.toString());
            bwRates1.newLine();
            sb.setLength(0);
            for (int i = 0; i < geneticsMapList.size(); i++) {
                sb.append(geneticsMapList.getDouble(i)).append(" ");
            }
            sb.deleteCharAt(sb.length()-1);
            bwRates1.write(sb.toString());
            bwRates1.newLine();
            bwRates1.flush();
            bwRates1.close();

            // pop2TaxaMap
            while ((line=brGenotype.readLine()).startsWith("##")){}
            headerList = PStringUtils.fastSplit(line);
            Map<String,IntList> popIndexMap = new HashMap<>();
            Set<String> popSet = RowTableTool.getColumnSet(taxaInfoFile.getAbsolutePath(), 2);
            List<String> popList = new ArrayList<>(popSet);
            Collections.sort(popList);
            for (int i = 0; i < popList.size(); i++) {
                popIndexMap.put(popList.get(i), new IntArrayList());
            }
            String taxonName;
            for (int i = 9; i < headerList.size(); i++) {
                taxonName = headerList.get(i);
                popName = taxonPopMap.get(taxonName);
                popIndexMap.get(popName).add(i);
            }

            // sample.names
            File sampleNamesFile = new File(subDirFile[0], "sample.names");
            BufferedWriter bwSampleName = IOTool.getWriter(sampleNamesFile);
            for (int i = 9; i < headerList.size(); i++) {
                sb.setLength(0);
                popName = taxonPopMap.get(headerList.get(i));
                sb.append(popName).append(" ").append(headerList.get(i)).append(" ");
                sb.append("0 0 0 2 -9");
                bwSampleName.write(sb.toString());
                bwSampleName.newLine();
            }
            bwSampleName.flush();
            bwSampleName.close();

            // snpfile and genofile
            File snpFile1 = new File(subDirFile[0], "snpfile.1");
            BufferedWriter bwSNPFile = IOTool.getWriter(snpFile1);
            Map<String,BufferedWriter> pop2BWMap = new HashMap<>();
            File genoFile;
            BufferedWriter bw;
            for (int i = 0; i < popList.size(); i++) {
                genoFile = new File(subDirFile[0], popList.get(i)+"genofile.1");
                bw = IOTool.getWriter(genoFile);
                pop2BWMap.put(popList.get(i), bw);
            }
            IntList indexesList;
            while ((line=brGenotype.readLine())!=null){
                temp =PStringUtils.fastSplit(line);
                sb.setLength(0);
                sb.append("         ").append(temp.get(2)).append(" ");
                sb.append(temp.get(0)).append(" ").append(0).append(" ");
                sb.append(temp.get(1)).append(" ");
                sb.append(temp.get(3)).append(" ");
                sb.append(temp.get(4)).append(" ");
                bwSNPFile.write(sb.toString());
                bwSNPFile.newLine();

                for (Map.Entry<String,IntList> entry: popIndexMap.entrySet()){
                    popName = entry.getKey();
                    indexesList = entry.getValue();
                    bw = pop2BWMap.get(popName);
                    sb.setLength(0);
                    for (int index : indexesList){
                        sb.append(temp.get(index));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bwSNPFile.flush();
            bwSNPFile.close();

            for (Map.Entry<String,BufferedWriter> entry : pop2BWMap.entrySet()){
                entry.getValue().flush();
                entry.getValue().close();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void transformVCF_to_Mosaic_inputFormat(String genotypeFileDir, String taxaInfoFile,
                                                          String recombinationMapDir,
                                                          String outDir_Mosaic){
        List<File> genotypeFiles = IOTool.getFileListInDirEndsWith(genotypeFileDir, "Ne.vcf");
        List<File> recombinationFiles = IOTool.getFileListInDirEndsWith(recombinationMapDir, "Map.txt");
        double[] admixture_proportion = new double[genotypeFiles.size()];
        int[] admixture_generation = new int[genotypeFiles.size()];
        double[] divergenceTime = new double[genotypeFiles.size()];
        List<String> temp;
        for (int i = 0; i < genotypeFiles.size(); i++) {
            temp = PStringUtils.fastSplit(genotypeFiles.get(i).getName(), "_");
            admixture_proportion[i] = Double.parseDouble(temp.get(4).replaceAll("proportion",""));
            admixture_generation[i] = Integer.parseInt(temp.get(5).replaceAll("genetration",""));
            divergenceTime[i] = Double.parseDouble(temp.get(7).replaceAll("divergence", "").replaceAll("Ne.vcf",""));
        }
        File[] subDir = new File[admixture_proportion.length];
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < admixture_proportion.length; i++) {
            sb.setLength(0);
            sb.append("simulate_two_way_admixture_proportion");
            sb.append(admixture_proportion[i]).append("_genetration");
            sb.append(admixture_generation[i]).append("_P1P2_divergence");
            sb.append(divergenceTime[i]).append("Ne");
            subDir[i] = new File(outDir_Mosaic, sb.toString());
            subDir[i].mkdir();
        }

        for (int i = 0; i < genotypeFiles.size(); i++) {
            Mosaic_utils.transformVCF_to_Mosaic_inputFormat(genotypeFiles.get(i), new File(taxaInfoFile),
                    recombinationFiles.get(i), subDir[i]);
        }

    }

    public static double[][][] build_mosaic_sh(String fastFilePath, String admixedPop,
                                  String mosaicRunDir, String logFile, int admixedIndividualNumber, int variantsNum,
                                               int nWayAdmixture, int coresNum, int chrID){
        List<File> mosaicDirs = IOTool.getDirListInDir(mosaicRunDir);
        String mosaicRpath = Simulation.class.getResource("mosaic.R").getPath();
        String extractMosaicLocalAncestryRpath = Simulation.class.getResource("extractMosaicLocalAncestry.R").getPath();
        double[][][] mosaicDir_taxa_variants_localAncestry = new double[mosaicDirs.size()][][];
        for (int i = 0; i < mosaicDir_taxa_variants_localAncestry.length; i++) {
            mosaicDir_taxa_variants_localAncestry[i] = new double[admixedIndividualNumber][];
            for (int j = 0; j < mosaicDir_taxa_variants_localAncestry[i].length; j++) {
                mosaicDir_taxa_variants_localAncestry[i][j] = new double[variantsNum];
                Arrays.fill(mosaicDir_taxa_variants_localAncestry[i][j], -1);
            }
        }
        IntStream.range(0, mosaicDirs.size()).forEach(i->{
            File currentWorkingDir = mosaicDirs.get(i);
            File mosaicInputDataDir = currentWorkingDir.listFiles(dir -> dir.getName().equals("mosaicData"))[0];
            try {
                // build and run mosaic command
                StringBuilder sb = new StringBuilder();
                sb.setLength(0);
                sb.append("/usr/local/bin/Rscript ").append(mosaicRpath).append(" ");
                sb.append(admixedPop).append(" ").append(mosaicInputDataDir.getPath()).append("/ ");
                sb.append("-f ").append(fastFilePath).append(" ");
                sb.append("-a ").append(nWayAdmixture).append(" ");
                sb.append("-m ").append(coresNum).append(" ");
                sb.append("-c ").append(chrID).append(":").append(chrID);
                ProcessBuilder processBuilder = new ProcessBuilder(sb.toString());
                processBuilder.directory(currentWorkingDir);
                processBuilder.redirectError(ProcessBuilder.Redirect.appendTo(new File(logFile)));
                Process process;
                String line;
                List<String> temp;
                process = processBuilder.start();
                int exitCode = process.waitFor();
                assert exitCode == 0 : currentWorkingDir.getName()+" mosaic command run failed";

                // extract mosaic local ancestry result
//                File mosaicResDir = currentWorkingDir.listFiles(dir -> dir.getName().startsWith("MOSAIC_RESULTS"))[0];
//                File[] rDataFiles = mosaicResDir.listFiles(dir -> dir.getName().endsWith("RData"));
//                String modelParametersRData = Arrays.stream(rDataFiles).filter(file -> !file.getName().startsWith("localanc")).toArray(File[]::new)[0].getAbsolutePath();
//                String localAncestryRData = Arrays.stream(rDataFiles).filter(file -> file.getName().startsWith("localanc")).toArray(File[]::new)[0].getAbsolutePath();
//                sb.setLength(0);
//                sb.append("Rscript ").append(extractMosaicLocalAncestryRpath).append(" ");
//                sb.append(modelParametersRData).append(" ").append(localAncestryRData).append(" ");
//                sb.append(mosaicInputDataDir);
//                processBuilder = new ProcessBuilder(sb.toString());
//                processBuilder.directory(currentWorkingDir);
//                processBuilder.redirectError(ProcessBuilder.Redirect.appendTo(new File(logFile)));
//                process = processBuilder.start();
//
//                BufferedReader br = new BufferedReader(new InputStreamReader(process.getInputStream()));
//                int variantsIndex = 0;
//                while ((line=br.readLine())!=null){
//                    temp =PStringUtils.fastSplit(line);
//                    for (int j = 0; j < temp.size(); j++) {
//                        mosaicDir_taxa_variants_localAncestry[i][j][variantsIndex] = Double.parseDouble(temp.get(j));
//                    }
//                    variantsIndex++;
//                }
//                br.close();
//                exitCode = process.waitFor();
//                assert exitCode == 0 : currentWorkingDir.getName()+" extract mosaic local ancestry failed";
            } catch (IOException | InterruptedException e) {
                e.printStackTrace();
            }

        });
        return mosaicDir_taxa_variants_localAncestry;
    }


}
