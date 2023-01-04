package daxing.v2.mosaic;

import com.google.common.collect.BiMap;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.HashBiMap;
import com.google.common.collect.Table;
import com.google.common.primitives.Ints;
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

public class Mosaic_runner {

    String subDirName = "mosaicData";

    /**
     * ChromosomeID\tGenotype\tRecombinationMap
     * 1\t/Users/xudaxing/chr.vcf\t/Users/xudaxing/recombinationMap.txt
     */
    String genotypeRecombinationMapFile;

    String taxaInfoFile;

    /**
     * Taxon\tPopulation\tPopulationName
     * tsk_0\t2\tC
     * tsk_1\t2\tC
     * tsk_2\t2\tC
     */
    Table<String, String, String> taxaInfoMap;

    String fastDirPath; // temp file for mosaic
    String logFile; // all log will be append to this logFile
    String outDir; // outDir
    IntList nWayAdmixture; // n way admixture
    int coresNumber; // core number, not threads number
    List<String> admixedPop; // admixed population name
    IntList chrIDList; // chrID of simulated genotypes
    List<File> genotypeList; // simulated genotype list
    List<File> recombinationFiles; // simulated recombination map list
    BiMap<File, File> genotypeRecombinationMap; // a bidirectional map between genotype and recombination map

    public Mosaic_runner(String parameterFile){
        try (BufferedReader brParameter = IOTool.getReader(parameterFile)) {
            String line;
            List<String> temp,tem;
            this.nWayAdmixture = new IntArrayList();
            this.admixedPop=new ArrayList<>();
            this.taxaInfoMap = HashBasedTable.create();
            while ((line=brParameter.readLine())!=null){
                temp = PStringUtils.fastSplit(line, ":");
                if (line.startsWith("GenotypeRecombinationMapPath")){
                    this.genotypeRecombinationMapFile = temp.get(1);
                    continue;
                }
                if (line.startsWith("TaxaInfoPath")){
                    this.taxaInfoFile = temp.get(1);
                    this.taxaInfoMap=RowTableTool.getTable(temp.get(1), 2);
                    continue;
                }
                if (line.startsWith("FastDirPath")){
                    this.fastDirPath=temp.get(1);
                    continue;
                }
                if (line.startsWith("LogFilePath")){
                    this.logFile=temp.get(1);
                    continue;
                }
                if (line.startsWith("OutDir")){
                    this.outDir=temp.get(1);
                    continue;
                }
                if (line.startsWith("nWayAdmixture")){
                    tem = PStringUtils.fastSplit(temp.get(1), ",");
                    for (String s : tem) {
                        this.nWayAdmixture.add(Integer.parseInt(s));
                    }
                    continue;
                }
                if (line.startsWith("CoresNumber")){
                    this.coresNumber=Integer.parseInt(temp.get(1));
                    continue;
                }
                if (line.startsWith("AdmixedPop")){
                    tem = PStringUtils.fastSplit(temp.get(1), ",");
                    this.admixedPop.addAll(tem);
                }
            }
            BufferedReader brMap = IOTool.getReader(this.genotypeRecombinationMapFile);
            brMap.readLine();
            this.chrIDList = new IntArrayList();
            this.genotypeList = new ArrayList<>();
            this.recombinationFiles = new ArrayList<>();
            this.genotypeRecombinationMap= HashBiMap.create();
            while ((line=brMap.readLine())!=null){
                temp =PStringUtils.fastSplit(line);
                this.chrIDList.add(Integer.parseInt(temp.get(0)));
                this.genotypeList.add(new File(temp.get(1)));
                this.recombinationFiles.add(new File(temp.get(2)));
                this.genotypeRecombinationMap.put(new File(temp.get(1)), new File(temp.get(2)));
            }
            brMap.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public int getNWayAdmixture(int index){
        return this.nWayAdmixture.getInt(index);
    }

    public String getAdmixedPopulation(int index){
        return this.admixedPop.get(index);
    }

    public int getChrID(int index){
        return this.chrIDList.getInt(index);
    }

    public File getGenotypeFile(int index){
        return this.genotypeList.get(index);
    }

    public File getRecombinationMap(int index){
        return this.recombinationFiles.get(index);
    }

    public int getSampleSizeOf(String population){
        int count = 0;
        for (Table.Cell<String,String,String> cell:this.taxaInfoMap.cellSet()){
            if (cell.getValue().equals(population)){
                count++;
            }
        }
        return count;
    }

    public int getAdmixedSampleSize(int index){
        return this.getSampleSizeOf(this.getAdmixedPopulation(index));
    }

    public double[][][] getLocalAncestry(){
        String[] outDirNames = this.genotypeList.stream().map(File::getName).map(s -> s.replaceAll(".vcf","")).toArray(String[]::new);
        File[] subDir_simulations = new File[outDirNames.length];
        for (int i = 0; i < outDirNames.length; i++) {
            subDir_simulations[i] = new File(this.outDir, outDirNames[i]);
            subDir_simulations[i].mkdir();
        }
        for (int i = 0; i < outDirNames.length; i++) {
            this.transformVCF_to_Mosaic_inputFormat(this.getGenotypeFile(i), this.getRecombinationMap(i), subDir_simulations[i]);
        }
        return this.runMosaicGetLocalAncestry(subDir_simulations);
    }

    private void transformVCF_to_Mosaic_inputFormat(File genotypeFile, File recombinationMap, File mosaicSimulationDir){
        String[] subDir = {this.subDirName};
        File[] subDirFile = new File[subDir.length];
        for (int i = 0; i < subDir.length; i++) {
            subDirFile[i] = new File(mosaicSimulationDir, subDir[i]);
            subDirFile[i].mkdir();
        }
        Map<String, String> taxonPopMap = RowTableTool.getMap(this.taxaInfoFile, 0, 2);
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
            Set<String> popSet = RowTableTool.getColumnSet(this.taxaInfoFile, 2);
            List<String> popList = new ArrayList<>(popSet);
            Collections.sort(popList);
            for (String value : popList) {
                popIndexMap.put(value, new IntArrayList());
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
            for (String s : popList) {
                genoFile = new File(subDirFile[0], s + "genofile.1");
                bw = IOTool.getWriter(genoFile);
                pop2BWMap.put(s, bw);
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

    /**
     *
     * @param subDir_simulations contain multiple simulation result dirs, such as different admixture proportion, generation
     *                     since admixture, divergence time between reference and source population et al. These
     *                     simulation results must be transformed to mosaic input format
     * @return local ancestry, dimension 1 is different simulation dirs, dimension 2 is taxon, dimension 3 is variants
     */
    private double[][][] runMosaicGetLocalAncestry(File[] subDir_simulations){
        String mosaicRpath = Simulation.class.getResource("mosaic.R").getPath();
        String extractMosaicLocalAncestryRpath = Simulation.class.getResource("extractMosaicLocalAncestry.R").getPath();
        double[][][] mosaicDir_taxa_variants_localAncestry = new double[subDir_simulations.length][][];
        // make sure logFile will be create from first call
        new File(this.logFile).delete();
        new File(this.logFile);
        IntStream.range(0, subDir_simulations.length).forEach(i->{
            File currentWorkingDir = subDir_simulations[i];
            File mosaicInputDataDir = currentWorkingDir.listFiles(dir -> dir.getName().equals(this.subDirName))[0];
            try {
                StringBuilder sb = new StringBuilder();
                List<String> command;
                ProcessBuilder processBuilder;
                Process process;
                int exitCode;

                // build and run mosaic command
                sb.setLength(0);
                sb.append("Rscript ").append(mosaicRpath).append(" ");
                sb.append(this.getAdmixedPopulation(i)).append(" ").append(mosaicInputDataDir.getPath()).append("/ ");
                sb.append("-f ").append(this.fastDirPath).append(" ");
                sb.append("-a ").append(this.getNWayAdmixture(i)).append(" ");
                sb.append("-m ").append(this.coresNumber).append(" ");
                sb.append("-c ").append(this.getChrID(i)).append(":").append(this.getChrID(i));
                command = PStringUtils.fastSplit(sb.toString(), " ");
                processBuilder = new ProcessBuilder(command);
                processBuilder.directory(currentWorkingDir);
                processBuilder.redirectErrorStream(true);
                processBuilder.redirectOutput(ProcessBuilder.Redirect.appendTo(new File(this.logFile)));
                process = processBuilder.start();
                exitCode = process.waitFor();
                assert exitCode == 0 : currentWorkingDir.getName()+" mosaic command run failed";

                // extract mosaic local ancestry result
                File mosaicResDir = currentWorkingDir.listFiles(dir -> dir.getName().startsWith("MOSAIC_RESULTS"))[0];
                File[] rDataFiles = mosaicResDir.listFiles(dir -> dir.getName().endsWith("RData"));
                String modelParametersRData = Arrays.stream(rDataFiles).filter(file -> !file.getName().startsWith("localanc")).toArray(File[]::new)[0].getAbsolutePath();
                String localAncestryRData = Arrays.stream(rDataFiles).filter(file -> file.getName().startsWith("localanc")).toArray(File[]::new)[0].getAbsolutePath();
                sb.setLength(0);
                sb.append("Rscript ").append(extractMosaicLocalAncestryRpath).append(" ");
                sb.append(modelParametersRData).append(" ").append(localAncestryRData).append(" ");
                sb.append(mosaicInputDataDir).append("/");
                command = PStringUtils.fastSplit(sb.toString(), " ");
                processBuilder = new ProcessBuilder(command);
                processBuilder.directory(currentWorkingDir);
                processBuilder.redirectError(ProcessBuilder.Redirect.appendTo(new File(this.logFile)));
                process = processBuilder.start();
                BufferedReader br = new BufferedReader(new InputStreamReader(process.getInputStream()));
                String line;
                List<String> temp;
                br.readLine();
                DoubleArrayList[] taxa_variants_localAncestry_list = new DoubleArrayList[this.getAdmixedSampleSize(i)];
                while ((line=br.readLine())!=null){
                    temp =PStringUtils.fastSplit(line);
                    for (int j = 0; j < temp.size(); j++) {
                        taxa_variants_localAncestry_list[j].add(Double.parseDouble(temp.get(j)));
                    }
                }
                br.close();
                exitCode = process.waitFor();
                assert exitCode == 0 : currentWorkingDir.getName()+" extract mosaic local ancestry failed";

                // assign double[][] to taxa_variant_localAncestry[][][]
                double[][] taxa_variant_localAncestry = new double[this.getAdmixedSampleSize(i)][];
                for (int j = 0; j < taxa_variants_localAncestry_list.length; j++) {
                    taxa_variant_localAncestry[j] = taxa_variants_localAncestry_list[j].toDoubleArray();
                }
                mosaicDir_taxa_variants_localAncestry[i] = taxa_variant_localAncestry;
            } catch (IOException | InterruptedException e) {
                e.printStackTrace();
            }

        });
        return mosaicDir_taxa_variants_localAncestry;
    }

    /**
     *
     * @param simulatedGenotypeDir input dir
     * @param simulatedRecombinationMapDir input dir
     * @param parameterFile output file
     * @param genotypeRecombinationMap output file
     * @param chrIDArray chrID
     * @param taxaInfoFile
     * @param fastDirPath
     * @param logFile
     * @param outDir
     * @param nWayAdmixture
     * @param coresNumber
     * @param admixedPop
     */
    public static void writeParameterFile(String simulatedGenotypeDir, String simulatedRecombinationMapDir,
                                          String parameterFile, String genotypeRecombinationMap, int[] chrIDArray,
                                          String taxaInfoFile, String fastDirPath, String logFile, String outDir,
                                          int[] nWayAdmixture, int coresNumber, String[] admixedPop){
        List<File> genotypeFiles = IOTool.getFileListInDirEndsWith(simulatedGenotypeDir, "Ne.vcf");
        List<File> recombinationFiles = IOTool.getFileListInDirEndsWith(simulatedRecombinationMapDir, "Map.txt");
        try (BufferedWriter bwParameter = IOTool.getWriter(parameterFile);
             BufferedWriter bwMap = IOTool.getWriter(genotypeRecombinationMap)) {
            StringBuilder sb = new StringBuilder();
            bwMap.write("ChromosomeID\tGenotype\tRecombinationMap\n");
            for (int i = 0; i < genotypeFiles.size(); i++) {
                sb.setLength(0);
                sb.append(chrIDArray[i]).append("\t");
                sb.append(genotypeFiles.get(i).getAbsolutePath()).append("\t");
                sb.append(recombinationFiles.get(i).getAbsolutePath());
                bwMap.write(sb.toString());
                bwMap.newLine();
            }
            sb.setLength(0);
            sb.append("## mosaic runner parameter file\n");
            sb.append("## File path\n");
            sb.append("GenotypeRecombinationMapPath:").append(genotypeRecombinationMap).append("\n");
            sb.append("TaxaInfoPath:").append(taxaInfoFile).append("\n");
            sb.append("FastDirPath:").append(fastDirPath).append("\n");
            sb.append("LogFilePath:").append(logFile).append("\n");
            sb.append("OutDir:").append(outDir).append("\n");
            sb.append("\n");
            sb.append("## Auguments of mosaic\n");
            sb.append("nWayAdmixture:").append(Ints.join(",", nWayAdmixture)).append("\n");
            sb.append("CoresNumber:").append(coresNumber).append("\n");
            sb.append("AdmixedPop:").append(String.join(",", admixedPop)).append("\n");
            bwParameter.write(sb.toString());
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


}
