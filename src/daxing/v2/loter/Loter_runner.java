package daxing.v2.loter;

import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;
import com.google.common.collect.Table;
import daxing.common.sh.CommandUtils;
import daxing.common.table.RowTableTool;
import daxing.common.utiles.IOTool;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class Loter_runner {

    BiMap<String, String> genotypeID2PathMap;
    String taxaInfoFile;

    /**
     * Taxon\tPopulation\tPopulationName
     * tsk_0\t2\tC
     * tsk_1\t2\tC
     * tsk_2\t2\tC
     */
    Table<String, String, String> taxaInfoMap;

    String logFile; // all log will be append to this logFile

    String outDir; // outDir

    String admixedPop; // must be one
    List<String> referencePopList; // at least two

    String[] subDir = {"loterGroup", "loterGroupVCF", "LAI"};
    File[] subDirFile;

    public Loter_runner(String parameterFile){
        this.referencePopList = new ArrayList<>();
        this.genotypeID2PathMap = HashBiMap.create();
        try (BufferedReader brParameter = IOTool.getReader(parameterFile)) {
            String line, lineGenotype;
            List<String> temp, tem;
            BufferedReader brGenotype;
            while ((line=brParameter.readLine())!=null){
                temp = PStringUtils.fastSplit(line, ":");
                if (line.startsWith("##")){}
                if (line.startsWith("GenotypePath")){
                    brGenotype = IOTool.getReader(temp.get(1));
                    brGenotype.readLine();
                    while ((lineGenotype = brGenotype.readLine())!=null){
                        tem = PStringUtils.fastSplit(lineGenotype);
                        this.genotypeID2PathMap.put(tem.get(0), tem.get(1));
                    }
                    brGenotype.close();
                    continue;
                }
                if (line.startsWith("TaxaInfoPath")){
                    this.taxaInfoFile = temp.get(1);
                    this.taxaInfoMap= RowTableTool.getTable(temp.get(1), 2);
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
                if (line.startsWith("admixedPop")){
                    this.admixedPop = temp.get(1);
                    continue;
                }
                if (line.startsWith("referencePop")){
                    tem = PStringUtils.fastSplit(temp.get(1), ",");
                    this.referencePopList.addAll(tem);
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        this.subDirFile = new File[this.subDir.length];
        for (int i = 0; i < this.subDir.length; i++) {
            this.subDirFile[i] = new File(this.outDir, this.subDir[i]);
            this.subDirFile[i].mkdir();
        }
//        this.preparePop();
//        this.splitPopGenotype();
//        this.run_loter();
    }

    /**
     *
     * @return admixed pop and reference pop
     */
    public List<String> getAllPop(){
        List<String> allPop = new ArrayList<>(this.referencePopList);
        allPop.add(this.admixedPop);
        return allPop;
    }

    public int getSampleSize(String popName){
        List<String> allPop = this.getAllPop();
        assert allPop.contains(popName) : "check popName";
        int count=0;
        for (Table.Cell<String,String,String> cell : this.taxaInfoMap.cellSet()){
            if (cell.getValue().equals(popName)){
                count++;
            }
        }
        return count;
    }

    public int getAdmixedPopSampleSize(){
        return this.getSampleSize(this.admixedPop);
    }

    public int getGenotypeCount(){
        return this.genotypeID2PathMap.size();
    }

    private void preparePop(){
        Set<String> popSet = RowTableTool.getColumnSet(this.taxaInfoFile, 2);
        Map<String, List<String>> pop2TaxaMap = new HashMap<>();
        Map<String, BufferedWriter> pop2BW = new HashMap<>();
        for (String pop : popSet){
            pop2TaxaMap.put(pop, new ArrayList<>());
        }
        for (Table.Cell<String,String,String> cell: this.taxaInfoMap.cellSet()){
            pop2TaxaMap.get(cell.getValue()).add(cell.getRowKey());
            pop2BW.put(cell.getValue(), IOTool.getWriter(new File(this.subDirFile[0], cell.getValue()+".txt")));
        }
        BufferedWriter bw;
        try {
            for (Map.Entry<String, List<String>> entry : pop2TaxaMap.entrySet()){
                bw = pop2BW.get(entry.getKey());
                for (String line : entry.getValue()){
                    bw.write(line);
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private void splitPopGenotype(){
        StringBuilder sb = new StringBuilder();
        List<String> commandList = new ArrayList<>();
        for (Map.Entry<String, String> entry : this.genotypeID2PathMap.entrySet()){
            for (String pop : this.getAllPop()){
                sb.setLength(0);
                sb.append("vcftools --gzvcf ").append(entry.getValue()).append(" ");
                sb.append("--keep ").append(new File(this.subDirFile[0], pop+".txt").getAbsolutePath()).append(" ");
                sb.append("--recode --out ").append(new File(this.subDirFile[1], entry.getKey()+"_"+pop).getAbsolutePath());
                commandList.add(sb.toString());
            }
        }
        int exitCode;
        for (String command: commandList){
            exitCode = CommandUtils.runOneCommand(command, this.outDir, new File(this.logFile));
            assert exitCode == 0 : command+" run failed";
        }
    }

    private void run_loter(){
        this.genotypeID2PathMap.entrySet().stream().forEach(entry->{
            StringBuilder sb = new StringBuilder();
            StringBuilder sbRef = new StringBuilder();
            sb.setLength(0);
            sb.append("source /Users/xudaxing/anaconda3/bin/activate loter && ");
            sb.append("loter_cli -r ");
            sbRef.setLength(0);
            for (String refPop : this.referencePopList){
                sbRef.append(new File(this.subDirFile[1], entry.getKey()+"_"+refPop+".recode.vcf").getAbsolutePath());
                sbRef.append(" ");
            }
            sb.append(sbRef).append(" -a ").append(new File(this.subDirFile[1], entry.getKey()+"_"+this.admixedPop+".recode.vcf").getAbsolutePath());
            sb.append(" -f vcf -o ").append(new File(this.subDirFile[2], entry.getKey()+".lai.txt"));
            int exitCode = CommandUtils.runMultipleCommands(sb.toString(), this.outDir, new File(this.logFile));
            assert exitCode == 0 : sb+" run failed";
        });
    }

    public double[][][] extractLocalAncestry(){
        List<File> laiFiles = new ArrayList<>();
        StringBuilder sb = new StringBuilder();
        File laiFile;
        for (Map.Entry<String, String> entry:this.genotypeID2PathMap.entrySet()){
            sb.setLength(0);
            sb.append(entry.getKey()).append(".lai.txt");
            laiFile = new File(this.subDirFile[2], sb.toString());
            laiFiles.add(laiFile);
        }
        double[][][] genotype_taxa_variants_localAncestry = new double[this.getGenotypeCount()][][];
        for (int i = 0; i < laiFiles.size(); i++) {
            genotype_taxa_variants_localAncestry[i] = this.extractLocalAncestry(laiFiles.get(i));
        }
        return genotype_taxa_variants_localAncestry;
    }

    private double[][] extractLocalAncestry(File laiFile){
        String line;
        List<String> temp;
        DoubleArrayList[] localAncestry = new DoubleArrayList[this.getAdmixedPopSampleSize()];
        for (int i = 0; i < localAncestry.length; i++) {
            localAncestry[i] = new DoubleArrayList();
        }
        int taxaIndex=0;
        try (BufferedReader br = IOTool.getReader(laiFile)) {
            while ((line=br.readLine())!=null){
                temp =PStringUtils.fastSplit(line, " ");
                for (int i = 0; i < temp.size(); i++) {
                    localAncestry[taxaIndex].add(Integer.parseInt(temp.get(i)) == 0 ? 0 : 4);
                }
                br.readLine();
                taxaIndex++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        double[][] localAnc = new double[this.getAdmixedPopSampleSize()][];
        for (int i = 0; i < localAnc.length; i++) {
            localAnc[i] = localAncestry[i].toDoubleArray();
        }
        return localAnc;
    }


}