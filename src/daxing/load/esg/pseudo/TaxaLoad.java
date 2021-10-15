package daxing.load.esg.pseudo;

import daxing.common.factors.WheatLineage;
import daxing.common.table.RowTableTool;
import daxing.common.utiles.ArrayTool;
import daxing.common.utiles.IOTool;
import daxing.common.factors.Ploidy;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;
import java.util.function.Predicate;

/**
 * 用于真实六(四)倍体的ABD亚基因组合成伪六(四)倍体
 */
public class TaxaLoad {

    List<DerivedCountInfo>[][] taxaDerivedCountInfo; //dim1: taxa, dim2: subgenome, List: gene
    String[] taxa;
    Ploidy ploidy;

    /**
     * landrace cultivar
     * @param allTaxaDerivedCountDir
     * @param germplasmInfoFile
     */
    public TaxaLoad(String allTaxaDerivedCountDir, String germplasmInfoFile, Ploidy ploidy){
        long start= System.nanoTime();
        List<File> fileList = IOTool.getVisibleFileListInDir(allTaxaDerivedCountDir);
        Predicate<File> p = this.getPseudoHexaploidPredict(germplasmInfoFile, ploidy);
        File[] files= fileList.stream().filter(p).toArray(File[]::new);
        this.buildFromDir(files, ploidy);
        System.out.println("Reading spent "+ Benchmark.getTimeSpanSeconds(start)+ " seconds");
    }

    private Predicate<File> getPseudoHexaploidPredict(String germplasmInfoFile, Ploidy ploidy){
        Map<String, String> taxonSubspecies_by6_TreeValidatedMap= RowTableTool.getMap(germplasmInfoFile,0, 10);
        Predicate<File> pLandrace =
                file -> taxonSubspecies_by6_TreeValidatedMap.get(PStringUtils.fastSplit(file.getName(),
                        ".").get(0)).equals("Landrace");
        Predicate<File> pCultivar =
                file -> taxonSubspecies_by6_TreeValidatedMap.get(PStringUtils.fastSplit(file.getName(),
                        ".").get(0)).equals("Cultivar");
        Predicate<File> pWE=file -> taxonSubspecies_by6_TreeValidatedMap.get(PStringUtils.fastSplit(file.getName(),
                ".").get(0)).equals("Wild_emmer");
        Predicate<File> pDE=file -> taxonSubspecies_by6_TreeValidatedMap.get(PStringUtils.fastSplit(file.getName(),
                ".").get(0)).equals("Domesticated_emmer");
        Predicate<File> pFTT=file -> taxonSubspecies_by6_TreeValidatedMap.get(PStringUtils.fastSplit(file.getName(),
                ".").get(0)).equals("Free_threshing_tetraploids");
        Predicate<File> pLandraceCultivar = pLandrace.or(pCultivar);
        Predicate<File> pWEDEFTT = pWE.or(pDE.or(pFTT));
        if (ploidy.equals(Ploidy.HEXAPLOID)){
            return pLandraceCultivar;
        }else if (ploidy.equals(Ploidy.tetraploid)){
            return pWEDEFTT;
        }
        return null;
    }

    private void buildFromDir(File[] allTaxaDerivedCountFiles, Ploidy ploidy){
        List<DerivedCountInfo>[][] derivedCountInfoList = new List[allTaxaDerivedCountFiles.length][];
        String[] taxa= new String[allTaxaDerivedCountFiles.length];
        for (int i = 0; i < derivedCountInfoList.length; i++) {
            derivedCountInfoList[i] = new List[ploidy.getSubgenomewNum()];
            for (int j = 0; j < derivedCountInfoList[i].length; j++) {
                derivedCountInfoList[i][j]= new ArrayList<>();
            }
        }
        try {
            BufferedReader br;
            String line;
            List<String> temp;
            int[] derivedCount;
            String geneName;
            WheatLineage wheatLineage;
            for (int i = 0; i < allTaxaDerivedCountFiles.length; i++) {
                taxa[i]=PStringUtils.fastSplit(allTaxaDerivedCountFiles[i].getName(), ".").get(0);;
                br = IOTool.getReader(allTaxaDerivedCountFiles[i]);
                br.readLine();
                while ((line = br.readLine())!=null){
                    temp = PStringUtils.fastSplit(line);
                    geneName = temp.get(0);
                    wheatLineage = WheatLineage.valueOf(temp.get(0).substring(8,9));
                    derivedCount = new int[9];
                    for (int j = 1; j < temp.size(); j++) {
                        derivedCount[j-1] = Integer.parseInt(temp.get(j));
                    }
                    derivedCountInfoList[i][wheatLineage.getIndex()].add(new DerivedCountInfo(geneName, derivedCount));
                }
                br.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        this.ploidy= ploidy;
        this.taxa=taxa;
        this.taxaDerivedCountInfo=derivedCountInfoList;
    }

    public List<DerivedCountInfo>[][] getTaxaDerivedCountInfo() {
        return taxaDerivedCountInfo;
    }

    public Ploidy getPloidy() {
        return ploidy;
    }

    public String[] getTaxa() {
        return taxa;
    }

    public int getTaxaNum(){
        return this.taxa.length;
    }

    public String getTaxon(int taxonIndex){
        return this.taxa[taxonIndex];
    }

    /**
     * 随机组合六倍体个体的三个亚基因组，形成伪六倍体
     */
    public void writePseudoHexaploid(int pseudoNum, String outDir){
        long start =System.nanoTime();
        int[][] randomABDTaxonIndex= new int[this.getPloidy().getSubgenomewNum()][]; // dim1: subgenome dim2: random number
        for (int i = 0; i < randomABDTaxonIndex.length; i++) {
            randomABDTaxonIndex[i]= ArrayTool.getRandomNumberArray(pseudoNum, 0, this.getTaxaNum());
        }
        try {
            BufferedWriter bw;
            String combinedTaxon;
            String[] taxonABD;
            String header="geneName\tnumSyn\tnumDerivedInSyn\tnumHeterInSyn\tnumNonsyn\tnumDerivedInNonsyn" +
                    "\tnumHeterInNonsyn\tnumHGDeleterious\tnumDerivedInHGDeleterious\tnumHeterInHGDeleterious";
            List<DerivedCountInfo> individualLoad;
            StringBuilder sb = new StringBuilder();
            int count=0;
            for (int i = 0; i < randomABDTaxonIndex[0].length; i++) {
                taxonABD = new String[this.getPloidy().getSubgenomewNum()];
                for (int j = 0; j < taxonABD.length; j++) {
                    taxonABD[j] = this.getTaxon(randomABDTaxonIndex[j][i]);
                    assert taxonABD[j] != null:"exist null";
                }
                combinedTaxon = String.join("_", taxonABD);
                bw = IOTool.getWriter(new File(outDir, combinedTaxon+".txt.gz"));
                bw.write(header);
                bw.newLine();
                individualLoad = new ArrayList<>();
                for (int j = 0; j < randomABDTaxonIndex.length; j++) {
                    individualLoad.addAll(this.getTaxaDerivedCountInfo()[randomABDTaxonIndex[j][i]][j]);
                }
                Collections.sort(individualLoad, Comparator.comparing(d -> d.geneName));
                for (DerivedCountInfo derivedCountInfo1 : individualLoad){
                    sb.setLength(0);
                    sb.append(derivedCountInfo1.toString());
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                count++;
                if (count%100==0){
                    System.out.println(i+" files had been written to "+ outDir);
                }
            }
            System.out.println("Total "+ count+ " files had been written to "+outDir);
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Writing spend "+Benchmark.getTimeSpanSeconds(start)+ " seconds");
    }

    public static class DerivedCountInfo {

        String geneName;
        /**
         * numSyn	numDerivedInSyn	numHeterInSyn
         * numNonsyn	numDerivedInNonsyn	numHeterInNonsyn
         * numHGDeleterious	numDerivedInHGDeleterious	numHeterInHGDeleterious
         */
        int[] derivedInfo;

        public DerivedCountInfo(String geneName, int[] derivedInfo){
            this.geneName= geneName;
            this.derivedInfo=derivedInfo;
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append(geneName).append("\t");
            for (int i = 0; i <this.derivedInfo.length; i++) {
                sb.append(derivedInfo[i]).append("\t");
            }
            sb.deleteCharAt(sb.length()-1);
            return sb.toString();
        }
    }
}
