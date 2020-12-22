package daxing.load.slidingWindow;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import daxing.common.IOTool;
import daxing.common.NumberTool;
import daxing.common.RowTableTool;
import daxing.common.Triads;
import daxing.load.ancestralSite.ChrSNPAnnoDB;
import daxing.load.ancestralSite.Standardization;
import daxing.load.neutralSiteLoad.DynamicSNPGenotypeDB;
import daxing.load.neutralSiteLoad.IndividualChrLoad;
import daxing.load.neutralSiteLoad.SNPGenotype;
import pgl.infra.table.RowTable;
import pgl.infra.utils.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class SlidingWindowGo {

    public static void getDerivedCount(String exonAnnoDir,String vcfDir, String ancestralDir, String pgfFile,
                                       String triadGeneNameFile, String outdir, String vmapIIGroupFile){
        RowTableTool<String> table=new RowTableTool<>(triadGeneNameFile);
        List<String> geneNames=table.getColumn(2);
        GeneWindowDB geneWindowDB =new GeneWindowDB(geneNames, pgfFile);
        List<File> exonAnnoFiles=IOUtils.getVisibleFileListInDir(exonAnnoDir);
        List<File> vcfFiles= IOUtils.getVisibleFileListInDir(vcfDir);
        List<File> ancestralFiles=IOUtils.getVisibleFileListInDir(ancestralDir);
        Map<String,File> vmapIITaxonoutDirMap=getTaxonOutDirMap(vmapIIGroupFile, outdir);
        int[][] indices= PArrayUtils.getSubsetsIndicesBySubsetSize(vcfFiles.size(), 3);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> integerList= Arrays.asList(subLibIndices);
            integerList.parallelStream().forEach(e-> getDerivedCount(exonAnnoFiles.get(e),vcfFiles.get(e),
                    ancestralFiles.get(e), geneWindowDB.getGeneWindow(e+1), vmapIITaxonoutDirMap));
        }
    }

    private static void getDerivedCount(File exonAnnoFile ,File vcfFile, File ancestralFile,
                                        List<GeneWindowDB.GeneWindow> geneWindows, Map<String,File> taxonOutDirMap){
        long start=System.nanoTime();
        Table<Integer,Integer, Character> ancestralTable=getAncestral(ancestralFile);
        String[] geneName=new String[geneWindows.size()];
        for (int i = 0; i < geneName.length; i++) {
            geneName[i]=geneWindows.get(i).getGeneName();
        }
        try (BufferedReader br = IOTool.getReader(vcfFile)) {
            String line;
            List<String> temp, taxonNames;
            while ((line=br.readLine()).startsWith("##")){}
            temp=PStringUtils.fastSplit(line);
            taxonNames=temp.subList(9, temp.size());
            IndividualChrLoad[] individualChrLoads=new IndividualChrLoad[taxonNames.size()];
            for (int i = 0; i < individualChrLoads.length; i++) {
                individualChrLoads[i]=new IndividualChrLoad(taxonNames.get(i), geneName);
            }
            int chr = -1, pos;
            char refBase, altBase, ancestralBase;
            DynamicSNPGenotypeDB dynamicSNPGenotypeDB=new DynamicSNPGenotypeDB();
            SNPGenotype snpGenotype=null;
            ChrSNPAnnoDB chrSNPAnnoDB=new ChrSNPAnnoDB(exonAnnoFile);
            int[][] derivedCount, nonsynDerivedCount;
            String outFileName=null, variantsType;
            for (int i = 0; i < geneWindows.size(); i++) {
                while ((line=br.readLine())!=null){
                    temp= PStringUtils.fastSplit(line);
                    chr=Integer.parseInt(temp.get(0));
                    pos=Integer.parseInt(temp.get(1));
                    if (!ancestralTable.contains(chr, pos)) continue;
                    refBase=temp.get(3).charAt(0);
                    altBase=temp.get(4).charAt(0);
                    ancestralBase=ancestralTable.get(chr, pos);
                    if ((ancestralBase!=refBase) && (ancestralBase!=altBase)) continue;
                    if (pos < geneWindows.get(i).getStartOfWindowCDSRange()) continue;
                    if (geneWindows.get(i).containCDSPos(chr, pos)){
                        variantsType=chrSNPAnnoDB.getVariantType(chr, pos);
                        dynamicSNPGenotypeDB.addSNPGenotype(SNPGenotype.getSNPGenotype(line, ancestralBase, variantsType,temp.subList(9, temp.size())));
                    }else {
                        variantsType=chrSNPAnnoDB.getVariantType(chr, pos);
                        snpGenotype=SNPGenotype.getSNPGenotype(line, ancestralBase, variantsType, temp.subList(9, temp.size()));
                        break;
                    }
                }
                if (dynamicSNPGenotypeDB.size()==0) continue;
                derivedCount=dynamicSNPGenotypeDB.countAllTaxonDerived();
                nonsynDerivedCount=dynamicSNPGenotypeDB.countAllTaxonNonsynDerived();
                for (int j = 0; j < derivedCount.length; j++) {
                    individualChrLoads[j].addGeneDerivedCount(i, derivedCount[j], nonsynDerivedCount[j]);
                }
                if (i<=geneWindows.size()-2){
                    dynamicSNPGenotypeDB.retainAll(geneWindows.get(i+1));
                    if (geneWindows.get(i+1).containCDSPos(snpGenotype.getChromosome(), snpGenotype.getPosition())){
                        dynamicSNPGenotypeDB.addSNPGenotype(snpGenotype);
                    }
                }
            }
            for (int j = 0; j < individualChrLoads.length; j++) {
                outFileName="chr"+PStringUtils.getNDigitNumber(3, chr)+"."+taxonNames.get(j)+".txt.gz";
                individualChrLoads[j].write(taxonOutDirMap.get(taxonNames.get(j)), outFileName);
            }
            System.out.println("chr"+PStringUtils.getNDigitNumber(3, chr)+" completed in "+ Benchmark.getTimeSpanMinutes(start)+" minutes");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static Table<Integer, Integer, Character> getAncestral(File ancestralFile){
        Table<Integer, Integer, Character> ancestralTable= HashBasedTable.create();
        try (BufferedReader br = IOTool.getReader(ancestralFile)) {
            br.readLine();
            String line;
            List<String> temp;
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                ancestralTable.put(Integer.parseInt(temp.get(0)), Integer.parseInt(temp.get(1)), temp.get(2).charAt(0));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return ancestralTable;
    }

    private static Map<String, File> getTaxonOutDirMap(String vmapIIGroupFile, String outDir){
        Map<String,File> taxonOutDirMap=new HashMap<>();
        RowTable<String> taxonTable=new RowTable<>(vmapIIGroupFile);
        List<String> taxonNames=taxonTable.getColumn(0);
        File taxonDir;
        for (int i = 0; i < taxonNames.size(); i++) {
            taxonDir=new File(outDir, taxonNames.get(i));
            taxonDir.mkdir();
            taxonOutDirMap.put(taxonNames.get(i), taxonDir);
        }
        return taxonOutDirMap;
    }

    public static void merge(String inputDir, String triadFile, String outDir){
        List<File> dirs=IOTool.getVisibleDir(inputDir);
        int[][] indices= PArrayUtils.getSubsetsIndicesBySubsetSize(dirs.size(), 21);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> integerList=Arrays.asList(subLibIndices);
            integerList.parallelStream().forEach(e->merge(dirs.get(e), triadFile, outDir));
        }
    }

    private static void merge(File inputDir, String triadFile, String outDir){
        Triads triads =new Triads(triadFile);
        List<File> files= IOUtils.getVisibleFileListInDir(inputDir.getAbsolutePath());
        RowTableTool<String> tableTool=new RowTableTool<>(files.get(0).getAbsolutePath());
        String outName= PStringUtils.fastSplit(files.get(0).getName(), ".").get(1);
        for (int j = 1; j < files.size(); j++) {
            tableTool.add(new RowTableTool<>(files.get(j).getAbsolutePath()));
        }
        List<String> triadIDList=new ArrayList<>();
        String geneName, triadID;
        for (int j = 0; j < tableTool.getRowNumber(); j++) {
            geneName=tableTool.getCell(j, 0);
            triadID= triads.getTraidID(geneName);
            triadIDList.add(triadID);
        }
        tableTool.insertColumn("TriadID", 0, triadIDList);
        Comparator<List<String>> c=Comparator.comparing(l->l.get(0));
        Comparator<List<String>> cc=c.thenComparing(l->l.get(1).substring(8,9));
        tableTool.sortBy(cc);
        tableTool.write(new File(outDir, outName+".txt.gz"), IOFileFormat.TextGzip);
    }

    public static void start(){
        String retainTriadDir="/Users/xudaxing/Documents/deleteriousMutation/002_analysis/014_deleterious/004_analysis/003_geneLoadIndividual_GeneHC/004_standardization/002_standizationByGeneLocal_100kb/002_addNeutralSiteLoad";
        String neutralSiteLoad="/Users/xudaxing/Documents/deleteriousMutation/002_analysis/014_deleterious/004_analysis/003_geneLoadIndividual_GeneHC/004_standardization/003_standizationByGeneLoacal_10Genes/001_individualLoad";
        String vmapIIGroupFile="/Users/xudaxing/Documents/deleteriousMutation/002_analysis/014_deleterious/vmapGroup.txt";
        String[] dirs={"002_addNeutralSiteLoad","003_standization"};
        String parentDir=new File(neutralSiteLoad).getParent();
        for (int i = 0; i < dirs.length; i++) {
            new File(parentDir, dirs[i]).mkdir();
        }
        insertColumnNeutralSiteLoad(retainTriadDir, neutralSiteLoad, vmapIIGroupFile, new File(parentDir, dirs[0]));
//        normalized(new File(parentDir, dirs[0]).getAbsolutePath(), new File(parentDir, dirs[1]).getAbsolutePath());
    }

    private static void insertColumnNeutralSiteLoad(String retainTriadDir, String neutralSiteLoadDir,
                                                    String vmapIIGroupFile, File outDir){
        List<File> neutralFiles=IOUtils.getVisibleFileListInDir(neutralSiteLoadDir);
        RowTableTool<String> vmapIITable=new RowTableTool<>(vmapIIGroupFile);
        Map<String,String> taxonGroupMap=vmapIITable.getHashMap(0,15);
        Predicate<File> landraceCultivar= f->taxonGroupMap.get(PStringUtils.fastSplit(f.getName(),".").get(0)).equals(
                "Landrace_Europe") || taxonGroupMap.get(PStringUtils.fastSplit(f.getName(),".").get(0)).equals("Cultivar");
        List<File> landraceCultivarFiles=neutralFiles.stream().filter(landraceCultivar).sorted().collect(Collectors.toList());
        List<File> files=IOUtils.getVisibleFileListInDir(retainTriadDir);
        String[] outNames= files.stream().map(File::getName).map(s->s.replaceAll(".txt.gz", ".neutralLoad.txt.gz")).toArray(String[]::new);
        (IntStream.range(0, files.size())).parallel().forEach(e->insertColumnNeutralSiteLoad(files.get(e),
                landraceCultivarFiles.get(e), new File(outDir, outNames[e])));
    }

    private static void insertColumnNeutralSiteLoad(File retainTriadFile, File neutralSiteLoadFile, File outFile){
        RowTableTool<String> neutralTable=new RowTableTool<>(neutralSiteLoadFile.getAbsolutePath());
        RowTableTool<String> retainTable=new RowTableTool<>(retainTriadFile.getAbsolutePath());
        Comparator<List<String>> c=Comparator.comparing(l->l.get(0));
        Comparator<List<String>> cc=c.thenComparing(l->l.get(1).substring(8,9));
        neutralTable.sortBy(cc);
        List<String> geneNameNeutral=neutralTable.getColumn(1);
        List<String> geneNameRetain=retainTable.getColumn(2);
        if (!geneNameNeutral.equals(geneNameRetain)){
            System.out.println("neutral geneNameList not equal retainTriadFileList name");
            System.exit(1);
        }
        List<String> numGeneLocal=neutralTable.getColumn(2);
        List<String> numDerivedInGeneLocal=neutralTable.getColumn(3);
        List<String> numNonsynGeneLocal=neutralTable.getColumn(4);
        List<String> numNonsynDerivedInGeneLocal=neutralTable.getColumn(5);
        retainTable.addColumn("numGeneWindow", numGeneLocal);
        retainTable.addColumn("numDerivedInGeneWindow", numDerivedInGeneLocal);
        retainTable.addColumn("numNonsynGeneWindow", numNonsynGeneLocal);
        retainTable.addColumn("numNonsynDerivedInGeneWindow",numNonsynDerivedInGeneLocal);
        retainTable.write(outFile,IOFileFormat.TextGzip);
    }

    private static void filter(String triadNeutralLoadFileInputDir, String outDir, int snpNumThresh, int neutralNumThresh){
        List<File> files=IOUtils.getVisibleFileListInDir(triadNeutralLoadFileInputDir);
        String[] outFileName= files.stream().map(File::getName).map(s->s.replaceAll(".txt.gz",".filtered.txt.gz")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->filter(files.get(e),new File(outDir, outFileName[e]), snpNumThresh, neutralNumThresh));
    }

    private static void filter(File triadNeutralLoadFileInputFile, File outFile, int snpNumThresh,
                               int neutralNumThresh) {
        try (BufferedReader br = IOTool.getReader(triadNeutralLoadFileInputFile);
             BufferedWriter bw = IOTool.getWriter(outFile)) {
            String header = br.readLine();
            bw.write(header);
            bw.newLine();
            String line;
            String[] lineABD;
            int[] snpNum, neutralNum;
            List<String> temp;
            int triadNum = 0;
            int retainedNum = 0;
            while ((line = br.readLine()) != null) {
                lineABD = new String[3];
                lineABD[0] = line;
                lineABD[1] = br.readLine();
                lineABD[2] = br.readLine();
                triadNum++;
                snpNum = new int[3];
                neutralNum = new int[3];
                for (int i = 0; i < lineABD.length; i++) {
                    temp = PStringUtils.fastSplit(lineABD[i]);
                    snpNum[i] = Integer.parseInt(temp.get(3)) + Integer.parseInt(temp.get(5));
                    neutralNum[i] = Integer.parseInt(temp.get(11)) - Integer.parseInt(temp.get(13));
                }
                if (snpNum[0] < snpNumThresh) continue;
                if (snpNum[1] < snpNumThresh) continue;
                if (snpNum[2] < snpNumThresh) continue;
                if (neutralNum[0] < neutralNumThresh) continue;
                if (neutralNum[1] < neutralNumThresh) continue;
                if (neutralNum[2] < neutralNumThresh) continue;
                retainedNum++;
                for (int i = 0; i < lineABD.length; i++) {
                    bw.write(lineABD[i]);
                    bw.newLine();
                }
            }
            System.out.println(triadNeutralLoadFileInputFile.getName() + ":\t" + triadNum + "\t" + retainedNum + "\t" + (double) retainedNum / triadNum);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void normalized(String inputDir, String outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir);
        String[] outNames= files.stream().map(File::getName).map(s->s.replaceAll(".txt.gz",".normalized.txt.gz")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->normalized(files.get(e), new File(outDir, outNames[e])));
    }

    private static void normalized(File inputFile, File outFile){
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            bw.write("TriadID\tnormalizedNumDerivedInNonsynA" +
                    "\tnormalizedNumDerivedInNonsynB\tnormalizedNumDerivedInNonsynD\tnonsynRegion" +
                    "\tnormalizedNumDerivedInHGDeleteriousA\tnormalizedNumDerivedInHGDeleteriousB" +
                    "\tnormalizedNumDerivedInHGDeleteriousD\tdelRegion");
            bw.newLine();
            br.readLine();
            String line, triadID, region;
            List<String>[] temp;
            int[] cdsLen;
            int[] neutralNum;
            double[] derivedNonsyn;
            double[] deleterious;
            StringBuilder sb;
            double[] normalizedDerivedNonSyn;
            double[] normalizedDerivedDel;
            while ((line=br.readLine())!=null){
                temp=new List[3];
                temp[0]=PStringUtils.fastSplit(line);
                temp[1]=PStringUtils.fastSplit(br.readLine());
                temp[2]=PStringUtils.fastSplit(br.readLine());
                cdsLen=new int[3];
                neutralNum=new int[3];
                derivedNonsyn=new double[3];
                deleterious=new double[3];
                normalizedDerivedNonSyn=new double[3];
                normalizedDerivedDel=new double[3];
                triadID=temp[0].get(0);
                for (int i = 0; i < cdsLen.length; i++) {
                    cdsLen[i]=Integer.parseInt(temp[i].get(1));
                    neutralNum[i]=Integer.parseInt(temp[i].get(11))-Integer.parseInt(temp[i].get(13));
                    derivedNonsyn[i]=Double.parseDouble(temp[i].get(6));
                    deleterious[i]=Double.parseDouble(temp[i].get(8));
                }
                for (int i = 0; i < 3; i++) {
                    normalizedDerivedNonSyn[i]=1000*100*derivedNonsyn[i]/(cdsLen[i]*neutralNum[i]);
                    normalizedDerivedDel[i]=1000*100*deleterious[i]/(cdsLen[i]*neutralNum[i]);
                }
                sb=new StringBuilder();
                sb.append(triadID).append("\t");
                for (int i = 0; i < 3; i++) {
                    sb.append(NumberTool.format(normalizedDerivedNonSyn[i], 5)).append("\t");
                }
                region= Standardization.getNearestPointIndex(normalizedDerivedNonSyn).getRegion();
                sb.append(region).append("\t");
                for (int i = 0; i < 3; i++) {
                    sb.append(NumberTool.format(normalizedDerivedDel[i], 5)).append("\t");
                }
                region=Standardization.getNearestPointIndex(normalizedDerivedDel).getRegion();
                sb.append(region);
                bw.write(sb.toString());
                bw.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
