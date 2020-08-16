package daxing.applets;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Multimap;
import com.google.common.collect.Table;
import com.google.common.collect.TreeMultimap;
import daxing.common.*;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.CombinatoricsUtils;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ScriptMethods {

    public static void getTopRows(File inputFile, int n, File outputFile){
        try{
            long start=System.nanoTime();
            BufferedReader br;
            BufferedWriter bw;
            if (inputFile.getName().endsWith("gz")){
                br=IOUtils.getTextGzipReader(inputFile.getAbsolutePath());
            }else {
                br=IOUtils.getTextReader(inputFile.getAbsolutePath());
            }
            if (outputFile.getName().endsWith(".gz")){
                bw=IOUtils.getTextGzipWriter(outputFile.getAbsolutePath());
            }
            else {
                bw=IOUtils.getTextWriter(outputFile.getAbsolutePath());
            }
            String line;
            int count=1;
            while ((line=br.readLine())!=null){
                if (count > n) break;
                bw.write(line);
                bw.newLine();
                count++;
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(outputFile.getName()+" is completed in "+Benchmark.getTimeSpanMinutes(start)+" minutes");
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void getTopRowsFromDir(String inputDir, int n, String outputDir){
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir);
        String[] names= files.stream().map(File::getName).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e-> ScriptMethods.getTopRows(files.get(e), n, new File(outputDir, names[e])));
    }

    /**
     *
     * @param inputFile 有header
     * @param rate
     * @param subsetFile 有header
     */
    public static void getSubsetFromFile(String inputFile, double rate, String subsetFile){
        long start=System.nanoTime();
        BufferedReader br=null;
        BufferedWriter bw = null;
        try {
            if (inputFile.endsWith("gz")){
                br=IOUtils.getTextGzipReader(inputFile);
                bw=IOUtils.getTextGzipWriter(subsetFile);
            }else {
                br=IOUtils.getTextReader(inputFile);
                bw=IOUtils.getTextWriter(subsetFile);
            }
            String header=br.readLine();
            bw.write(header);
            bw.newLine();
            String line;
            double r=-1d;
            int count=0;
            int total=0;
            while ((line=br.readLine())!=null){
                total++;
                r=Math.random();
                if (r>rate) continue;
                count++;
                bw.write(line);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println("samping "+NumberTool.parse(count)+"("+NumberTool.parse(total)+") row from "
                    +new File(inputFile).getName()+" into "+new File(subsetFile).getName()+" in "
                    +Benchmark.getTimeSpanMinutes(start)+" minutes");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void getSubsetFromDir(String inputDir, double rate, String outDir){
        long start=System.nanoTime();
        File[] input=IOUtils.listRecursiveFiles(new File(inputDir));
        Predicate<File> p=File::isHidden;
        String[] files= Arrays.stream(input).filter(p.negate()).map(File::getAbsolutePath).toArray(String[]::new);
        String[] filesName=Arrays.stream(input).filter(p.negate()).map(File::getName)
                .map(str->str.replaceAll("vcf", "subset.vcf")).toArray(String[]::new);
        IntStream.range(0, files.length).forEach(e->ScriptMethods.getSubsetFromFile(files[e], rate, new File(outDir,
                filesName[e]).getAbsolutePath()));
        System.out.println(outDir+" subset is completed in "+Benchmark.getTimeSpanHours(start)+" hours");
    }

    public static void getSubsetFromDir(String inputDir, String outDir){
        ScriptMethods.getSubsetFromDir(inputDir, 0.001, outDir);
    }

    public static void calculateLD(String ingputDir, int binWidth_kb, int threshForDistance_Mb, String outDir){
        File[] files= IOUtils.listRecursiveFiles(new File(ingputDir));
        Predicate<File> h= File::isHidden;
        File[] f=Arrays.stream(files).filter(h.negate()).toArray(File[]::new);
        String[] outNames=Arrays.stream(f).map(File::getName).toArray(String[]::new);
        IntStream.range(0, f.length).parallel().forEach(e->{
            calculateLD(f[e], binWidth_kb, threshForDistance_Mb, new File(outDir, outNames[e]));
        });
    }

    public static void calculateLD(File ingputFile, int binWidth_kb, int threshForDistance_Mb, File outFile){
        try (BufferedReader br = IOUtils.getTextReader(ingputFile.getAbsolutePath());
             BufferedWriter bw=IOUtils.getTextWriter(outFile.getAbsolutePath())) {
            String line;
            List<String> temp;
            int distance=Integer.MIN_VALUE;
            double r2=Double.MIN_VALUE;
            int thresh=binWidth_kb*1000;
            int limit=threshForDistance_Mb*1000000;
            int kb=binWidth_kb;
            DescriptiveStatistics r2Stats=new DescriptiveStatistics();
            br.readLine();
            bw.write("window_kb\tnumberInWindow\tmeanOfR2\n");
            StringBuilder sb;
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                distance=Integer.parseInt(temp.get(0));
                if (distance > limit) break;
                r2=Double.parseDouble(temp.get(1));
                if (distance > thresh){
                    sb=new StringBuilder();
                    sb.append(kb).append("\t").append(r2Stats.getN()).append("\t").append(r2Stats.getMean());
                    bw.write(sb.toString());
                    bw.newLine();
                    thresh+=binWidth_kb*1000;
                    kb+=binWidth_kb;
                    r2Stats=new DescriptiveStatistics();
                }
                r2Stats.addValue(r2);
            }
            sb=new StringBuilder();
            sb.append(kb).append("\t").append(r2Stats.getN()).append("\t").append(r2Stats.getMean());
            bw.write(sb.toString());
            bw.newLine();
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }

    }

    public static void getAncestral(String inputDir, String ancestralDir){
        List<File> files= IOUtils.getVisibleFileListInDir(inputDir);
        String[] outNames=files.stream().map(File::getName).map(str->str.replaceAll("\\.txt",".ancestral.txt"))
                .toArray(String[]::new);
        String[] subgenome=files.stream().map(File::getName).map(str->str.substring(0,1)).toArray(String[]::new);
        IntStream.range(0, files.size()).forEach(e->getAncestral(files.get(e),new File(ancestralDir, outNames[e]),
                subgenome[e]));
    }

    public static void getAncestral(File inputFile, File ancestralFile, String subgenome){
        String[] acgt={"A","C","G","T"};
        try (BufferedReader bufferedReader = IOTool.getReader(inputFile);
             BufferedWriter bufferedWriter=IOTool.getTextWriter(ancestralFile)) {
            String line;
            List<String> temp, secer, hv;
            int chrID, pos, indexOfSecer, indexOfHv;
            bufferedWriter.write("chr\tpos\tancestral\n");
            StringBuilder sb;
            while ((line=bufferedReader.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                chrID= RefV1Utils.getChrID(temp.get(0).substring(3,4)+subgenome, Integer.parseInt(temp.get(2)));
                pos=RefV1Utils.getPosOnChrID(temp.get(0).substring(3,4)+subgenome, Integer.parseInt(temp.get(2)));
                secer=PStringUtils.fastSplit(temp.get(4),",");
                hv=PStringUtils.fastSplit(temp.get(5), ",");
                indexOfSecer=secer.indexOf("1");
                indexOfHv=hv.indexOf("1");
                if (indexOfHv!=indexOfSecer) continue;
                sb=new StringBuilder();
                sb.append(chrID).append("\t").append(pos).append("\t").append(acgt[indexOfHv]);
                bufferedWriter.write(sb.toString());
                bufferedWriter.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void retainVmapIExonGeno(String pgfFile, String geneHCFile, String genoInputDir, String outDir){
        PGF pgf=new PGF(pgfFile);
        System.out.println(pgf.getGeneNumber());
        Set<String> hcGenes=RowTableTool.getColumnSet(geneHCFile,0);
        Predicate<PGF.Gene> hcGenePredict=gene -> hcGenes.contains(gene.getGeneName());
        pgf.removeIf(hcGenePredict.negate());
        System.out.println(pgf.getGeneNumber());
        List<File> genoFiles=IOUtils.getVisibleFileListInDir(genoInputDir);
        try {
            BufferedReader bufferedReader;
            BufferedWriter bufferedWriter;
            String line, outFileName, header;
            List<String> temp;
            int chrID, vcfPos, geneIndex;
            for (int i = 0; i < genoFiles.size(); i++) {
                outFileName=genoFiles.get(i).getName().replaceAll("geno","exon.geon");
                bufferedReader=IOTool.getReader(genoFiles.get(i));
                bufferedWriter=IOTool.getTextGzipWriter(new File(outDir, outFileName));
                header=bufferedReader.readLine();
                bufferedWriter.write(header);
                bufferedWriter.newLine();
                while ((line=bufferedReader.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    chrID= RefV1Utils.getChrID(temp.get(0), Integer.parseInt(temp.get(1)));
                    vcfPos=RefV1Utils.getPosOnChrID(temp.get(0), Integer.parseInt(temp.get(1)));
                    geneIndex=pgf.getGeneIndex(chrID, vcfPos);
                    if (geneIndex<0) continue;
                    if (!pgf.isWithinThisGeneExon(geneIndex, chrID, vcfPos)) continue;
                    bufferedWriter.write(line);
                    bufferedWriter.newLine();
                }
                bufferedReader.close();
                bufferedWriter.flush();
                bufferedWriter.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     *
     * @param vcfDir
     * @param vcfComplementBedDir
     * @param subgenome "AB" or "D"
     * @param groupFile
     * DomesticatedEmmer	CItr14822
     * DomesticatedEmmer	CItr14824
     * DomesticatedEmmer	CItr14916
     * FreeThreshTetraploid	CItr14892
     * FreeThreshTetraploid	CItr7798
     * @param group
     * @param shOutFile
     */
    public static void smc(String vcfDir, String vcfComplementBedDir, String outDir, String subgenome, String groupFile,
                           String distTaxonFile, String[] group, String shOutFile, String genomeFile, String logDir){
        List<File> vcfFiles=IOUtils.getFileListInDirEndsWith(vcfDir, "gz");
        List<File> vcfComplementBedFiles=IOUtils.getFileListInDirEndsWith(vcfComplementBedDir, "gz");
        Predicate<File> subgenomeP = null;
        TIntArrayList subgenomeChrs = null;
        if (subgenome.equals("AB")){
            TIntArrayList ab=new TIntArrayList(WheatLineage.ablineage());
            subgenomeP=f->ab.contains(StringTool.getNumFromString(f.getName()));
            subgenomeChrs=new TIntArrayList(WheatLineage.ablineage());
        }else if (subgenome.equals("D")){
            TIntArrayList d=new TIntArrayList(WheatLineage.valueOf("D").getChrID());
            subgenomeP=f->d.contains(StringTool.getNumFromString(f.getName()));
            subgenomeChrs=new TIntArrayList(WheatLineage.valueOf("D").getChrID());
        }else {
            System.out.println("error your parameter subgenome");
            System.exit(1);
        }
        List<File> subgenomeVcfFiles=vcfFiles.stream().filter(subgenomeP).collect(Collectors.toList());
        List<File> subgenomeVcfComplementBedFiles=vcfComplementBedFiles.stream().filter(subgenomeP).collect(Collectors.toList());
        Multimap<String, String> map= TreeMultimap.create();
        Map<String, String> popDistTaxon=new HashMap<>();
        Map<Integer, Integer> chrSizeMap=new HashMap<>();
        try (BufferedReader bufferedReader = IOTool.getReader(groupFile);
             BufferedReader brDistTaxon=IOTool.getReader(distTaxonFile);
             BufferedReader bufferedReader1=IOTool.getReader(genomeFile);
             BufferedWriter bufferedWriter =IOTool.getTextWriter(shOutFile)) {
            String line;
            List<String> temp;
            bufferedReader.readLine();
            while ((line=bufferedReader.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                map.put(temp.get(0), temp.get(1));
            }
            brDistTaxon.readLine();
            while ((line=brDistTaxon.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                popDistTaxon.put(temp.get(0), temp.get(1));
            }
            while ((line=bufferedReader1.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                chrSizeMap.put(Integer.parseInt(temp.get(0)), Integer.parseInt(temp.get(1)));
            }
            StringBuilder sb;
            List<String> individuals;
            String distTaxon;
            File outFile = null;
            for (int i = 0; i < subgenomeVcfFiles.size(); i++) {
                for (int j = 0; j < group.length; j++) {
                    individuals=new ArrayList<>(map.get(group[j]));
                    distTaxon=popDistTaxon.get(group[j]);
                    sb=new StringBuilder();
                    sb.append("nohup smc++ vcf2smc --cores 1 -m ").append(subgenomeVcfComplementBedFiles.get(i)).append(" ");
                    sb.append(subgenomeVcfFiles.get(i)).append(" -l ").append(chrSizeMap.get(subgenomeChrs.get(i))).append(" -d ");
                    sb.append(distTaxon).append(" ").append(distTaxon).append(" ");
                    outFile=new File(outDir, subgenomeVcfFiles.get(i).getName().substring(0,6)+"."+group[j]+".smc.gz");
                    sb.append(outFile.getAbsolutePath()).append(" ").append(subgenomeChrs.get(i)).append(" ");
                    sb.append(group[j]).append(":");
                    for (int k = 0; k < individuals.size(); k++) {
                        sb.append(individuals.get(k)).append(",");
                    }
                    sb.deleteCharAt(sb.length()-1);
                    sb.append(" >"+logDir+"/"+subgenomeVcfFiles.get(i).getName().substring(0,6)+"_"+group[j]+".log");
                    sb.append(" ").append("2>&1");
                    bufferedWriter.write(sb.toString());
                    bufferedWriter.newLine();
                }
            }
            bufferedWriter.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     *
     * @param vcfDir
     * @param vcfComplementBedDir
     * @param subgenome "AB" or "D"
     * @param groupFile
     * DomesticatedEmmer	CItr14822
     * DomesticatedEmmer	CItr14824
     * DomesticatedEmmer	CItr14916
     * FreeThreshTetraploid	CItr14892
     * FreeThreshTetraploid	CItr7798
     * @param
     * @param shOutDir
     */
    public static void smc_split(String vcfDir, String vcfComplementBedDir, String outDir, String subgenome,
                            String groupFile, String distTaxonFile, String groupA, String groupB,
                                 String shOutDir, String genomeFile, String logDir){
        List<File> vcfFiles=IOUtils.getFileListInDirEndsWith(vcfDir, "gz");
        List<File> vcfComplementBedFiles=IOUtils.getFileListInDirEndsWith(vcfComplementBedDir, "gz");
        Predicate<File> subgenomeP = null;
        TIntArrayList subgenomeChrs = null;
        if (subgenome.equals("AB")){
            TIntArrayList ab=new TIntArrayList(WheatLineage.ablineage());
            subgenomeP=f->ab.contains(StringTool.getNumFromString(f.getName()));
            subgenomeChrs=new TIntArrayList(WheatLineage.ablineage());
        }else if (subgenome.equals("D")){
            TIntArrayList d=new TIntArrayList(WheatLineage.valueOf("D").getChrID());
            subgenomeP=f->d.contains(StringTool.getNumFromString(f.getName()));
            subgenomeChrs=new TIntArrayList(WheatLineage.valueOf("D").getChrID());
        }else {
            System.out.println("error your parameter subgenome");
            System.exit(1);
        }
        List<File> subgenomeVcfFiles=vcfFiles.stream().filter(subgenomeP).collect(Collectors.toList());
        List<File> subgenomeVcfComplementBedFiles=vcfComplementBedFiles.stream().filter(subgenomeP).collect(Collectors.toList());
        Multimap<String, String> map= TreeMultimap.create();
        Map<String, String> popDistTaxon=new HashMap<>();
        Map<String, Double> depthDistTaxon=new HashMap<>();
        Map<Integer, Integer> chrSizeMap=new HashMap<>();
        try (BufferedReader bufferedReader = IOTool.getReader(groupFile);
             BufferedReader brDistTaxon=IOTool.getReader(distTaxonFile);
             BufferedReader bufferedReader1=IOTool.getReader(genomeFile);
             BufferedWriter bufferedWriter =IOTool.getTextWriter(new File(shOutDir, groupA+"."+groupB+".sh"))) {
            String line;
            List<String> temp;
            bufferedReader.readLine();
            while ((line=bufferedReader.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                map.put(temp.get(0), temp.get(1));
            }
            brDistTaxon.readLine();
            while ((line=brDistTaxon.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                popDistTaxon.put(temp.get(0), temp.get(1));
                depthDistTaxon.put(temp.get(1), Double.parseDouble(temp.get(2)));
            }
            while ((line=bufferedReader1.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                chrSizeMap.put(Integer.parseInt(temp.get(0)), Integer.parseInt(temp.get(1)));
            }
            StringBuilder sb;
            List<String> individualsA, individualsB;
            String distTaxonA, distTaxonB, distTaxon;
            File outFile = null;
            for (int i = 0; i < subgenomeVcfFiles.size(); i++) {
                individualsA=new ArrayList<>(map.get(groupA));
                individualsB=new ArrayList<>(map.get(groupB));
                distTaxonA=popDistTaxon.get(groupA);
                distTaxonB=popDistTaxon.get(groupB);
                double depthDistTaxonA=depthDistTaxon.get(distTaxonA);
                double depthDistTaxonB=depthDistTaxon.get(distTaxonB);
                distTaxon= depthDistTaxonA > depthDistTaxonB ? distTaxonA : distTaxonB;
                sb=new StringBuilder();
                sb.append("nohup smc++ vcf2smc --cores 1 -m ").append(subgenomeVcfComplementBedFiles.get(i)).append(" ");
                sb.append(subgenomeVcfFiles.get(i)).append(" -l ").append(chrSizeMap.get(subgenomeChrs.get(i))).append(" -d ");
                sb.append(distTaxon).append(" ").append(distTaxon).append(" ");
                outFile=new File(outDir, subgenomeVcfFiles.get(i).getName().substring(0,6)+"."+groupA+"."+groupB+".smc" + ".txt");
                sb.append(outFile.getAbsolutePath()).append(" ").append(subgenomeChrs.get(i)).append(" ");
                sb.append(groupA).append(":");
                for (int k = 0; k < individualsA.size(); k++) {
                    sb.append(individualsA.get(k)).append(",");
                }
                sb.deleteCharAt(sb.length()-1);
                sb.append(" ").append(groupB).append(":");
                for (int j = 0; j < individualsB.size(); j++) {
                    sb.append(individualsB.get(j)).append(",");
                }
                sb.deleteCharAt(sb.length()-1);
                sb.append(" >"+logDir+"/"+subgenomeVcfFiles.get(i).getName().substring(0,6)+"_"+groupA+"."+groupB+".log");
                sb.append(" ").append("2>&1");
                bufferedWriter.write(sb.toString());
                bufferedWriter.newLine();
            }
            bufferedWriter.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void bulidSMC(String vcfDir, String bedDir, String outDir, String groupFile, String distTaxon,
                                String outFileAB, String outFileD, String genomeFile, String logDir){
//        vcfDir="/data1/home/daxing/vmap2.1Data/002_vmap2.1RefToAncestral";
//        bedDir="/data1/home/daxing/vmap2.1Data/001_vmap2.1_complementChr";
//        outDir="/Users/xudaxing/Desktop/out";
//        groupFile="/Users/xudaxing/Desktop/groupSMC_R1.txt";
//        distTaxon="/Users/xudaxing/Desktop/popDistTaxonR1.txt";
//        outFileAB="/Users/xudaxing/Desktop/resAB.sh";
//        outFileD="/Users/xudaxing/Desktop/resD.sh";
//        genomeFile="/Users/xudaxing/Desktop/genome.txt";
        String[] group_AB={"WildEmmer","DomesticatedEmmer","FreeThreshTetraploid","Landrace", "Cultivar"};
        String[] group_D={"Ae.tauschii","Landrace", "Cultivar"};
//        logDir="/Users/xudaxing/Desktop/log";
        ScriptMethods.smc(vcfDir, bedDir, outDir,"AB", groupFile, distTaxon, group_AB, outFileAB, genomeFile, logDir);
        ScriptMethods.smc(vcfDir, bedDir, outDir,"D", groupFile, distTaxon, group_D, outFileD, genomeFile, logDir);
    }

    public static void bulidSMC_split(String vcfDir, String bedDir, String outDir, String groupFile, String distTaxon,
                                      String outDirAB, String outDirD, String genomeFile, String logDir){
//        vcfDir="/Users/xudaxing/Desktop/vcf";
//        bedDir="/Users/xudaxing/Desktop/bed";
//        outDir="/Users/xudaxing/Desktop/out";
//        groupFile="/Users/xudaxing/Desktop/groupSMC_R1.txt";
//        distTaxon="/Users/xudaxing/Desktop/popDistTaxonR1.txt";
//        outDirAB="/Users/xudaxing/Desktop/resAB";
//        outDirD="/Users/xudaxing/Desktop/resD";
//        genomeFile="/Users/xudaxing/Desktop/genome.txt";
        String[] group_AB={"WildEmmer","DomesticatedEmmer","FreeThreshTetraploid","Landrace", "Cultivar"};
        String[] group_D={"Ae.tauschii","Landrace", "Cultivar"};
//        logDir="/Users/xudaxing/Desktop/log";
        String groupA, groupB;
        Iterator<int[]> abIterator= CombinatoricsUtils.combinationsIterator(group_AB.length, 2);
        Iterator<int[]> dIterator= CombinatoricsUtils.combinationsIterator(group_D.length, 2);
        while (abIterator.hasNext()){
            int[] combination=abIterator.next();
            groupA=group_AB[combination[0]];
            groupB=group_AB[combination[1]];
            ScriptMethods.smc_split(vcfDir, bedDir, outDir, "AB", groupFile, distTaxon, groupA, groupB, outDirAB,
                    genomeFile, logDir);
            ScriptMethods.smc_split(vcfDir, bedDir, outDir, "AB", groupFile, distTaxon, groupB, groupA, outDirAB,
                    genomeFile, logDir);
        }
        while (dIterator.hasNext()){
            int[] combination=dIterator.next();
            groupA=group_D[combination[0]];
            groupB=group_D[combination[1]];
            ScriptMethods.smc_split(vcfDir, bedDir, outDir, "D", groupFile, distTaxon, groupA, groupB, outDirD,
                    genomeFile, logDir);
            ScriptMethods.smc_split(vcfDir, bedDir, outDir, "D", groupFile, distTaxon, groupB, groupA, outDirD,
                    genomeFile, logDir);
        }
    }

    /**
     * 拆分sh命令
     * @param inputSH
     * @param outDir
     * @param numThreads
     * @param numCommands
     */
    public static void splitSh(String inputSH, String outDir, int numThreads, int numCommands){
        String[] subDir={"scriptAll", "scriptOne"};
        File[] files=new File[subDir.length];
        for (int i = 0; i < subDir.length; i++) {
            files[i]=new File(outDir, subDir[i]);
            files[i].mkdir();
        }
        BufferedWriter[] bws=new BufferedWriter[numThreads];
        for (int i = 0; i < bws.length; i++) {
            bws[i]=IOTool.getTextWriter(new File(files[0], "a"+PStringUtils.getNDigitNumber(3,i)+".sh"));
        }
        try (BufferedReader br = IOTool.getReader(inputSH);
             BufferedWriter bw =IOTool.getTextWriter(new File(files[1], "oneScript.sh"))) {
            String line;
            int count=0;
            int index=0;
            while ((line=br.readLine())!=null){
                count++;
                if (count > numCommands){
                    bws[index].flush();
                    index++;
                    count=1;
                    bws[index].write(line);
                    bws[index].newLine();
                    continue;
                }
                bws[index].write(line);
                bws[index].newLine();
            }
            bws[index].flush();
            List<File> files1=IOUtils.getFileListInDirEndsWith(files[0].getAbsolutePath(), "sh");
            StringBuilder sb;
            for (int i = 0; i < files1.size(); i++) {
                sb=new StringBuilder();
                sb.append("nohup sh ").append(outDir).append("/").append("scriptAll/").append(files1.get(i).getName());
                sb.append(" >"+outDir+"/scriptOne/"+files1.get(i).getName()+".log 2>&1").append(" &");
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void mergeSMC(String inputDir, String outDir, int binNumAB, int binNumD){
        mergeSMC_AB(inputDir, outDir, binNumAB);
        mergeSMC_D(inputDir, outDir, binNumD);
    }

    public static void mergeSMC_D(String inputDir, String outDir, int binNum){
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir);
        TIntArrayList dindex=new TIntArrayList(WheatLineage.valueOf("D").getChrID());
        Predicate<File> d=f->dindex.contains(StringTool.getNumFromString(f.getName()));
        List<File> dFiles=files.stream().filter(d).collect(Collectors.toList());
        TIntArrayList[] num=new TIntArrayList[2];
        num[0]=new TIntArrayList();
        num[1]=new TIntArrayList();
        int count=0;
        int index=0;
        boolean flag=false;
        for (int i = 0; i < dFiles.size(); i++) {
            count++;
            if (count>binNum){
                flag=!flag;
                index=flag ? 1 : 0;
                count=1;
            }
            num[index].add(i);
        }
        IntStream.range(0, num[0].size()).forEach(e->mergeSMC(dFiles.get(num[0].get(e)), dFiles.get(num[1].get(e)), outDir));
    }

    public static void mergeSMC_AB(String inputDir, String outDir, int binNum){
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir);
        TIntArrayList dindex=new TIntArrayList(WheatLineage.valueOf("D").getChrID());
        Predicate<File> d=f->dindex.contains(StringTool.getNumFromString(f.getName()));
        List<File> dFiles=files.stream().filter(d.negate()).collect(Collectors.toList());
        TIntArrayList[] num=new TIntArrayList[2];
        num[0]=new TIntArrayList();
        num[1]=new TIntArrayList();
        int count=0;
        int index=0;
        boolean flag=false;
        for (int i = 0; i < dFiles.size(); i++) {
            count++;
            if (count>binNum){
                flag=!flag;
                index=flag ? 1 : 0;
                count=1;
            }
            num[index].add(i);
        }
        IntStream.range(0, num[0].size()).forEach(e->mergeSMC(dFiles.get(num[0].get(e)), dFiles.get(num[1].get(e)), outDir));
    }

    public static void mergeSMC(File file1, File file2, String outDir){
        int chrID=StringTool.getNumFromString(file1.getName().substring(0,6));
        String filename=file1.getName().substring(6);
        String chr=RefV1Utils.getChromosome(chrID, 1);
        try (BufferedWriter bw = IOTool.getTextWriter(new File(outDir, "chr" + chr + filename))) {
            long numLine1=Files.lines(file1.toPath()).count();
            int num1=0;
            BufferedReader br1=IOTool.getReader(file1);
            BufferedReader br2=IOTool.getReader(file2);
            String header=br1.readLine();
            num1++;
            bw.write(header);
            bw.newLine();
            String line;
            List<String> temp1 = null, temp2=null;
            while ((line=br1.readLine())!=null && num1 < numLine1-1){
                num1++;
                bw.write(line);
                bw.newLine();
            }
            temp1=PStringUtils.fastSplit(line, " ");
            br1.close();
            br2.readLine();
            line=br2.readLine();
            temp2=PStringUtils.fastSplit(line, " ");
            int count=0;
            if (temp1.subList(1,temp1.size()).equals(temp2.subList(1,temp2.size()))){
                count=Integer.parseInt(temp1.get(0))+Integer.parseInt(temp2.get(0));
                bw.write(count+" "+String.join(" ", temp1.subList(1,4)));
                bw.newLine();
            }else {
                bw.write(line);
                bw.newLine();
            }
            while ((line=br2.readLine())!=null){
                bw.write(line);
                bw.newLine();
            }
            br2.close();
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void vmap2_exon_ML_Parsimony(String exonSNP_annoFile, String outEst, String outpasimony){
        try (BufferedReader bufferedReader = IOTool.getReader(exonSNP_annoFile);
             BufferedWriter bufferedWriter= IOTool.getTextWriter(outEst);
             BufferedWriter bufferedWriter1 =IOTool.getTextWriter(outpasimony)) {
            bufferedReader.readLine();
            String line, estAncestral, pasimonyAncestral;
            List<String> temp;
            StringBuilder sb;
            bufferedWriter.write("Chr-Pos-Ancestral-est");
            bufferedWriter1.write("Chr-Pos-Ancestral-pasimony");
            bufferedWriter.newLine();
            bufferedWriter1.newLine();
            while ((line=bufferedReader.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                estAncestral=temp.get(22);
                pasimonyAncestral=temp.get(31);
                if (estAncestral.equals("NA") && pasimonyAncestral.equals("NA")) continue;
                sb=new StringBuilder();
                if (!estAncestral.equals("NA")){
                    sb.append(temp.get(1)).append("-").append(temp.get(2)).append("-").append(estAncestral);
                    bufferedWriter.write(sb.toString());
                    bufferedWriter.newLine();
                }
                if (!pasimonyAncestral.equals("NA")){
                    sb=new StringBuilder();
                    sb.append(temp.get(1)).append("-").append(temp.get(2)).append("-").append(pasimonyAncestral);
                    bufferedWriter1.write(sb.toString());
                    bufferedWriter1.newLine();
                }
            }
            bufferedWriter.flush();
            bufferedWriter1.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void split_est_to_RefChr(String inputDir, String outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir);
        BufferedWriter[] bws;
        Map<String, BufferedWriter> chrBRmaps;
        String[] abd={"A","B","D"};
        BufferedReader br;
        BufferedWriter bw;
        for (int i = 0; i < abd.length; i++) {
            chrBRmaps=new HashMap<>();
            bws=new BufferedWriter[7];
            for (int j = 0; j < bws.length; j++) {
                int chr=j+1;
                bws[j]=IOTool.getTextWriter(new File(outDir, "chr"+chr+abd[i]+".est.txt"));
                chrBRmaps.put(chr+abd[i], bws[j]);
            }
            br=IOTool.getReader(files.get(i));
            try {
                String line, refChr;
                List<String> temp;
                while ((line=br.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    refChr=temp.get(0).substring(3)+abd[i];
                    bw=chrBRmaps.get(refChr);
                    temp.set(0, refChr);
                    bw.write(String.join("\t", temp));
                    bw.newLine();
                }
                br.close();
                for(Map.Entry<String, BufferedWriter> entry: chrBRmaps.entrySet()){
                    entry.getValue().flush();
                    entry.getValue().close();
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    private static Table<String, Integer, String> getTable(String inputFile){
        Table<String, Integer, String> table = HashBasedTable.create();
        try (BufferedReader bufferedReader = IOTool.getReader(inputFile)) {
            String line;
            List<String> temp;
            int chr, pos, refPos;
            String refChr;
            bufferedReader.readLine();
            while ((line=bufferedReader.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                if (temp.get(2).equals(temp.get(3))) continue;
                chr=Integer.parseInt(temp.get(0));
                pos=Integer.parseInt(temp.get(1));
                refChr=RefV1Utils.getChromosome(chr, pos);
                refPos=RefV1Utils.getPosOnChromosome(chr, pos);
                table.put(refChr, refPos, temp.get(2)+"-"+temp.get(3));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return table;
    }

    public static void test(String estByRefChrDir, String inputFile, String outFileML, String outFile){
        Table<String, Integer, String> table=getTable(inputFile);
        List<File> files=IOUtils.getVisibleFileListInDir(estByRefChrDir);
        BufferedReader br;
        BufferedWriter bw=IOTool.getTextWriter(outFileML);
        BufferedWriter bw1=IOTool.getTextWriter(outFile);
        try {
            String line, chr;
            List<String> temp;
            int pos;
            StringBuilder sb;
            bw.write("Chr\tPos\tAncestral-est_Ancestral-pasimony");
            bw.newLine();
            for (int i = 0; i < files.size(); i++) {
                br=IOTool.getReader(files.get(i));
                br.readLine();
                while ((line=br.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    chr=temp.get(0);
                    pos=Integer.parseInt(temp.get(2));
                    if (!table.contains(chr, pos)) continue;
                    sb=new StringBuilder();
                    sb.append(chr).append("\t").append(pos).append("\t").append(table.get(chr, pos));
                    bw.write(sb.toString());
                    bw.newLine();
                    bw1.write(line);
                    bw1.newLine();
                }
                br.close();
            }
            bw.flush();
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     *
     * @param softPath
     * @param genoFileDir
     * @param outDir
     * @param popsFilePath
     * @param p2File
     * @param p3File
     * @param outSHFile
     * @param ifDsubgenome
     */
    public static void fdSH(String softPath, String genoFileDir, String outDir, String popsFilePath, String p2File,
                            String p3File, String outSHFile, boolean ifDsubgenome){
        List<File> genoFiles=IOUtils.getVisibleFileListInDir(genoFileDir);
        Predicate<File> d=file -> file.getName().substring(4,5).equals("D");
        Predicate<File> ab=file -> file.getName().substring(4,5).equals("A") || file.getName().substring(4,5).equals("B");
        List<File> dGenoFiles=genoFiles.stream().filter(d).collect(Collectors.toList());
        List<File> abGenoFiles=genoFiles.stream().filter(ab).collect(Collectors.toList());
        List<File> filteredGeno=ifDsubgenome ? dGenoFiles : abGenoFiles;
        List<String> p2List=getColumn(p2File, 0);
        List<String> p3List=getColumn(p3File, 0);
        StringBuilder sb=new StringBuilder();
        try (BufferedWriter bw = IOTool.getTextWriter(outSHFile)) {
            String genoFile;
            for (int i = 0; i < filteredGeno.size(); i++) {
                for (int j = 0; j < p2List.size(); j++) {
                    for (int k = 0; k < p3List.size(); k++) {
                        genoFile=filteredGeno.get(i).getName().replaceAll(".geno.gz$","");
                        sb.setLength(0);
                        sb.append("python ").append(softPath).append(" --windType sites -g ").append(filteredGeno.get(i).getAbsolutePath());
                        sb.append("  -f phased -o ").append(outDir).append("/").append(genoFile).append("_");
                        sb.append(p2List.get(j)).append("_").append(p3List.get(k));
                        sb.append(".csv -w 100 -m 3 --overlap 50 -P1 IndianDwarfWheat -P2 ").append(p2List.get(j));
                        sb.append(" -P3").append(p3List.get(k)).append(" -O ancestral -T 1 --popsFile ").append(popsFilePath);
                        sb.append(" --writeFailedWindows");
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static List<String> getColumn(String pFile, int columnIndex){
        List<String> lines=new ArrayList<>();
        try (BufferedReader br = IOTool.getReader(pFile)) {
            String line;
            List<String> temp;
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                lines.add(temp.get(columnIndex));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return lines;
    }
}
