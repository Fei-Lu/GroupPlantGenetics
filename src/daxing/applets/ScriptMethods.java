package daxing.applets;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import daxing.common.*;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
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
                           String[] group, String shOutFile, String genomeFile){
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
        Multimap<String, String> map= ArrayListMultimap.create();
        Map<Integer, Integer> chrSizeMap=new HashMap<>();
        try (BufferedReader bufferedReader = IOTool.getReader(groupFile);
             BufferedReader bufferedReader1=IOTool.getReader(genomeFile);
             BufferedWriter bufferedWriter =IOTool.getTextWriter(shOutFile)) {
            String line;
            List<String> temp;
            while ((line=bufferedReader.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                map.put(temp.get(0), temp.get(1));
            }
            while ((line=bufferedReader1.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                chrSizeMap.put(Integer.parseInt(temp.get(0)), Integer.parseInt(temp.get(1)));
            }
            StringBuilder sb;
            List<String> individuals;
            File outFile = null;
            String logDir="/data4/home/aoyue/vmap2/daxing/analysis/016_demographic/log";
            for (int i = 0; i < subgenomeVcfFiles.size(); i++) {
                for (int j = 0; j < group.length; j++) {
                    individuals=new ArrayList<>(map.get(group[j]));
                    Collections.sort(individuals);
                    for (int k = 0; k < individuals.size(); k++) {
                        sb=new StringBuilder();
                        sb.append("nohup smc++ vcf2smc -m ").append(subgenomeVcfComplementBedFiles.get(i)).append(" ");
                        sb.append(subgenomeVcfFiles.get(i)).append(" -l ").append(chrSizeMap.get(subgenomeChrs.get(i))).append(" -d ");
                        sb.append(individuals.get(k)).append(" ").append(individuals.get(k)).append(" ");
                        outFile=new File(outDir, subgenomeVcfFiles.get(i).getName().substring(0,6)+"."+group[j]+"."+individuals.get(k)+".smc.txt");
                        sb.append(outFile.getAbsolutePath()).append(" ").append(subgenomeChrs.get(i)).append(" ");
                        sb.append(group[j]).append(":");
                        sb.append(String.join(",", individuals));
                        sb.append(" >"+logDir+"/"+subgenomeVcfFiles.get(i).getName().substring(0,6)+"_"+group[j]+"_"+individuals.get(k)+".log");
                        sb.append(" ").append("2>&1");
                        bufferedWriter.write(sb.toString());
                        bufferedWriter.newLine();
                    }
                }
            }
            bufferedWriter.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void bulidSMC(){
        String vcfDir="/Users/xudaxing/Desktop/vcf";
        String bedDir="/Users/xudaxing/Desktop/bed";
        String outDir="/Users/xudaxing/Desktop/ww";
        String subgenome="AB";
        String groupFile="/Users/xudaxing/Desktop/Ne.txt";
        String outFile="/Users/xudaxing/Desktop/resAB.sh";
        String genomeFile="/Users/xudaxing/Desktop/genome.txt";
        String[] group_AB={"WildEmmer_N","DomesticatedEmmer","FreeThreshTetraploid","Landrace_Europe", "Cultivar"};
        String[] group_D={"Ae.tauschii"};
        if (subgenome.equals("AB")){
            ScriptMethods.smc(vcfDir, bedDir, outDir,subgenome, groupFile, group_AB, outFile, genomeFile);
        }else if (subgenome.equals("D")){
            ScriptMethods.smc(vcfDir, bedDir, outDir,subgenome, groupFile, group_D, outFile, genomeFile);
        }else {
            System.out.println("error, check your subgenome parameter");
            System.exit(1);
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
                sb.append("nohup sh ../").append("scriptAll/").append(files1.get(i).getName()).append(" &");
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
