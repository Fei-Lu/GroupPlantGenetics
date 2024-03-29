package daxing.applets;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Multimap;
import com.google.common.collect.Table;
import com.google.common.collect.TreeMultimap;
import daxing.common.chrrange.ChrRange;
import daxing.common.factors.WheatLineage;
import daxing.common.genotype.GenoGrid;
import daxing.common.table.RowTableTool;
import daxing.common.utiles.IOTool;
import daxing.common.utiles.StringTool;
import daxing.common.wheat.PGF;
import daxing.common.wheat.Triads;
import daxing.load.complementary.TriadsBlock;
import daxing.load.complementary.TriadsBlockUtils;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.CombinatoricsUtils;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.range.Range;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.NumberFormat;
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
     * @param rate such as 0.01
     * @param subsetFile 有header
     */
    public static void getSubsetFromFile(String inputFile, double rate, String subsetFile){
        long start=System.nanoTime();
        BufferedReader br;
        BufferedWriter bw;
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
            double r;
            int count=0;
            int total=0;
            NumberFormat numberFormat=NumberFormat.getInstance();
            numberFormat.setGroupingUsed(true);
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
            System.out.println("samping "+numberFormat.format(count)+"("+numberFormat.format(total)+") row from "
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
                .map(str->str.replaceAll(".txt.gz", ".subset.001.txt.gz")).toArray(String[]::new);
        IntStream.range(0, files.length).parallel().forEach(e->ScriptMethods.getSubsetFromFile(files[e], rate, new File(outDir,
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
        IntStream.range(0, f.length).parallel().forEach(e-> calculateLD(f[e], binWidth_kb, threshForDistance_Mb, new File(outDir, outNames[e])));
    }

    public static void calculateLD(File ingputFile, int binWidth_kb, int threshForDistance_Mb, File outFile){
        try (BufferedReader br = IOUtils.getTextReader(ingputFile.getAbsolutePath());
             BufferedWriter bw=IOUtils.getTextWriter(outFile.getAbsolutePath())) {
            String line;
            List<String> temp;
            int distance;
            double r2;
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
             BufferedWriter bufferedWriter=IOTool.getWriter(ancestralFile)) {
            String line;
            List<String> temp, secer, hv;
            int chrID, pos, indexOfSecer, indexOfHv;
            bufferedWriter.write("chr\tpos\tancestral\n");
            StringBuilder sb;
            while ((line=bufferedReader.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                chrID= RefV1Utils.getChrID(temp.get(0).charAt(3)+subgenome, Integer.parseInt(temp.get(2)));
                pos=RefV1Utils.getPosOnChrID(temp.get(0).charAt(3)+subgenome, Integer.parseInt(temp.get(2)));
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
        Set<String> hcGenes= RowTableTool.getColumnSet(geneHCFile,0);
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
            for (File genoFile : genoFiles) {
                outFileName = genoFile.getName().replaceAll("geno", "exon.geon");
                bufferedReader = IOTool.getReader(genoFile);
                bufferedWriter = IOTool.getWriter(new File(outDir, outFileName));
                header = bufferedReader.readLine();
                bufferedWriter.write(header);
                bufferedWriter.newLine();
                while ((line = bufferedReader.readLine()) != null) {
                    temp = PStringUtils.fastSplit(line);
                    chrID = RefV1Utils.getChrID(temp.get(0), Integer.parseInt(temp.get(1)));
                    vcfPos = RefV1Utils.getPosOnChrID(temp.get(0), Integer.parseInt(temp.get(1)));
                    geneIndex = pgf.getGeneIndex(chrID, vcfPos);
                    if (geneIndex < 0) continue;
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
     * @param vcfDir vcfDir
     * @param vcfComplementBedDir vcfComplementBedDir
     * @param outDir outDir
     * @param subgenome "AB" or "D"
     * @param groupFile groupFile
     * DomesticatedEmmer	CItr14822
     * DomesticatedEmmer	CItr14824
     * DomesticatedEmmer	CItr14916
     * FreeThreshTetraploid	CItr14892
     * FreeThreshTetraploid	CItr7798
     * @param distTaxonFile distTaxonFile
     * @param groupA groupA
     * @param groupB groupB
     * @param shOutDir shOutDir
     * @param genomeFile genomeFile
     * @param logDir logDir
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
             BufferedWriter bufferedWriter =IOTool.getWriter(new File(shOutDir, groupA+"."+groupB+".sh"))) {
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
            File outFile;
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
                for (String s : individualsA) {
                    sb.append(s).append(",");
                }
                sb.deleteCharAt(sb.length()-1);
                sb.append(" ").append(groupB).append(":");
                for (String s : individualsB) {
                    sb.append(s).append(",");
                }
                sb.deleteCharAt(sb.length()-1);
                sb.append(" >").append(logDir).append("/").append(subgenomeVcfFiles.get(i).getName(), 0, 6).append("_").append(groupA).append(".").append(groupB).append(".log");
                sb.append(" ").append("2>&1");
                bufferedWriter.write(sb.toString());
                bufferedWriter.newLine();
            }
            bufferedWriter.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
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
     * @param inputSH inputSH
     * @param outDir outDir
     * @param numThreads numThreads
     * @param numCommands numCommands
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
            bws[i]=IOTool.getWriter(new File(files[0], "a"+PStringUtils.getNDigitNumber(3,i)+".sh"));
        }
        try (BufferedReader br = IOTool.getReader(inputSH);
             BufferedWriter bw =IOTool.getWriter(new File(files[1], "oneScript.sh"))) {
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
            for (File file : files1) {
                sb = new StringBuilder();
                sb.append("nohup sh ").append(outDir).append("/").append("scriptAll/").append(file.getName());
                sb.append(" >").append(outDir).append("/scriptOne/").append(file.getName()).append(".log 2>&1").append(" &");
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void vmap2_exon_ML_Parsimony(String exonSNP_annoFile, String outEst, String outpasimony){
        try (BufferedReader bufferedReader = IOTool.getReader(exonSNP_annoFile);
             BufferedWriter bufferedWriter= IOTool.getWriter(outEst);
             BufferedWriter bufferedWriter1 =IOTool.getWriter(outpasimony)) {
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
                bws[j]=IOTool.getWriter(new File(outDir, "chr"+chr+abd[i]+".est.txt"));
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
        BufferedWriter bw=IOTool.getWriter(outFileML);
        BufferedWriter bw1=IOTool.getWriter(outFile);
        try {
            String line, chr;
            List<String> temp;
            int pos;
            StringBuilder sb;
            bw.write("Chr\tPos\tAncestral-est_Ancestral-pasimony");
            bw.newLine();
            for (File file : files) {
                br = IOTool.getReader(file);
                br.readLine();
                while ((line = br.readLine()) != null) {
                    temp = PStringUtils.fastSplit(line);
                    chr = temp.get(0);
                    pos = Integer.parseInt(temp.get(2));
                    if (!table.contains(chr, pos)) continue;
                    sb = new StringBuilder();
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
     * @param softPath softPath
     * @param genoFileDir genoFileDir
     * @param outDir outDir
     * @param popsFilePath popsFilePath
     * @param p2File p2File
     * @param p3File p3File
     * @param outSHFile outSHFile
     * @param ifDsubgenome ifDsubgenome
     */
    public static void fdSH(String softPath, String genoFileDir, String outDir, String popsFilePath, String p2File,
                            String p3File, String outSHFile, boolean ifDsubgenome){
        List<File> genoFiles=IOUtils.getVisibleFileListInDir(genoFileDir);
        Predicate<File> d=file -> file.getName().charAt(4) == 'D';
        Predicate<File> ab=file -> file.getName().charAt(4) == 'A' || file.getName().charAt(4) == 'B';
        List<File> dGenoFiles=genoFiles.stream().filter(d).collect(Collectors.toList());
        List<File> abGenoFiles=genoFiles.stream().filter(ab).collect(Collectors.toList());
        List<File> filteredGeno=ifDsubgenome ? dGenoFiles : abGenoFiles;
        List<String> p2List=getColumn(p2File, 0);
        List<String> p3List=getColumn(p3File, 0);
        StringBuilder sb=new StringBuilder();
        try (BufferedWriter bw = IOTool.getWriter(outSHFile)) {
            String genoFile;
            for (File file : filteredGeno) {
                for (String s : p2List) {
                    for (String value : p3List) {
                        genoFile = file.getName().replaceAll(".geno.gz$", "");
                        sb.setLength(0);
                        sb.append("python ").append(softPath).append(" --windType sites -g ").append(file.getAbsolutePath());
                        sb.append("  -f phased -o ").append(outDir).append("/").append(genoFile).append("_");
                        sb.append(s).append("_").append(value);
                        sb.append(".csv -w 100 -m 3 --overlap 50 -P1 IndianDwarfWheat -P2 ").append(s);
                        sb.append(" -P3").append(value).append(" -O ancestral -T 1 --popsFile ").append(popsFilePath);
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

    /**
     *
     * @param triadsBlockChrRange important gene
     * @param pgfFile pgfFile
     * @param geneListFile geneListFile
     * @param outFile add triads and gene name information
     */
    public static void geneList(String triadsBlockChrRange, String pgfFile, String geneListFile, String outFile){
        TriadsBlock[] triadsBlocksSortByA= TriadsBlockUtils.readFromTriadsBlockChrRange(triadsBlockChrRange);
        TriadsBlock[] triadsBlocksSortByB=TriadsBlockUtils.readFromTriadsBlockChrRange(triadsBlockChrRange);
        TriadsBlock[] triadsBlocksSortByD=TriadsBlockUtils.readFromTriadsBlockChrRange(triadsBlockChrRange);
        Comparator<TriadsBlock> comparatorA=Comparator.comparing(triadsBlock -> triadsBlock.getChrRanges()[0]);
        Comparator<TriadsBlock> comparatorB=Comparator.comparing(triadsBlock -> triadsBlock.getChrRanges()[1]);
        Comparator<TriadsBlock> comparatorD=Comparator.comparing(triadsBlock -> triadsBlock.getChrRanges()[2]);
        Arrays.sort(triadsBlocksSortByA, comparatorA);
        Arrays.sort(triadsBlocksSortByB, comparatorB);
        Arrays.sort(triadsBlocksSortByD, comparatorD);
        List<TriadsBlock[]> triadsBlocksList = new ArrayList<>();
        List<Comparator<TriadsBlock>> comparatorList=new ArrayList<>();
        triadsBlocksList.add(triadsBlocksSortByA);
        triadsBlocksList.add(triadsBlocksSortByB);
        triadsBlocksList.add(triadsBlocksSortByD);
        comparatorList.add(comparatorA);
        comparatorList.add(comparatorB);
        comparatorList.add(comparatorD);
        PGF pgf=new PGF(pgfFile);
        pgf.sortGeneByGeneRange();
        try (BufferedReader br = IOTool.getReader(geneListFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            String line, chr, len;
            int start, end, subIndex, triadsBlockHit, triadsBlockIndexUp, triadsBlockIndexDown, geneIndex;
            int posStartOnChrID, posEndOnChrID;
            List<String> temp, tem;
            line=br.readLine();
            bw.write(line);
            bw.newLine();
            ChrRange chrRange;
            String[] abd={"A","B","D"};
            TriadsBlock triadsBlock;
            TriadsBlock[] triadsBlockArray;
            Set<String> triadsIDSet;
            Set<String> geneNameSet;
            List<String> triadsIDList;
            List<String> geneNameList;
            while ((line=br.readLine())!=null){
                tem= PStringUtils.fastSplit(line);
                temp=tem.stream().map(String::trim).collect(Collectors.toList());
                chr=temp.get(2);
                start=Integer.parseInt(temp.get(3).trim());
                end=Integer.parseInt(temp.get(4).trim())+1;
                len=temp.get(5).trim();
                chrRange=new ChrRange(chr, start, end);
                if (chr.equalsIgnoreCase("Un")){
                    temp.set(3, String.valueOf(start));
                    temp.set(4, String.valueOf(end));
                    temp.set(5, len);
                    bw.write(String.join("\t", temp));
                    bw.newLine();
                    continue;
                }
                subIndex=Arrays.binarySearch(abd, chr.substring(1,2));
                triadsBlock=new TriadsBlock(chrRange);
                triadsBlockArray=triadsBlocksList.get(subIndex);
                triadsBlockHit=Arrays.binarySearch(triadsBlockArray, triadsBlock, comparatorList.get(subIndex));
                triadsBlockIndexUp= triadsBlockHit < 0 ? -triadsBlockHit-2 : triadsBlockHit;
                triadsBlockIndexDown= triadsBlockHit < 0 ? -triadsBlockHit-1 : triadsBlockHit+1;
                triadsIDSet=new HashSet<>();
                geneNameSet=new HashSet<>();
                for (int i = triadsBlockIndexUp; i > -1; i--) {
                    if (triadsBlockArray[i].getChrRanges()[subIndex].isOverlapped(chrRange)) break;
                    triadsIDSet.add(triadsBlockArray[i].getTriadsID());
                }
                for (int i = triadsBlockIndexDown; i < triadsBlockArray.length; i++) {
                    if (triadsBlockArray[i].getChrRanges()[subIndex].isOverlapped(chrRange)) break;
                    triadsIDSet.add(triadsBlockArray[i].getTriadsID());
                }
                posStartOnChrID= RefV1Utils.getPosOnChrID(chr, start);
                posEndOnChrID=RefV1Utils.getPosOnChrID(chr, end);
                for (int i = posStartOnChrID; i <posEndOnChrID; i++) {
                    geneIndex=pgf.getGeneIndex(chrRange.getChrID(), i);
                    if (geneIndex < 0) continue;
                    geneNameSet.add(pgf.getGene(geneIndex).getGeneName());
                }
                triadsIDList=new ArrayList<>(triadsIDSet);
                geneNameList=new ArrayList<>(geneNameSet);
                Collections.sort(triadsIDList);
                Collections.sort(geneNameList);
                temp.set(3, String.valueOf(start));
                temp.set(4, String.valueOf(end));
                temp.set(5, len);
                temp.set(8, String.join(",", triadsIDList));
                temp.set(9, String.join(",", geneNameList));
                bw.write(String.join("\t", temp));
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void geneTriadsList(String triadsFile, String pgfFile, String geneListFile, String outFile){
        Triads triads =new Triads(triadsFile);
        PGF pgf=new PGF(pgfFile);
        pgf.sortGeneByGeneRange();
        try (BufferedReader br = IOTool.getReader(geneListFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            String line, chr;
            int start, end;
            line=br.readLine();
            bw.write(line);
            bw.newLine();
            List<String> temp;
            ChrRange chrRange;
            int geneIndex;
            Set<String> geneNameSet;
            List<String> geneNameList;
            StringBuilder sb=new StringBuilder();
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                chr=temp.get(2);
                start=Integer.parseInt(temp.get(3).trim());
                end=Integer.parseInt(temp.get(4).trim())+1;
                chrRange=new ChrRange(chr, start, end);
                geneNameSet=new HashSet<>();
                for (int i = chrRange.getStartOnChrID(); i < chrRange.getEndOnChrID(); i++) {
                    geneIndex=pgf.getGeneIndex(chrRange.getChrID(), i);
                    if (geneIndex < 0) continue;
                    geneNameSet.add(pgf.getGene(geneIndex).getGeneName());
                }
                geneNameList=new ArrayList<>(geneNameSet);
                Collections.sort(geneNameList);

            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void reCollectedUncompletedScript(String shDir, String shLogDir, String outDir){
        List<File> shFiles = IOTool.getFileListInDirEndsWith(shDir, "sh");
        List<File> shLogFiles = IOTool.getFileListInDirEndsWith(shLogDir, "log");
        String[] outNames=shFiles.stream().map(File::getName).map(s -> s.replaceAll("sh","uncompleted_sh")).toArray(String[]::new);
        IntStream.range(0, shFiles.size()).parallel().forEach(e->{
            try (BufferedReader brSh = IOTool.getReader(shFiles.get(e));
                 BufferedReader brShLog=IOTool.getReader(shLogFiles.get(e));
                 BufferedWriter bw = IOTool.getWriter(new File(outDir, outNames[e]))) {
                String line;
                int cnt=0;
                int count=0;
                while((line= brShLog.readLine())!=null){
                    if (line.startsWith("Done")){
                        cnt++;
                    }
                }
                while ((line= brSh.readLine())!=null){
                    count++;
                    if (count > cnt){
                        bw.write(line);
                        bw.newLine();
                    }
                }
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }

    public static void searchRunnedCommand(String inputResDir, String runnedNameOutFile){
        List<File> files = IOTool.getFileListInDirEndsWith(inputResDir, "csv");
        Predicate<File> filePredicate = file -> file.length() > 200;
        List<String> runnedFiles=files.stream().filter(filePredicate).map(File::getName).collect(Collectors.toList());
        try (BufferedWriter bw = IOTool.getWriter(runnedNameOutFile)) {
            bw.write(String.join("\n", runnedFiles));
            bw.newLine();
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void filterVMAP2(String inputDir, String coordinateRangeFile, double rate, String outDir){
        List<File> fileList = IOTool.getFileListInDirEndsWith(inputDir, ".gz");
        List<ChrRange> chrRanges=extractWindow(new File(coordinateRangeFile));
        Collections.sort(chrRanges);
        List<ChrRange> startRangeList=filterStart(chrRanges);
        List<ChrRange> endRangeList=filterEnd(chrRanges);
        String[] outNames= fileList.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz", "_"+rate+".vcf.gz")).toArray(String[]::new);
        IntStream.range(0, fileList.size()).parallel().forEach(e->{
            try (BufferedReader br = IOTool.getReader(fileList.get(e));
                 BufferedWriter bw =IOTool.getWriter(new File(outDir, outNames[e]))) {
                String line;
                List<String> temp;
                while ((line=br.readLine()).startsWith("##")){
                    bw.write(line);
                    bw.newLine();
                }
                bw.write(line);
                bw.newLine();
                int pos, chrID, index_start, index_end;
                ChrRange chrRange;
                double r;
                while ((line=br.readLine())!=null){
                    temp = PStringUtils.fastSplit(line.substring(0,30));
                    chrID = Integer.parseInt(temp.get(0));
                    pos =Integer.parseInt(temp.get(1));
                    chrRange = ChrRange.changeToChrRange(new Range(chrID, pos, pos+1));
                    index_start = Collections.binarySearch(startRangeList, chrRange);
                    index_end = Collections.binarySearch(endRangeList, chrRange);
                    if (index_start >= 0 || index_end >=0){
                        bw.write(line);
                        bw.newLine();
                    }else {
                        r=Math.random();
                        if (r>rate) continue;
                        bw.write(line);
                        bw.newLine();
                    }
                }
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }

    public static List<ChrRange> filterStart(List<ChrRange> chrRangeList){
        List<ChrRange> chrRanges=new ArrayList<>();
        for (ChrRange chrRange : chrRangeList){
            chrRanges.add(new ChrRange(chrRange.getChr(), chrRange.getStart(), chrRange.getStart()+1));
        }
        return chrRanges;
    }

    public static List<ChrRange> filterEnd(List<ChrRange> chrRangeList){
        List<ChrRange> chrRanges=new ArrayList<>();
        for (ChrRange chrRange : chrRangeList){
            chrRanges.add(new ChrRange(chrRange.getChr(), chrRange.getEnd()-1, chrRange.getEnd()));
        }
        return chrRanges;
    }

    public static List<ChrRange> extractWindow(File fdres){
        List<ChrRange> chrRanges=new ArrayList<>();
        try (BufferedReader br = IOTool.getReader(fdres)) {
            String line;
            List<String> temp;
            ChrRange chrRange;
            br.readLine();
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                chrRange=new ChrRange(temp.get(0), Integer.parseInt(temp.get(1)), Integer.parseInt(temp.get(2))+1);
                chrRanges.add(chrRange);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return chrRanges;
    }

    public static void getMissingRateFrom() throws IOException {
        Integer a = 5;
        String inputDir="/Users/xudaxing/Desktop/MissingRate/001_vcf";
        String svFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/004_burdenPerIndividualOrPerGene/002_1B1R/sv_taxa_GERP1.5_PPH0.5.txt";
        String nonSVFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/004_burdenPerIndividualOrPerGene/002_1B1R/nonsv_taxa_GERP1.5_PPH0.5.txt";
        String outDir="/Users/xudaxing/Desktop/MissingRate/002_missing";
        List<File> fileList = IOTool.getFileListInDirEndsWith(inputDir, ".gz");
        String[] outNames= fileList.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz",".missingRate.txt.gz")).toArray(String[]::new);
        List<String> svTaxa= Files.readAllLines(Paths.get(svFile));
        List<String> nonsvTaxa=Files.readAllLines(Paths.get(nonSVFile));
        IntStream.range(0, fileList.size()).forEach(e->{
            GenoGrid genoGrid = new GenoGrid(fileList.get(e).getAbsolutePath(), GenoIOFormat.VCF_GZ);
            try (BufferedWriter bw = IOTool.getWriter(new File(outDir, outNames[e]))) {
                StringBuilder sb = new StringBuilder();
                bw.write("ChrID\tPos\tMissingRate1B1R\tMissingRateNon1B1R");
                bw.newLine();
                int[] svIndices=new int[svTaxa.size()];
                int[] nonsvIndices= new int[nonsvTaxa.size()];
                for (int i = 0; i < svTaxa.size(); i++) {
                    svIndices[i]=genoGrid.getTaxonIndex(svTaxa.get(i));
                }
                for (int i = 0; i < nonsvTaxa.size(); i++) {
                    nonsvIndices[i]=genoGrid.getTaxonIndex(nonsvTaxa.get(i));
                }
                GenoGrid genoGridSV = genoGrid.getSubsetGenotypeByTaxon(svIndices);
                GenoGrid genoGridNonSV = genoGrid.getSubsetGenotypeByTaxon(nonsvIndices);
                double missingNumSV, missingNumNonSV;
                double totalCountSV, totalCountNonSV;
                double missingSV, missingNonSV;
                NumberFormat numberFormat = NumberFormat.getInstance();
                numberFormat.setGroupingUsed(false);
                numberFormat.setMaximumFractionDigits(5);
                for (int i = 0; i < genoGrid.getSiteNumber(); i++) {
                    missingNumSV= genoGridSV.getMissingNumberBySite(i);
                    totalCountSV = genoGridSV.getTaxaNumber();
                    missingNumNonSV= genoGridNonSV.getMissingNumberBySite(i);
                    totalCountNonSV = genoGridNonSV.getTaxaNumber();
                    missingSV= missingNumSV/totalCountSV;
                    missingNonSV= missingNumNonSV/totalCountNonSV;
                    sb.setLength(0);
                    sb.append(genoGrid.getChromosome(i)).append("\t").append(genoGrid.getPosition(i)).append("\t");
                    sb.append(numberFormat.format(missingSV)).append("\t").append(numberFormat.format(missingNonSV));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }

    public static void transform(String inputFile, String outFile){
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            String line;
            List<String> temp;
            line=br.readLine();
            bw.write(line);
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            int chrID, start, end, refStart, refEnd;
            String refChr;
            while ((line=br.readLine())!=null){
                temp =PStringUtils.fastSplit(line);
                chrID = Integer.parseInt(temp.get(0));
                start = Integer.parseInt(temp.get(2));
                end = Integer.parseInt(temp.get(3));
                refStart = RefV1Utils.getPosOnChromosome(chrID, start);
                refEnd = RefV1Utils.getPosOnChromosome(chrID, end);
                refChr = RefV1Utils.getChromosome(chrID, start);
                sb.setLength(0);
                sb.append(refChr).append("\t").append(temp.get(1)).append("\t");
                sb.append(refStart).append("\t").append(refEnd).append("\t");
                sb.append(String.join("\t", temp.subList(4,temp.size())));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * ancestral for AB: H. vulgare and Ae. tauschii
     * ancestral for D: H. vulgare and Triticum urartu
     * @param inputDir from kang ancestral VCF
     * @param outDir
     */
    public static void extractVMap1_ancestral(String inputDir, String outDir){
        List<File> files = IOTool.getFileListInDirEndsWith(inputDir, ".vcf.gz");
        String[] outNames = files.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz",".ancestral.txt.gz")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->{
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw = IOTool.getWriter(new File(outDir, outNames[e]))) {
                String line, chr, genotype1, genotype2;
                List<String> temp, altAlleleList;
                int ancestralAlleleIndex;
                StringBuilder sb = new StringBuilder();
                sb.append("Chr\tPos\tAncestral");
                bw.write(sb.toString());
                bw.newLine();
                chr = files.get(e).getName().substring(0,2);
                while ((line=br.readLine()).startsWith("##")) continue;
                List<String> abChrs = WheatLineage.abLineage();
                int ancestralIndex=abChrs.contains(chr) ? 12 : 11;
                List<String> alleleList;
                while ((line=br.readLine())!=null){
                    temp =PStringUtils.fastSplit(line);
                    genotype1 = temp.get(ancestralIndex);
                    genotype2 = temp.get(15);
                    if (genotype1.equals(".") || genotype2.equals(".")) continue;
                    if (!genotype1.equals(genotype2)) continue;
                    altAlleleList= PStringUtils.fastSplit(temp.get(4),",");
                    alleleList = new ArrayList<>();
                    alleleList.add(temp.get(3));
                    alleleList.addAll(altAlleleList);
                    ancestralAlleleIndex = Integer.parseInt(genotype1);
                    sb.setLength(0);
                    sb.append(chr).append("\t").append(temp.get(1)).append("\t");
                    sb.append(alleleList.get(ancestralAlleleIndex));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }

    public static void splitLoterVCF(String v2Dir, String outDir){
        List<File> v2Files = IOTool.getFileListInDirEndsWith(v2Dir, ".vcf.gz");
        for (int i = 0; i < v2Files.size(); i++) {
            splitLoterVCF(v2Files.get(i), outDir);
        }
    }

    public static void splitLoterVCF(File v2File, String outDir){
        int variantsNum = getVariantsNjum(v2File);
        int flag = (variantsNum/10) + 1;
        BufferedWriter[] bws = new BufferedWriter[10];
        String indexID, outNamePrefix;
        for (int i = 0; i < bws.length; i++) {
            indexID = PStringUtils.getNDigitNumber(2, i);
            outNamePrefix = v2File.getName().replaceAll(".vcf.gz","");
            bws[i] = IOTool.getWriter(new File(outDir, outNamePrefix+"_"+indexID+".vcf.gz"));
        }
        try (BufferedReader br = IOTool.getReader(v2File)) {
            String line;
            List<String> headerList;
            headerList = new ArrayList<>();
            while ((line= br.readLine()).startsWith("##")){
                headerList.add(line);
            }
            headerList.add(line);
            for (int i = 0; i < bws.length; i++) {
                for (int j = 0; j < headerList.size(); j++) {
                    bws[i].write(headerList.get(j));
                    bws[i].newLine();
                }
            }
            int index = 0;
            int count = 0;
            int flagCum = flag;
            while ((line=br.readLine())!=null){
                count++;
                if (count>flagCum){
                    index++;
                    flagCum += flag;
                }
                bws[index].write(line);
                bws[index].newLine();
            }
            for (int i = 0; i < bws.length; i++) {
                bws[i].flush();
                bws[i].close();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static int getVariantsNjum(File v2File){
        int variatnsNum = 0;
        try (BufferedReader br = IOTool.getReader(v2File)) {
            String line;
            while ((line=br.readLine()).startsWith("##")){}
            while ((line=br.readLine())!=null){
                variatnsNum++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return variatnsNum;
    }

}
