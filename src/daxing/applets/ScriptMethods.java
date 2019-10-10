package daxing.applets;

import daxing.common.*;
import daxing.filterSNP.Cells;
import daxing.filterSNP.DepthInfo;
import daxing.filterSNP.Dot;
import format.position.ChrPos;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import utils.Benchmark;
import utils.IOUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

public class ScriptMethods {

    /**
     * lastz /data1/home/daxing/reference/triticum_aestivum/bychr/chr005.fa /data1/home/daxing/reference/triticum_urartu/bychr/chr1A.fa --notransition --step=20  --nogapped --ambiguous=iupac --format=maf > /data1/home/daxing/ancestralAllele/MAF/ta_tu/wheatTauschii/Ta_chr005_vs_Au_chr1A.maf 2>/data1/home/daxing/ancestralAllele/MAF/ta_tu/wheatTauschiiLog/Ta_chr005_vs_Au_chr1ALog.txt &
     *
     * @param wheatInputDir
     * @param outgroupInputDir
     * @param outMAFDir
     * @param outSH
     */
    public static void getLastz(String wheatInputDir, String outgroupInputDir, String outMAFDir, String logFileDir, String outSH){
        File[] files1= IOUtils.listRecursiveFiles(new File(wheatInputDir));
        int[] d=ArrayTool.getWheatLineageOf(WheatLineage.D);
        List<Integer> l= Arrays.stream(d).boxed().collect(Collectors.toList());
        Predicate<File> dPredicate= e-> l.contains(StringTool.getNumFromString(e.getName()));
        File[] files2=IOUtils.listRecursiveFiles(new File(outgroupInputDir));
        Predicate<File> p= e->e.getName().endsWith("fa");
        File[] f1= Arrays.stream(files1).filter(p).filter(dPredicate).toArray(File[]::new);
        File[] f2=Arrays.stream(files2).filter(p).toArray(File[]::new);
        try(BufferedWriter bw=IOUtils.getTextWriter(outSH)){
            StringBuilder sb;
            sb=new StringBuilder();
            for (int i = 0; i < f1.length; i++) {
                for (int j = 0; j < f2.length; j++) {
                    sb.append("lastz ").append(f1[i].getAbsolutePath()).append(" ").append(f2[j].getAbsolutePath())
                            .append(" --notransition --step=20  --nogapped --ambiguous=iupac --format=maf > ")
                            .append(outMAFDir).append("/Ta_").append(f1[i].getName(), 0, 6).append("_vs_Au_")
                            .append(f2[j].getName(), 0, 5).append(".maf 2>").append(logFileDir+"/Ta_")
                            .append(f1[i].getName(), 0, 6).append("_vs_Au_").append(f2[j].getName(), 0, 5)
                            .append("Log.txt &");
                    bw.write(sb.toString());
                    bw.newLine();
                    sb=new StringBuilder();
                }
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

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
            if (inputFile.getName().endsWith(".gz")){
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
        File[] files=IOUtils.listRecursiveFiles(new File(inputDir));
        String[] names= Arrays.stream(files).map(File::getName).toArray(String[]::new);
        IntStream.range(0, files.length).parallel().forEach(e->{
            ScriptMethods.getTopRows(files[e], n, new File(outputDir, names[e]));
        });
    }

    // chr pos depth sd输入文件，输出目录
    public static void getCellDensity(String inputFile, String outputDir){
        DepthInfo depthInfo=new DepthInfo(inputFile);
        List<Dot> dotList=depthInfo.getDotList();
        Cells cells=new Cells(30F, 20F, 100);
        int[] indexOfDepthSD;
        for (int i = 0; i < dotList.size(); i++) {
            indexOfDepthSD= cells.binarySearch(dotList.get(i));
            cells.getCell(indexOfDepthSD).add(dotList.get(i));
        }
        cells.write(outputDir);
    }

    // 根据输出的postition目录，筛选累计密度达到指定百分比的格子所对应的ChrPos信息
    public static void getHighCumulativePos(String inputDirOfPosition, String outputFile, int num){
        File[] files=IOUtils.listRecursiveFiles(new File(inputDirOfPosition));
        List<String> chrPosList=new ArrayList<>();
        BufferedReader[] brs=new BufferedReader[num+1];
        try{
            for (int i = 0; i < brs.length; i++) {
                brs[i]=IOUtils.getTextReader(files[i].getAbsolutePath());
                String line;
                brs[i].readLine();
                List<String> lines=new ArrayList<>();
                while ((line=brs[i].readLine())!=null){
                    lines.add(line);
                }
                System.out.println(i+"\t"+lines.size());
                chrPosList.addAll(lines);
                brs[i].close();
            }
        }catch (Exception e){
            e.printStackTrace();
        }

        BufferedWriter bw=IOUtils.getTextWriter(outputFile);
//        int[] randoms= ArrayTool.getRandomNonrepetitionArray(5000, 0, chrPosList.size());
//        Arrays.sort(randoms);
        try {
            bw.write("CHR"+"\t"+"POS"+"\t"+"AverageDepth"+"\t"+"SD"+"\n");
//            int index;
            for (int i = 0; i < chrPosList.size(); i++) {
//                index=Arrays.binarySearch(randoms, i);
//                if (index<0) continue;
                bw.write(chrPosList.get(i));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void getSubsetForGraph(String inputFile, String outputFile, int num){
        try(BufferedReader br=IOUtils.getTextReader(inputFile);
            BufferedWriter bw=IOUtils.getTextWriter(outputFile)){
            String line;
            List<String> lineList=new ArrayList<>();
            String header=br.readLine();
            bw.write(header);
            bw.newLine();
            while ((line=br.readLine())!=null){
                lineList.add(line);
            }
            int[] randoms= ArrayTool.getRandomNonrepetitionArray(num, 0, lineList.size());
            Arrays.sort(randoms);
            int index;
            for (int i = 0; i < lineList.size(); i++) {
                index=Arrays.binarySearch(randoms, i);
                if (index<0) continue;
                bw.write(lineList.get(i));
                bw.newLine();
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }

    }

    public static void addDepthGroup(String databaseFile, String queryFile, int queryNum, String outFile){
        try(BufferedReader brQuery=IOUtils.getTextReader(queryFile);
            BufferedWriter bw=IOUtils.getTextWriter(outFile)) {
            ChrPos[] database;
            database=Files.newBufferedReader(Paths.get(databaseFile)).lines().skip(1).map(PStringUtils::fastSplit).map(l->new ChrPos(Short.parseShort(l.get(0)), Integer.parseInt(l.get(1)))).sorted().toArray(ChrPos[]::new);
            String line;
            short chr;
            int pos;
            List<String> lineList;
            List<String> lines=new ArrayList<>();
            String header=brQuery.readLine();
            bw.write(header+"\t"+"Group");
            bw.newLine();
            StringBuilder sb=new StringBuilder();
            int[] chrs= Arrays.stream(database).map(ChrPos::getChromosome)
                    .mapToInt(Integer::valueOf).toArray();
            int[] poss=Arrays.stream(database).map(ChrPos::getPosition)
                    .mapToInt(Integer::intValue).toArray();
            int indexChr;
            int indexPOS;
            while ((line=brQuery.readLine())!=null){
                lineList=PStringUtils.fastSplit(line);
                chr=Short.parseShort(lineList.get(0));
                pos=Integer.parseInt(lineList.get(1));
                indexChr= Arrays.binarySearch(chrs, chr);
                indexPOS=Arrays.binarySearch(poss, pos);
                if (indexChr>0 && indexPOS>0){
                    sb.append(line).append("\t").append(queryNum);
                    lines.add(sb.toString());
                }else {
                    sb.append(line).append("\t").append(0);
                    lines.add(sb.toString());
                }
                sb=new StringBuilder();
            }
            brQuery.close();
            for (int i = 0; i < lines.size(); i++) {
                bw.write(lines.get(i));
                bw.newLine();
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void getHapPos(File posAlleleDir, File hapPosDir){
        File[] files=IOUtils.listRecursiveFiles(posAlleleDir);
        BufferedReader[] brs=new BufferedReader[files.length];
        BufferedWriter[] bws=new BufferedWriter[files.length];
        for (int i = 0; i < files.length; i++) {
            brs[i]=IOUtils.getTextGzipReader(files[i].getAbsolutePath());
            bws[i]=IOUtils.getTextGzipWriter(files[i].getAbsolutePath());
        }
        List<ChrPos> chrPosList;
        StringBuilder sb=new StringBuilder();
        try{
            for (int i = 0; i < brs.length; i++) {
                chrPosList=brs[i].lines().skip(1).map(PStringUtils::fastSplit)
                        .map(e->new ChrPos(Short.parseShort(e.get(0)),Integer.parseInt(e.get(1)))).collect(Collectors.toList());
                sb=new StringBuilder();
                for (int j = 0; j < chrPosList.size(); j++) {
                    sb.append(chrPosList.get(j).getChromosome()).append("\t").append(chrPosList.get(j).getPosition());
                    bws[i].write(sb.toString());
                    bws[i].newLine();
                }
                bws[i].flush();
                bws[i].close();
            }
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void caculateChr(String inputFile, String outFile){
        try(BufferedReader br=IOUtils.getTextReader(inputFile);
            BufferedWriter bw=IOUtils.getTextWriter(outFile)){
            String header, line;
            List<String> lineList;
            header=br.readLine();
            bw.write(header);
            bw.newLine();
            int[] chrID;
            int[] snpNum;
            int[] biallelicNum;
            int[] indelNum;
            int[] insertionNum;
            int[] delectionNum;
            while ((line=br.readLine())!=null){
                chrID=new int[2];
                snpNum=new int[2];
                biallelicNum=new int[2];
                indelNum=new int[2];
                insertionNum=new int[2];
                delectionNum=new int[2];
                StringBuilder sb=new StringBuilder();
                lineList=PStringUtils.fastSplit(line);
                chrID[0]=Integer.parseInt(lineList.get(0));
                snpNum[0]=Integer.parseInt(lineList.get(1));
                biallelicNum[0]=Integer.parseInt(lineList.get(2));
                indelNum[0]=Integer.parseInt(lineList.get(3));
                insertionNum[0]=Integer.parseInt(lineList.get(4));
                delectionNum[0]=Integer.parseInt(lineList.get(5));
                line=br.readLine();
                lineList=PStringUtils.fastSplit(line);
                chrID[1]=Integer.parseInt(lineList.get(0));
                snpNum[1]=Integer.parseInt(lineList.get(1));
                biallelicNum[1]=Integer.parseInt(lineList.get(2));
                indelNum[1]=Integer.parseInt(lineList.get(3));
                insertionNum[1]=Integer.parseInt(lineList.get(4));
                delectionNum[1]=Integer.parseInt(lineList.get(5));
                sb.append(VCF.chrToChrMap.get(chrID[0])).append("\t").append(snpNum[0]+snpNum[1]).append("\t").append(biallelicNum[0]+biallelicNum[1]).append("\t").append(indelNum[0]+indelNum[1]).append("\t").append(insertionNum[0]+insertionNum[1]).append("\t").append(delectionNum[0]+delectionNum[1]);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void checkOutgroupAllele(String inputFile, String outDir){
        long start=System.nanoTime();
        String line;
        List<String> lineList;
        BufferedReader br;
        String header;
        TIntHashSet chrSet=new TIntHashSet();
        TIntArrayList chrList=new TIntArrayList();
        TIntArrayList posList=new TIntArrayList();
        List<String> refList=new ArrayList<>();
        List<String> outgroupList=new ArrayList<>();
        try{
            br=IOUtils.getTextReader(inputFile);
            header=br.readLine();
            while ((line=br.readLine())!=null){
                lineList=PStringUtils.fastSplit(line);
                chrSet.add(Integer.parseInt(lineList.get(0)));
                chrList.add(Integer.parseInt(lineList.get(0)));
                posList.add(Integer.parseInt(lineList.get(1)));
                refList.add(lineList.get(2));
                outgroupList.add(lineList.get(3));
            }
            br.close();
            System.out.println(inputFile+" reading is completed");
            TIntArrayList chrs=new TIntArrayList(chrSet);
            int[] chrNum=chrs.toArray();
            Arrays.sort(chrNum);
            BufferedWriter[] bws=new BufferedWriter[chrNum.length];
            for (int i = 0; i < chrNum.length; i++) {
                bws[i]=IOUtils.getTextWriter(new File(outDir, chrNum[i]+".txt").getAbsolutePath());
            }
            for (int i = 0; i < chrNum.length; i++) {
                bws[i].write(header);
                bws[i].newLine();
            }
            int index=Integer.MIN_VALUE;
            StringBuilder sb=new StringBuilder();
            for (int i = 0; i < chrList.size(); i++) {
                index=Arrays.binarySearch(chrNum, chrList.get(i));
                sb.append(chrList.get(i)).append("\t").append(posList.get(i)).append("\t")
                        .append(refList.get(i)).append("\t").append(outgroupList.get(i));
                bws[index].write(sb.toString());
                bws[index].newLine();
                sb=new StringBuilder();
            }
            for (int i = 0; i < bws.length; i++) {
                bws[i].flush();
                bws[i].close();
            }
            System.out.println(inputFile+" completed in "+ Benchmark.getTimeSpanMinutes(start)+" minutes");
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void mergeOutgroupAllele(String inputDir1, String inputDir2, String outDir){
        File[] inputFile1=IOUtils.listRecursiveFiles(new File(inputDir1).getAbsoluteFile());
        File[] inputFile2=IOUtils.listRecursiveFiles(new File(inputDir2).getAbsoluteFile());
        Predicate<File> p=File::isHidden;
        List<String> fileName1=Arrays.stream(inputFile1).filter(p.negate()).map(File::getName).collect(Collectors.toList());
        List<String> fileName2=Arrays.stream(inputFile2).filter(p.negate()).map(File::getName).collect(Collectors.toList());
        fileName1.retainAll(fileName2);
        BufferedReader[]  brs1=new BufferedReader[fileName1.size()];
        BufferedReader[]  brs2=new BufferedReader[fileName1.size()];
        BufferedWriter[]  bws =new BufferedWriter[fileName1.size()];
        for (int i = 0; i < fileName1.size(); i++) {
            brs1[i]=IOUtils.getTextReader(new  File(inputDir1, fileName1.get(i)).getAbsolutePath());
            brs2[i]=IOUtils.getTextReader(new File(inputDir2, fileName1.get(i)).getAbsolutePath());
            bws[i]= IOUtils.getTextWriter(new File(outDir, fileName1.get(i)).getAbsolutePath());
        }
        try {
            String line, header;
            List<String> lineList;
            TIntHashSet chrs;
            TIntArrayList posList;
            List<String> refList;
            List<String> outgroupList;
            StringBuilder sb=new StringBuilder();
            for (int i = 0; i < fileName1.size(); i++) {
                chrs=new TIntHashSet();
                posList=new TIntArrayList();
                refList=new ArrayList<>();
                outgroupList=new ArrayList<>();
                header=brs1[i].readLine();
                while ((line=brs1[i].readLine())!=null){
                    lineList=PStringUtils.fastSplit(line);
                    chrs.add(Integer.parseInt(lineList.get(0)));
                    posList.add(Integer.parseInt(lineList.get(1)));
                    refList.add(lineList.get(2));
                    outgroupList.add(lineList.get(3));
                }
                brs1[i].close();
                brs2[i].readLine();
                while ((line=brs2[i].readLine())!=null){
                    lineList=PStringUtils.fastSplit(line);
                    chrs.add(Integer.parseInt(lineList.get(0)));
                    posList.add(Integer.parseInt(lineList.get(1)));
                    refList.add(lineList.get(2));
                    outgroupList.add(lineList.get(3));
                }
                brs2[i].close();
                if (chrs.size()>1){
                    System.out.println(new File(inputDir1, fileName1.get(i)).getAbsolutePath()+"\t"+new File(inputDir2, fileName1.get(i)).getAbsolutePath()+" file error");
                    System.exit(1);
                }
                bws[i].write(header);
                bws[i].newLine();
                for (int j = 0; j < posList.size(); j++) {
                    sb.append(chrs.iterator().next()).append("\t").append(posList.get(j)).append("\t").append(refList.get(j)).append("\t").append(outgroupList.get(j));
                    bws[i].write(sb.toString());
                    bws[i].newLine();
                    sb=new StringBuilder();
                }
                bws[i].flush();
                bws[i].close();
            }
        }catch (Exception e){
            e.printStackTrace();
        }

    }

    //Chr Pos 将vcf chr pos 转换为 ref chr pos
    public static void mergeVCFPosIntoRefPos(String inDir, String outDir, String chrConvertionRuleFile){
        ChrConvertionRule chrConvertionRule=new ChrConvertionRule(Paths.get(chrConvertionRuleFile));
        File[] input=IOUtils.listRecursiveFiles(new File(inDir));
        Predicate<File> p=File::isHidden;
        File[] files=Arrays.stream(input).filter(p.negate()).toArray(File[]::new);
        int[] chrID=IntStream.range(1, 43).toArray();
        String[] subnames=Arrays.stream(files).map(File::getName).map(str->str.substring(6)).toArray(String[]::new);
        IntStream.iterate(0, n->n+2).limit(21).forEach(e->{
            ScriptMethods.mergeVCFPosIntoRefPos(files[e], files[e+1], new File(outDir, VCF.chrToChrMap.get(chrID[e])+subnames[e]), chrConvertionRule);
        });
    }

    public static void mergeVCFPosIntoRefPos(File infile1, File infile2, File outfile, ChrConvertionRule chrConvertionRule){
        try{
            List<ChrPos> chrPosList1=IOUtils.getTextGzipReader(infile1.getAbsolutePath()).lines().skip(1).map(PStringUtils::fastSplit)
                    .map(l->new ChrPos(Short.parseShort(l.get(0)), Integer.parseInt(l.get(1)))).collect(Collectors.toList());
            List<ChrPos> chrPosList2=IOUtils.getTextGzipReader(infile2.getAbsolutePath()).lines().skip(1).map(PStringUtils::fastSplit)
                    .map(l->new ChrPos(Short.parseShort(l.get(0)), Integer.parseInt(l.get(1)))).collect(Collectors.toList());
            chrPosList1.addAll(chrPosList2);
            Collections.sort(chrPosList1);
            BufferedWriter bw=IOUtils.getTextGzipWriter(outfile.getAbsolutePath());
            bw.write("Chr\tPos\n");
            StringBuilder sb;
            ChrPos chrPos;
            String chr;
            int refPos;
            for (int i = 0; i < chrPosList1.size(); i++) {
                chrPos=chrPosList1.get(i);
                chr=chrConvertionRule.getRefChrFromVCFChr(chrPos.getChromosome());
                refPos=chrConvertionRule.getRefPosFromVCFChrPos(chrPos);
                sb=new StringBuilder();
                sb.append(chr).append("\t").append(refPos);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void addChrPosToGerpRes(File gerpResFile, ChrConvertionRule chrConvertionRule, String chrStr, File outFileDir){
        int[] chrID=chrConvertionRule.getVCFChrFromRefChr(chrStr);
        String[] outFileNames=new String[2];
        outFileNames[0]="chr"+PStringUtils.getNDigitNumber(3, chrID[0])+".wheat.gerp++";
        outFileNames[1]="chr"+PStringUtils.getNDigitNumber(3, chrID[1])+".wheat.gerp++";
        BufferedWriter[] bws=new BufferedWriter[2];
        bws[0]=IOUtils.getTextWriter(new File(outFileDir, outFileNames[0]).getAbsolutePath());
        bws[1]=IOUtils.getTextWriter(new File(outFileDir, outFileNames[1]).getAbsolutePath());
        try (BufferedReader br = IOUtils.getTextReader(gerpResFile.getAbsolutePath())) {
            long numOfLine=IOUtils.getTextReader(gerpResFile.getAbsolutePath()).lines().count();
            int chrSize=chrConvertionRule.getChrSize(chrStr);
            if (numOfLine!=chrSize){
                System.out.println("The number of rows in "+gerpResFile.getName()+" is "+numOfLine+"\t"+"The "+chrStr+" size is "+chrSize);
                return;
            }
            String header="Chr\tPos\tGerpNeutralRate\tGerpScore\n";
            bws[0].write(header);
            bws[1].write(header);
            String line;
            List<String> lineList;
            int index=1;
            StringBuilder sb;
            short chr;
            int pos;
            ChrPos chrPos;
            int chrInitialValue=chrConvertionRule.getVCFChrFromRefChrPos(chrStr, 1);
            for (int i = 0; i < 2; i++) {
                while ((line=br.readLine())!=null){
                    lineList=PStringUtils.fastSplit(line);
                    chrPos=chrConvertionRule.getVCFChrPosFromRefChrPos(chrStr, index);
                    chr=chrPos.getChromosome();
                    if (chrInitialValue!=chr){
                        bws[i].flush();
                        bws[i].close();
                        i=1;
                        chrInitialValue=chr;
                    }
                    pos=chrPos.getPosition();
                    sb=new StringBuilder();
                    sb.append(chr).append("\t").append(pos).append("\t").append(lineList.get(0)).append("\t").append(lineList.get(1));
                    bws[i].write(sb.toString());
                    bws[i].newLine();
                    index++;
                }
                bws[i].flush();
                bws[i].close();
            }
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void addChrPosToGerpRes(String gerpResFileInDir, ChrConvertionRule chrConvertionRule, String outDir){
        File[] input=IOUtils.listRecursiveFiles(new File(gerpResFileInDir));
        File[] f=IOUtils.listFilesEndsWith(input, "gerp++");
        Predicate<File> p=File::isHidden;
        File[] files=Arrays.stream(f).filter(p.negate()).limit(7).toArray(File[]::new);
        int[] num= Arrays.stream(files).map(File::getName).map(e->PStringUtils.fastSplit(e, ".")).map(e->e.get(1))
                .map(s->s.substring(3)).mapToInt(Integer::parseInt).toArray();
        String[] str=Arrays.stream(files).map(File::getName).map(e->PStringUtils.fastSplit(e, ".")).map(e->e.get(0))
                .map(s->s.substring(5)).toArray(String[]::new);
        String[] chrs=new String[num.length];
        for (int i = 0; i < num.length; i++) {
            chrs[i]=num[i]+str[i];
        }
        IntStream.range(0, files.length).parallel().forEach(e->{
            ScriptMethods.addChrPosToGerpRes(files[e], chrConvertionRule, chrs[e], new File(outDir));
        });

    }

    public static void addAncestralAlleleToDB(String dbInputDir, String outDir){

    }

    public static void addAncestralAlleleToDB(File dbInputFile, File ancestralAlleleFile, File chrStr, File outFile){
        try(BufferedReader brDB=IOUtils.getTextGzipReader(dbInputFile.getAbsolutePath());
            BufferedWriter bw=IOUtils.getTextGzipWriter(outFile.getAbsolutePath())){

            String header=brDB.readLine();
            bw.write(header);
            String line;
            String lineList;
            while ((line=brDB.readLine())!=null){

            }
        }catch (Exception e){
            e.printStackTrace();
        }
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
            List<String> lines=new ArrayList<>();
            while ((line=br.readLine())!=null){
                lines.add(line);
            }
            br.close();
            int subsetLineNum = (int) Math.round(lines.size()*rate);
            int[] randomLinesIndex=ArrayTool.getRandomNonrepetitionArray(subsetLineNum, 0, lines.size());
            Arrays.sort(randomLinesIndex);
            for (int i = 0; i < randomLinesIndex.length; i++) {
                bw.write(lines.get(randomLinesIndex[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println(subsetFile+" is completed in "+Benchmark.getTimeSpanSeconds(start)+"s");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void getSubsetFromDir(String inputDir, double rate, String outDir){
        long start=System.nanoTime();
        File[] input=IOUtils.listRecursiveFiles(new File(inputDir));
        Predicate<File> p=File::isHidden;
        String[] files= Arrays.stream(input).filter(p.negate()).map(File::getAbsolutePath).toArray(String[]::new);
        String[] filesName=Arrays.stream(input).filter(p.negate()).map(File::getName).toArray(String[]::new);
        IntStream.range(0, files.length).forEach(e->ScriptMethods.getSubsetFromFile(files[e], rate, new File(outDir, filesName[e]).getAbsolutePath()));
        System.out.println(outDir+" subset is completed in "+Benchmark.getTimeSpanMinutes(start)+" minutes");
    }

    public static void getSubsetFromDir(String inputDir, String outDir){
        ScriptMethods.getSubsetFromDir(inputDir, 0.01, outDir);
    }

    public static void caculateAncestralAllele(File ancestralAlleFile, File dbFile, File outFile){
        long start = System.nanoTime();
        try (BufferedReader br1 = IOUtils.getTextReader(ancestralAlleFile.getAbsolutePath());
             BufferedReader br2=IOUtils.getTextGzipReader(dbFile.getAbsolutePath());
             BufferedWriter bw =IOUtils.getTextGzipWriter(outFile.getAbsolutePath())) {
            TIntHashSet chrSet=new TIntHashSet();
            int chr = -1;
            int chrDB=-1;
            TIntArrayList posList=new TIntArrayList();
            List<String> ancestralAlleleList=new ArrayList<>();
            String temp;
            List<String> tempList;
            int i=0;
            br1.readLine();
            String header=br2.readLine();
            StringBuilder sb=new StringBuilder();
            sb.append(header).append("\t").append("AncestralAllele").append("\t").append("Daf");
            bw.write(sb.toString());
            bw.newLine();
            while ((temp=br1.readLine())!=null){
                i++;
                tempList = PStringUtils.fastSplit(temp);
                chrSet.add(Integer.parseInt(tempList.get(0)));
                posList.add(Integer.parseInt(tempList.get(1)));
                ancestralAlleleList.add(tempList.get(3));
            }
            if (chrSet.size()>1){
                System.out.println("ancestralAlleFile error, program exit");
                System.exit(1);
            }else {
                chr=chrSet.iterator().next();
            }
            int index=Integer.MIN_VALUE;
            double daf=-1d;
            double maf=-1d;
            String majorAllele=null;
            String minorAllele=null;
            double count=0;
            while ((temp = br2.readLine())!=null){
                tempList = PStringUtils.fastSplit(temp);
                chrDB=Integer.parseInt(tempList.get(0));
                majorAllele=tempList.get(4);
                minorAllele=tempList.get(5);
                maf=Double.parseDouble(tempList.get(6));
                if (chr!=chrDB){
                    System.out.println("chr of dbFile and ancestralAlleFile is different ");
                    System.exit(1);
                }
                index=posList.binarySearch(Integer.parseInt(tempList.get(1)));
                sb=new StringBuilder();
                if (index < 0){
                   sb.append(temp).append("\t").append("NA").append("\t").append("NA");
                   bw.write(sb.toString());
                   bw.newLine();
                }else {
                    if (ancestralAlleleList.get(index).equals(majorAllele)){
                        daf=maf;
                        sb=sb.append(temp).append("\t").append(ancestralAlleleList.get(index)).append("\t").append(daf);
                    }else if (ancestralAlleleList.get(index).equals(minorAllele)){
                        daf=1-maf;
                        sb=sb.append(temp).append("\t").append(ancestralAlleleList.get(index)).append("\t").append(daf);
                    }else {
                        count++;
                        sb=sb.append(temp).append("\t").append(ancestralAlleleList.get(index)).append("\t").append("NA");
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            System.out.print(ancestralAlleFile.getName()+" "+(int) count+"/"+ i +" ("+String.format("%.4f", (count/i)*100)+"%)"+" ancestral allele neither belong to major allele nor minor allele");
            System.out.println(", completed in "+Benchmark.getTimeSpanMinutes(start)+" minutes");
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void caculateAncestralAllele(String ancestralAlleFileDir, String dbFileDir, String outFileDir){
        File[] files1=IOUtils.listRecursiveFiles(new File(ancestralAlleFileDir));
        File[] files2=IOUtils.listRecursiveFiles(new File(dbFileDir));
        Predicate<File> p=File::isHidden;
        File[] f1=Arrays.stream(files1).filter(p.negate()).toArray(File[]::new);
        File[] f2=Arrays.stream(files2).filter(p.negate()).toArray(File[]::new);
        String[] outName= Arrays.stream(f2).map(File::getName).map(str->str.replaceAll("SIFT", "Anc")).toArray(String[]::new);
        int[] dIndex=WheatLineage.getWheatLineageOf(WheatLineage.D);
        Arrays.stream(dIndex).forEach(e->{
            ScriptMethods.caculateAncestralAllele(f1[e-1], f2[e-1], new File(outFileDir, outName[e-1]));
        });
    }

    public static void groups(String inputFile, String outputFile){
        try {
            BufferedReader br=IOUtils.getTextReader(inputFile);
            BufferedWriter bw =IOUtils.getTextWriter(outputFile);
            String line;
            List<String> lineList;
            TDoubleArrayList daf=new TDoubleArrayList();
            br.readLine();
            while ((line=br.readLine())!=null){
                lineList=PStringUtils.fastSplit(line);
                daf.add(Double.parseDouble(lineList.get(1)));
            }
            br.close();
            double binsize=1d/100;
            double[] groups= DoubleStream.iterate(binsize, n->n+binsize).limit(100).toArray();
            int[] count=new int[groups.length];
            int index=-1;
            for (int i = 0; i < daf.size(); i++) {
                index=Arrays.binarySearch(groups, daf.get(i));
                if (index>-1){
                    count[index+1]++;
                }else {
                    index=-index-1;
                    count[index]++;
                }
            }
            double[] rate=ArrayTool.getElementPercent(count);
            bw.write("group");
            bw.newLine();
            for (int i = 0; i < count.length; i++) {
                bw.write(String.valueOf(count[i])+"\t"+String.valueOf(rate[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

}
