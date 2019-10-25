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
        IntStream.range(0, files.length).parallel().forEach(e->ScriptMethods.getSubsetFromFile(files[e], rate, new File(outDir,
                filesName[e]).getAbsolutePath()));
        System.out.println(outDir+" subset is completed in "+Benchmark.getTimeSpanHours(start)+" hours");
    }

    public static void getSubsetFromDir(String inputDir, String outDir){
        ScriptMethods.getSubsetFromDir(inputDir, 0.001, outDir);
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

    public static void groups(String inputFile, String outputFile, int groupsNumber){
        try {
            BufferedReader br=IOUtils.getTextReader(inputFile);
            BufferedWriter bw =IOUtils.getTextWriter(outputFile);
            String line;
            List<String> lineList;
            TDoubleArrayList dafSyn=new TDoubleArrayList();
            TDoubleArrayList dafNon=new TDoubleArrayList();
            String type;
            double siftScore=-1;
            double daf=-1;
            br.readLine();
            while ((line=br.readLine())!=null){
                lineList=PStringUtils.fastSplit(line);
                if (lineList.get(0).equals("NA")) continue;
                if (lineList.get(1).equals("NA")) continue;
                if (lineList.get(0).equals("Variant_type")) continue;
                type=lineList.get(0);
                siftScore=Double.parseDouble(lineList.get(1));
                daf=Double.parseDouble(lineList.get(2));
                if (lineList.get(0).equals("NONSYNONYMOUS") && siftScore<0.05){
                    dafNon.add(daf);
                }else if (lineList.get(0).equals("SYNONYMOUS")){
                    dafSyn.add(daf);
                }
            }
            br.close();
            double binsize=1d/groupsNumber;
            double[] groups= DoubleStream.iterate(binsize, n->n+binsize).limit(groupsNumber).toArray();
            int[] countSyn=new int[groups.length];
            int[] countNon=new int[groups.length];
            int index=-1;
            for (int i = 0; i < dafSyn.size(); i++) {
                index=Arrays.binarySearch(groups, dafSyn.get(i));
                if (index>-1){
                    countSyn[index+1]++;
                }else {
                    index=-index-1;
                    countSyn[index]++;
                }
            }
            for (int i = 0; i < dafNon.size(); i++) {
                index=Arrays.binarySearch(groups, dafNon.get(i));
                if (index>-1){
                    countNon[index+1]++;
                }else {
                    index=-index-1;
                    countNon[index]++;
                }
            }
            double[] rateSyn=ArrayTool.getElementPercent(countSyn);
            double[] rateNon=ArrayTool.getElementPercent(countNon);
            bw.write("type"+"\t"+"group"+"\t"+"Count"+"\t"+"rate");
            bw.newLine();
            for (int i = 0; i < groups.length; i++) {
                bw.write("SYNONYMOUS"+"\t"+groups[i] +"\t"+countSyn[i]+"\t"+ rateSyn[i]);
                bw.newLine();
                bw.write("NONSYNONYMOUS"+"\t"+groups[i] +"\t"+countNon[i]+"\t"+ rateNon[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public static void groups(String inputFile, String outputFile, int groupsNumber, boolean retain36){
        try {
            BufferedReader br=IOUtils.getTextReader(inputFile);
            BufferedWriter bw =IOUtils.getTextWriter(outputFile);
            String line;
            List<String> lineList;
            TDoubleArrayList dafSyn=new TDoubleArrayList();
            TDoubleArrayList dafNon=new TDoubleArrayList();
            String type;
            double siftScore=-1;
            double daf=-1;
            br.readLine();
            while ((line=br.readLine())!=null){
                lineList=PStringUtils.fastSplit(line);
                if (lineList.get(0).equals("NA")) continue;
                if (lineList.get(1).equals("NA")) continue;
                if (lineList.get(0).equals("Variant_type")) continue;
                type=lineList.get(0);
                siftScore=Double.parseDouble(lineList.get(1));
                daf=Double.parseDouble(lineList.get(2));
                if (lineList.get(0).equals("NONSYNONYMOUS") && siftScore<0.05){
                    dafNon.add(daf);
                }else if (lineList.get(0).equals("SYNONYMOUS")){
                    dafSyn.add(daf);
                }
            }
            br.close();
            double binsize=1d/groupsNumber;
            double[] groups= DoubleStream.iterate(binsize, n->n+binsize).limit(groupsNumber).toArray();
            int[] countSyn=new int[groups.length];
            int[] countNon=new int[groups.length];
            int index=-1;
            for (int i = 0; i < dafSyn.size(); i++) {
                index=Arrays.binarySearch(groups, dafSyn.get(i));
                if (index>-1){
                    countSyn[index+1]++;
                }else {
                    index=-index-1;
                    countSyn[index]++;
                }
            }
            for (int i = 0; i < dafNon.size(); i++) {
                index=Arrays.binarySearch(groups, dafNon.get(i));
                if (index>-1){
                    countNon[index+1]++;
                }else {
                    index=-index-1;
                    countNon[index]++;
                }
            }
            double[] rateSyn=ArrayTool.getElementPercent(countSyn);
            double[] rateNon=ArrayTool.getElementPercent(countNon);
            bw.write("type"+"\t"+"group"+"\t"+"Count"+"\t"+"rate");
            bw.newLine();
            List<Integer> flag=IntStream.range(18, 83).boxed().collect(Collectors.toList());
            for (int i = 0; i < groups.length; i++) {
                if (flag.contains(i)) continue;
                bw.write("SYNONYMOUS"+"\t"+groups[i] +"\t"+countSyn[i]+"\t"+ rateSyn[i]);
                bw.newLine();
                bw.write("NONSYNONYMOUS"+"\t"+groups[i] +"\t"+countNon[i]+"\t"+ rateNon[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public static void getSubPoplationFromAnnotationDB(File subPopulationChrPosFile, File annotationFile, File outFile){
        try (BufferedReader br1=IOUtils.getTextReader(subPopulationChrPosFile.getAbsolutePath());
             BufferedReader br2 = IOUtils.getTextGzipReader(annotationFile.getAbsolutePath());
             BufferedWriter bw=IOUtils.getTextGzipWriter(outFile.getAbsolutePath())) {
            String line;
            List<String> linesList;
            TIntHashSet chrSet=new TIntHashSet();
            TIntArrayList posList=new TIntArrayList();
            while ((line=br1.readLine()).startsWith("##")){}
            while ((line = br1.readLine())!=null){
                linesList=PStringUtils.fastSplit(line);
                chrSet.add(Integer.parseInt(linesList.get(0)));
                posList.add(Integer.parseInt(linesList.get(1)));
            }
            if (chrSet.size()>1){
                System.out.println("please check your subPopulationChrPosFile: "+subPopulationChrPosFile.getName()+" , it has more than one chromosome, program exist");
                System.exit(1);
            }
            String header=br2.readLine();
            bw.write(header);
            bw.newLine();
            int chr=-1;
            int pos=-1;
            int index=Integer.MIN_VALUE;
            while ((line=br2.readLine())!= null){
                linesList=PStringUtils.fastSplit(line);
                chr=Integer.parseInt(linesList.get(0));
                pos=Integer.parseInt(linesList.get(1));
                if (chr!= chrSet.iterator().next()){
                    System.out.println("please check your "+subPopulationChrPosFile.getName()
                            +" and "+annotationFile.getName()+", their chromosomes are not equal, program exist");
                    System.exit(1);
                }
                index=posList.binarySearch(pos);
                if (index<0) continue;
                bw.write(line);
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void getSubPoplationFromAnnotationDB(String subPopulationChrPosFileDir, String annotationFileDir, String outFileDir){
        File[] file1=IOUtils.listRecursiveFiles(new File(subPopulationChrPosFileDir));
        File[] file2=IOUtils.listRecursiveFiles(new File(annotationFileDir));
        Predicate<File> p=File::isHidden;
        File[] f1= Arrays.stream(file1).filter(p.negate()).toArray(File[]::new);
        File[] f2=Arrays.stream(file2).filter(p.negate()).toArray(File[]::new);
        String[] outFileName=Arrays.stream(f2).map(File::getName).map(str->str.replaceAll(".txt.gz$", "abd.subpop.txt.gz")).toArray(String[]::new);
        IntStream.range(0, f1.length).forEach(e-> ScriptMethods.getSubPoplationFromAnnotationDB(f1[e], f2[e], new File(outFileDir, outFileName[e])));
    }

    public static void getSubPopVcfFromMergedVCF(File mergedVcf, File subPopVcf1, File subPopVcf2){
        try (BufferedReader br = IOUtils.getTextReader(mergedVcf.getAbsolutePath());
            BufferedWriter bw1=IOUtils.getTextWriter(subPopVcf1.getAbsolutePath());
            BufferedWriter bw2=IOUtils.getTextWriter(subPopVcf2.getAbsolutePath())) {
            String line;
            List<String> lineList;
            List<String> infoList;
            List<String> subPopList1;
            List<String> subPopList2;
            List<String> l1;
            List<String> l2;
            StringBuilder sb = new StringBuilder();
            while ((line=br.readLine()).startsWith("##")){
                sb.append(line);
                sb.append("\n");
            }
            bw1.write(sb.toString());
            bw2.write(sb.toString());
            bw1.newLine();
            bw2.newLine();
            lineList=PStringUtils.fastSplit(line);
            infoList=lineList.stream().limit(9).collect(Collectors.toList());
            subPopList1=lineList.stream().skip(9).limit(419).collect(Collectors.toList());
            subPopList2=lineList.stream().skip(9+419).collect(Collectors.toList());
            l1=new ArrayList<>();
            l2=new ArrayList<>();
            l1.addAll(infoList);
            l1.addAll(subPopList1);
            l2.addAll(infoList);
            l2.addAll(subPopList2);
            bw1.write(l1.stream().collect(Collectors.joining("\t")));
            bw1.newLine();
            bw2.write(l2.stream().collect(Collectors.joining("\t")));
            bw2.newLine();
            double maf=Double.MIN_VALUE;
            List<String> temp;
            String tem;
            while ((line=br.readLine())!=null){
                lineList=PStringUtils.fastSplit(line);
                infoList=lineList.stream().limit(9).collect(Collectors.toList());
                subPopList1=lineList.stream().skip(9).limit(419).collect(Collectors.toList());
                subPopList2=lineList.stream().skip(9+419).collect(Collectors.toList());
                l1=new ArrayList<>();
                l2=new ArrayList<>();
                l1.addAll(infoList);
                l1.addAll(subPopList1);
                l2.addAll(infoList);
                l2.addAll(subPopList2);
                maf=VCF.calculateMaf(l1.stream().collect(Collectors.joining("\t")));
                temp=PStringUtils.fastSplit(l1.get(7), ";");
                tem="MAF="+maf;
                temp.set(4, tem);
                l1.set(7, temp.stream().collect(Collectors.joining(";")));
                bw1.write(l1.stream().collect(Collectors.joining("\t")));
                bw1.newLine();
                l2.set(7, temp.stream().collect(Collectors.joining(";")));
                bw2.write(l2.stream().collect(Collectors.joining("\t")));
                bw2.newLine();
            }
            bw1.flush();
            bw2.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void selectChrPosWithDaf(File annotationDB, File chrposOutFile){
        try (BufferedReader br = IOUtils.getTextGzipReader(annotationDB.getAbsolutePath());
            BufferedWriter bw=IOUtils.getTextWriter(chrposOutFile.getAbsolutePath())) {
            String line;
            List<String> lineList;
            lineList=PStringUtils.fastSplit(br.readLine());
            StringBuilder sb=new StringBuilder();
            sb.append(lineList.get(0)).append("\t").append(lineList.get(1));
            bw.write(sb.toString());
            bw.newLine();
            while ((line=br.readLine())!=null){
                sb=new StringBuilder();
                lineList=PStringUtils.fastSplit(line);
                if (lineList.get(13).equals("NA")) continue;
                sb.append(lineList.get(0)).append("\t").append(lineList.get(1));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void selectChrPosWithDaf(String annotationDBDir, String chrposOutFileDir){
        File[] files1=IOUtils.listRecursiveFiles(new File(annotationDBDir));
        Predicate<File> p=File::isHidden;
        File[] f1= Arrays.stream(files1).filter(p.negate()).toArray(File[]::new);
        String[] outNames=Arrays.stream(f1).map(File::getName).map(str->str.replaceAll("txt.gz$", "chrpos")).toArray(String[]::new);
        IntStream.range(0, f1.length).forEach(e->ScriptMethods.selectChrPosWithDaf(f1[e], new File(chrposOutFileDir, outNames[e])));
    }

    public static void caculateMaf(File inputChrposFile, File vcfFile, File outChrPosMaf){
        try (BufferedReader br1 = IOUtils.getTextReader(inputChrposFile.getAbsolutePath());
             BufferedReader br2=IOUtils.getTextGzipReader(vcfFile.getAbsolutePath())) {
             BufferedWriter bw=IOUtils.getTextWriter(outChrPosMaf.getAbsolutePath());
            bw.write("Chr"+"\t"+"Pos"+"\t"+"Maf"+"\n");
            String line;
             List<String> lineList;
             TIntHashSet chrSet=new TIntHashSet();
             TIntArrayList posList=new TIntArrayList();
             br1.readLine();
             while ((line=br1.readLine())!=null){
                 lineList=PStringUtils.fastSplit(line);
                 chrSet.add(Integer.parseInt(lineList.get(0)));
                 posList.add(Integer.parseInt(lineList.get(1)));
             }
             if(chrSet.size()>1) {
                 System.out.println("file error, check your "+inputChrposFile.getName()+" program exit");
                 System.exit(1);
             }
             int chr=-1;
             int pos=-1;
             double daf=-1d;
             int index=Integer.MIN_VALUE;
             StringBuilder sb;
             while ((line=br2.readLine()).startsWith("##")){}
             while ((line=br2.readLine())!=null){
                 lineList=PStringUtils.fastSplit(line);
                 chr=Integer.parseInt(lineList.get(0));
                 pos=Integer.parseInt(lineList.get(1));
                 if (chr!=chrSet.iterator().next()) {
                     System.out.println("please check your "+inputChrposFile.getName()+" and "+vcfFile.getName()+" , they " +
                             "have different chromosomes, program exit");
                     System.exit(1);
                 }
                 index=posList.binarySearch(pos);
                 if (index<0) continue;
                 daf=VCF.calculateMaf(line);
                 sb=new StringBuilder();
                 sb.append(chr).append("\t").append(pos).append("\t").append(daf);
                 bw.write(sb.toString());
                 bw.newLine();
             }
             bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void caculateMaf(String inputChrposFileDir, String vcfFileDir, String outChrPosDafDir){
        File[] files1=IOUtils.listRecursiveFiles(new File(inputChrposFileDir));
        int[] d_lineage=WheatLineage.getWheatLineageOf(WheatLineage.D);
        List<String> f2=VCF.getAllVcfInputPath(vcfFileDir, d_lineage);
        Predicate<File> p=File::isHidden;
        File[] f1=Arrays.stream(files1).filter(p.negate()).toArray(File[]::new);
        String[] outName=f2.stream().map(File::new).map(File::getName).map(str->str.substring(0, 6)+".chrpos.daf.txt").toArray(String[]::new);
        IntStream.range(0,f1.length).parallel().forEach(e->ScriptMethods.caculateMaf(f1[e], new File(f2.get(e)), new File(outChrPosDafDir, outName[e])));
    }

    public static  void reaplaceDafInAnnotationDB(File subPopChrposMafFile, File subPopAnnotationDBFile, File replacedSubPopAnnotationDBFile){
        try (BufferedReader br1 = IOUtils.getTextReader(subPopChrposMafFile.getAbsolutePath());
             BufferedReader br2=IOUtils.getTextGzipReader(subPopAnnotationDBFile.getAbsolutePath());
             BufferedWriter bw=IOUtils.getTextWriter(replacedSubPopAnnotationDBFile.getAbsolutePath())) {
            String line;
            List<String> lineList;
            TIntHashSet chrSet=new TIntHashSet();
            TIntArrayList posList=new TIntArrayList();
            TDoubleArrayList mafList=new TDoubleArrayList();
            br1.readLine();
            while ((line=br1.readLine())!=null){
                lineList=PStringUtils.fastSplit(line);
                chrSet.add(Integer.parseInt(lineList.get(0)));
                posList.add(Integer.parseInt(lineList.get(1)));
                mafList.add(Double.parseDouble(lineList.get(2)));
            }
            if (chrSet.size()>1){
                System.out.println("please check your "+subPopChrposMafFile.getName()+" program exit");
                System.exit(1);
            }
            String header=br2.readLine();
            bw.write(header);
            bw.newLine();
            int chr=-1;
            int pos=-1;
            int index=Integer.MIN_VALUE;
            String majorAllele;
            String minorAllele;
            String ancestralAllele;
            double daf=-1;
            StringBuilder sb;
            while ((line=br2.readLine())!=null){
                lineList=PStringUtils.fastSplit(line);
                chr=Integer.parseInt(lineList.get(0));
                pos=Integer.parseInt(lineList.get(1));
                majorAllele=lineList.get(4);
                minorAllele=lineList.get(5);
                ancestralAllele=lineList.get(12);
                if (chr!=chrSet.iterator().next()) {
                    System.out.println("please check your input file "+subPopChrposMafFile.getName()+" and "
                            +subPopAnnotationDBFile.getName()+" program exit");
                    System.exit(1);
                }
                index=posList.binarySearch(pos);
                sb=new StringBuilder();
                if (index<0) continue;
                if (ancestralAllele=="NA") continue;
                if (majorAllele.equals(ancestralAllele)){
                    daf=NumberTool.format(mafList.get(index), 5);
                    sb.append(lineList.stream().limit(6).collect(Collectors.joining("\t"))).append("\t")
                            .append(NumberTool.format(mafList.get(index),5)).append("\t");
                    sb.append(lineList.stream().skip(7).limit(6).collect(Collectors.joining("\t"))).append("\t");
                    sb.append(daf);
                    bw.write(sb.toString());
                    bw.newLine();
                }else if (minorAllele.equals(ancestralAllele)){
                    daf=NumberTool.format(1-mafList.get(index), 5);
                    sb.append(lineList.stream().limit(6).collect(Collectors.joining("\t"))).append("\t")
                            .append(NumberTool.format(mafList.get(index),5)).append("\t");
                    sb.append(lineList.stream().skip(7).limit(6).collect(Collectors.joining("\t"))).append("\t");
                    sb.append(daf);
                    bw.write(sb.toString());
                    bw.newLine();
                }else {
                    sb.append(lineList.stream().limit(6).collect(Collectors.joining("\t"))).append("\t")
                            .append(mafList.get(index)).append("\t");
                    sb.append(lineList.stream().skip(7).limit(6).collect(Collectors.joining("\t"))).append("\t");
                    sb.append("NA");
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static  void reaplaceDafInAnnotationDB(String subPopChrposMafFileDir, String subPopAnnotationDBFileDir, String replacedSubPopAnnotationDBFileDir){
        File[] files1=IOUtils.listRecursiveFiles(new File(subPopChrposMafFileDir));
        File[] files2=IOUtils.listRecursiveFiles(new File(subPopAnnotationDBFileDir));
        Predicate<File> p=File::isHidden;
        File[] f1=Arrays.stream(files1).filter(p.negate()).toArray(File[]::new);
        File[] f2=Arrays.stream(files2).filter(p.negate()).toArray(File[]::new);
        String[] outNames=Arrays.stream(f2).map(File::getName).toArray(String[]::new);
        IntStream.range(0,f2.length).parallel().forEach(e->
                ScriptMethods.reaplaceDafInAnnotationDB(f1[e], f2[e], new File(replacedSubPopAnnotationDBFileDir, outNames[e])));
    }

    public static void ww(String input, String out){
        try (BufferedReader br = IOUtils.getTextReader(input);
             BufferedWriter bw=IOUtils.getTextWriter(out)) {
            double maf=-1d;
            bw.write("maf"+"\n");
            String line;
            while ((line=br.readLine()).startsWith("##")){}
            while ((line=br.readLine())!=null){
                maf=VCF.calculateMaf(line);
                if (maf==0)continue;
                if (maf==-1) continue;
                bw.write(String.valueOf(maf));
                bw.newLine();
            }
        }catch (Exception e){
            e.printStackTrace();
        }
    }
}
